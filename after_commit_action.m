% After commit, create a report and the benchmark inside the output folder.
DESTINATION = fullfile(pwd, 'outputs');

% Retrieving git hash
[status, commitHash] = system('git rev-parse HEAD');
commitHash = strtrim(commitHash);

if status ~= 0
    error('Failed to retrieve git hash: %d', status)
end

cprintf('blue', 'Current commit: %s\n', commitHash);

destinationFolder = fullfile(DESTINATION, commitHash);

if exist(destinationFolder, 'dir')
    error('Report has been generated for commit %s', commitHash);
end

cprintf('blue', 'Creating directory: %s\n', destinationFolder);
mkdir(destinationFolder);

fprintf('Running unittest and generate report: %s\n', destinationFolder);

import matlab.unittest.TestRunner;
import matlab.unittest.TestSuite;
import matlab.unittest.plugins.TestReportPlugin;
runner = TestRunner.withNoPlugins;
reportFile = fullfile(destinationFolder, 'test_report.pdf');
plugin = TestReportPlugin.producingPDF(reportFile);
runner.addPlugin(plugin);
suite = matlab.unittest.TestSuite.fromFolder(pwd);
result = runner.run(suite);

cprintf('blue', 'Running Benchmarks.\n');
perf = runperf;
cprintf('blue', 'Writing Benchmark to "perf.mat".\n');
save(fullfile(destinationFolder, 'perf.mat'), 'perf');

cprintf('blue', 'Done.\n');



