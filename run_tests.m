function result = run_tests(report)
suite = matlab.unittest.TestSuite.fromFolder(pwd);

if report
    import matlab.unittest.TestRunner;
    import matlab.unittest.TestSuite;
    import matlab.unittest.plugins.TestReportPlugin;
    runner = TestRunner.withNoPlugins;
    reportFile = fullfile('outputs', ['test_report_' datestr(datetime, 'yyyy-mm-dd_HH_MM_SS')]);
    plugin = TestReportPlugin.producingHTML(reportFile);
    runner.addPlugin(plugin);
    result = runner.run(suite);
else
    suite = matlab.unittest.TestSuite.fromFolder(pwd);
    result = suite.run;
end
end