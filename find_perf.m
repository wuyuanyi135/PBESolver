function [perf, directory] = find_perf(key, outputDir)
arguments
    key;
    outputDir = 'outputs';
end

if isnumeric(key)
    % if key is a non-positive number:
    % 0: current commit
    % -n: the previous nth commit.
    % if the perf for the commit is missing, an error is thrown.
    intKey = int32(key);
    if intKey > 0
        error('Key cannot be a positive number');
    end
    commitHash = get_hash_by_offset(-intKey);
    
elseif isstr(key)
    % if input is a string, it should either be
    %   'head' | 'head~N' (case insensitive)
    %   'latest' (last entry in the output directory)
    %   '@[tag]'
    %   '[substring of the hash]'
    if startsWith(key, 'head', 'IgnoreCase', true)
        commitHash = get_hash_by_offset(key);
    elseif strcmp(key, 'latest')
        S = dir(outputDir);
        S = S(3:end);
        S = S([S.isdir]);
        [~,idx] = sort([S.datenum], 'descend');
        S = S(idx);
        commitHash = S(1).name;
    elseif startsWith(key, '@')
        tag = key(2:end);
        commitHash = get_hash_by_tag_or_substring(tag);
    else
        % query by substring of the hash
        commitHash = get_hash_by_tag_or_substring(key);
    end
else
    error('Unknown key type');
end
recordDir = fullfile(outputDir, commitHash);
if ~exist(recordDir, 'dir')
    error('%s does not have test entry', commitHash);
end

matFile = fullfile(recordDir, 'perf.mat');
data = load(matFile);
perf = data.perf;
directory = recordDir;
end

function commitHash = get_hash_by_offset(offset)
if isinteger(offset)
    cmd = sprintf('git rev-parse HEAD~%d', offset);
else
    cmd = sprintf('git rev-parse %s', offset);
end
[status, commitHash] = system(cmd);
if status ~= 0
    error('Failed to query commit id with head offset');
end
commitHash = strtrim(commitHash);
end

function commitHash = get_hash_by_tag_or_substring(tag)
cmd = sprintf('git rev-list -n 1 %s', tag);
[status, commitHash] = system(cmd);
if status ~= 0
    if status == 128
        error('%s is ambigious id.', tag);
    end
    error('Failed to query commit id with tag: %d', status);
end
commitHash = strtrim(commitHash);    
end
