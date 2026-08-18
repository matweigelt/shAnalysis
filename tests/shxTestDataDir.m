function d = shxTestDataDir()
%SHXTESTDATADIR Resolve the test-fixture folder for all test suites.
%
%   D = shxTestDataDir() returns, in order of precedence:
%     1. getenv('SHX_TESTDATA_FOLDER')          if set and a folder
%     2. getpref('shAnalysis','TestDataFolder') if set and a folder
%     3. fullfile(<tests dir>, 'test_data')     the in-repo default
%   The first two follow the SHX_SERIES_FOLDER / setpref convention of
%   setup_shAnalysis (env first, preference second), so a machine can
%   keep the fixtures outside the repository (e.g. on a data drive)
%   while CI keeps using the in-repo copy unchanged.
%
%   Outputs
%     d  (1 x :) char  absolute path of the fixture folder (existence
%        is NOT enforced here - suites assume/verify per test so a
%        missing folder produces their usual, readable diagnostics)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-18, 12:20 UTC.

d = getenv('SHX_TESTDATA_FOLDER');
if ~isempty(d) && isfolder(d)
    return
end
if ispref('shAnalysis', 'TestDataFolder')
    d = getpref('shAnalysis', 'TestDataFolder');
    if ~isempty(d) && isfolder(d)
        return
    end
end
d = fullfile(fileparts(mfilename('fullpath')), 'test_data');
end
