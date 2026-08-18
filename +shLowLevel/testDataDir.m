function d = testDataDir()
%TESTDATADIR Resolve the sample/fixture data folder (relocatable).
%
%   D = shLowLevel.testDataDir() returns, in order of precedence:
%     1. getenv('SHX_TESTDATA_FOLDER')          if set and a folder
%     2. getpref('shAnalysis','TestDataFolder') if set and a folder
%     3. <toolbox>/tests/test_data              the in-repo default
%   Same contract as tests/shxTestDataDir (which the test suites use);
%   this packaged twin exists because TOOLBOX code also consults the
%   shipped samples - readDDK's "DDK3" name form and the demos - and
%   must not depend on the tests folder being on the path. The v3.25
%   relocation missed exactly that: with the fixtures moved to a data
%   drive, readDDK("DDK3") failed on the acceptance machine while all
%   suite-internal paths worked (caught by the 2026-08-18 acceptance
%   run, testDataFolderAndDDKNames).
%
%   Outputs
%     d  (1 x :) char  folder path (existence NOT enforced for the
%        fallback - callers keep their own diagnostics)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-18, 13:55 UTC.

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
d = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'tests', 'test_data');
end
