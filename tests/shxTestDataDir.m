function d = shxTestDataDir()
%SHXTESTDATADIR Resolve the test-fixture folder for all test suites.
%
%   D = shxTestDataDir() delegates to shLowLevel.testDataDir (v3.26.2)
%   so toolbox code, demos, and suites share ONE resolution order:
%   SHX_TESTDATA_FOLDER env var, TestDataFolder preference,
%   tests/test_data fallback. This thin wrapper stays because the
%   suites call it unqualified and were released this way in v3.25.0.
%
%   Outputs
%     d  (1 x :) char  fixture folder path
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-18, 13:55 UTC.

d = shLowLevel.testDataDir();
end
