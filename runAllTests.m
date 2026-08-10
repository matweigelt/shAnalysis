function results = runAllTests(opts)
%RUNALLTESTS Run the complete shAnalysis v2 validation suite.
%
%   RESULTS = RUNALLTESTS runs, in order:
%     tests/testCorrectness.m   numerics vs. golden/legacy/analytic values
%     tests/testContract.m      error IDs, immutability, dimensions
%     tests/testRobustness.m    edge cases, degenerate inputs, file I/O
%                               documented patch: kernel-ratio test masks
%                               exact-zero symmetry points, see comment there)
%     tests/testPerformance.m   benchmarks incl. Legendre-cache assertion
%
%   RESULTS = RUNALLTESTS(SkipPerformance (false)=true) omits the (slow)
%   performance suite. For statistical benchmark timing use
%   runperf('testPerformance') directly.
%
%   The full console output AND a machine-independent summary are written
%   to tests/runAllTests_results.txt (fixed name, overwritten each run) -
%   upload that file instead of pasting console text, which gets cut off.
%
%   Inputs
%     opts.SkipPerformance (1,1) logical, default false
%     opts.LogFile         text, default fullfile(tests,'runAllTests_results.txt')
%   Outputs
%     results  (1,N) matlab.unittest.TestResult  full result array
%              (pass/fail/duration per test; log written to LogFile)
%
%   Claude (Fable 5), 2026-08-07.
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    opts.SkipPerformance (1,1) logical = false
    opts.LogFile {mustBeTextScalar} = ""
end

root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'tests'));

logFile = char(opts.LogFile);
if isempty(logFile)
    logFile = fullfile(root, 'tests', 'runAllTests_results.txt');
end
if isfile(logFile), delete(logFile); end
diary(logFile);
cleanupDiary = onCleanup(@() diary('off'));
v = shLowLevel.version();
fprintf('%s v%s test run - %s\n', v.Name, v.Version, char(datetime('now', ...
    Format = 'yyyy-MM-dd HH:mm:ss')));
fprintf('MATLAB %s on %s\n', version, computer);

suites = {'testCorrectness', 'testContract', 'testRobustness'};
if ~opts.SkipPerformance
    suites{end+1} = 'testPerformance';
end

results = matlab.unittest.TestResult.empty;
for k = 1:numel(suites)
    fprintf('\n=== %s ===\n', suites{k});
    results = [results, runtests(suites{k})]; %#ok<AGROW>
end

fprintf('\n================ summary ================\n');
fprintf('  total:  %d\n  passed: %d\n  failed: %d\n  incomplete: %d\n', ...
    numel(results), nnz([results.Passed]), nnz([results.Failed]), ...
    nnz([results.Incomplete]));
if any([results.Failed])
    fprintf('  FAILED tests:\n');
    fprintf('    %s\n', results([results.Failed]).Name);
end

% compact per-test table at the end of the log (robust against console
% truncation: the file always holds the complete run)
fprintf('\n---- per-test results ----\n');
for k = 1:numel(results)
    if results(k).Failed, st = 'FAIL';
    elseif results(k).Incomplete, st = 'INC ';
    else, st = 'ok  ';
    end
    fprintf('%s  %7.2fs  %s\n', st, results(k).Duration, results(k).Name);
end
diary('off');
fprintf('\nFull log written to:\n  %s\n', logFile);
end
