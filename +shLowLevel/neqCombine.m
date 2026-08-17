function [x, out] = neqCombine(neqs, opts)
%NEQCOMBINE VCE-weighted combination on the normal-equation level.
%
%   [X, OUT] = shLowLevel.neqCombine(NEQS) combines K normal equations
%   into one solution the COST-G way - on the normal-equation level,
%   BEFORE solving - with one variance component per contribution,
%   estimated by Foerstner/Koch iteration with partial redundancies:
%       N     = sum_i (1/s_i^2) N_i,  b = sum_i (1/s_i^2) b_i
%       x     = N^-1 b
%       Om_i  = ltpl_i - 2 x' b_i + x' N_i x     (= v_i' P_i v_i)
%       r_i   = nobs_i - (1/s_i^2) tr(N^-1 N_i)  (partial redundancy)
%       s_i^2 = Om_i / r_i
%   The redundancy invariant sum(r_i) = sum(nobs_i) - P holds at every
%   iteration and is asserted. This is the rigorous sibling of
%   shLowLevel.combineCenters, which combines on the SOLUTION level
%   with approximate noise shapes; here the full information content
%   of each contribution enters (Kvas 2019, Sec. 2.1-2.2 formalism).
%
%   CONSISTENT BACKGROUNDS REQUIRED, stated loudly: every N_i, b_i and
%   ltpl_i must refer to the SAME background models and parametrization
%   (same shLowLevel.shIndex ordering). Mixing backgrounds combines
%   apples with oranges and no algebra will warn you - align first.
%
%   Inputs
%     neqs  (1 x K) struct  one contribution each, in EITHER form:
%           N (P x P double), b (P x 1 double), and for VCE also
%           ltpl (1 x 1 double, l'Pl of the contribution) and
%           nobs (1 x 1 double, its observation count); OR directly a
%           shLowLevel.readSINEX result with kind = 'NEQ' (fields M, x,
%           stats): then N = M, b = M*x, ltpl/nobs from stats (the
%           +SOLUTION/STATISTICS block; v3.20 readSINEX parses it).
%           An optional field name labels the contribution.
%
%   Options
%     Weights ([])      (K x 1) double  FIXED weights w_i = 1/s_i^2:
%           skips the VCE iteration entirely (then ltpl/nobs are not
%           needed); [] runs VCE, which REQUIRES ltpl and nobs on every
%           contribution and errors loudly otherwise - no silent noise
%           assumptions
%     MaxIter (20)      (1 x 1) VCE iteration cap
%     Tol (1e-6)        (1 x 1) relative change of all s_i^2 to stop
%     Covariance (false) also return Cx = N^-1 (P x P inverse - costly
%           for large P, off by default)
%
%   Outputs
%     x    (P x 1) double  combined solution
%     out  (1 x 1) struct  sigma2 (K x 1) final variance factors (the
%          fixed 1./Weights when given), w (K x 1) = 1./sigma2,
%          redundancy (K x 1) partial redundancies (NaN in fixed mode),
%          nIter (1 x 1), converged (1 x 1 logical), N (P x P) combined
%          normal matrix, Cx (P x P, Covariance=true only), names
%          (K x 1 string)
%
%   Example  % two centers, same month, same background
%     a = shLowLevel.readSINEX(fileA, Index = idx);
%     b = shLowLevel.readSINEX(fileB, Index = idx);
%     [x, out] = shLowLevel.neqCombine([a, b]);
%     disp(out.sigma2')          % who carries how much noise power
%
%   Numerics pre-validated in Python (tools/dev/validate_neqcombine.py,
%   5 checks: fixed == stacked GLS, factor recovery, redundancy
%   invariant, combined beats singles, equal-noise symmetry).
%
%   Reference: Koch (1999) variance component estimation; Foerstner
%   (1979); Kvas, TU Graz PhD thesis (2019), Secs. 2.1-2.2; the
%   COST-G combination concept (Jaeggi et al.).
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-17, 22:40 UTC.

arguments
    neqs (1,:) struct
    opts.Weights (:,1) double = []
    opts.MaxIter (1,1) double {mustBeInteger, mustBePositive} = 20
    opts.Tol (1,1) double {mustBePositive} = 1e-6
    opts.Covariance (1,1) logical = false
end

K = numel(neqs);
if K < 1
    error('shLowLevel:neqCombine:empty', 'Need at least one contribution.');
end
fixed = ~isempty(opts.Weights);
if fixed && numel(opts.Weights) ~= K
    error('shLowLevel:neqCombine:badWeights', ...
        'Weights has %d entries for %d contributions.', numel(opts.Weights), K);
end

% ---- normalize contributions to N, b (, ltpl, nobs)
Nc = cell(1, K); bc = cell(1, K);
ltpl = NaN(K, 1); nobs = NaN(K, 1);
names = strings(K, 1);
for i = 1:K
    e = neqs(i);
    if isfield(e, 'N') && ~isempty(e.N)
        Nc{i} = e.N; bc{i} = e.b(:);
        if isfield(e, 'ltpl') && ~isempty(e.ltpl), ltpl(i) = e.ltpl; end
        if isfield(e, 'nobs') && ~isempty(e.nobs), nobs(i) = e.nobs; end
    elseif isfield(e, 'M') && ~isempty(e.M)
        if ~isfield(e, 'kind') || ~strcmp(e.kind, 'NEQ')
            error('shLowLevel:neqCombine:notNEQ', ...
                'Contribution %d: readSINEX struct without kind=''NEQ''.', i);
        end
        Nc{i} = e.M; bc{i} = e.M * e.x(:);
        if isfield(e, 'stats') && ~isempty(e.stats)
            if isfield(e.stats, 'wsos'), ltpl(i) = e.stats.wsos; end
            if isfield(e.stats, 'nobs'), nobs(i) = e.stats.nobs; end
        end
    else
        error('shLowLevel:neqCombine:badInput', ...
            'Contribution %d carries neither N/b nor a readSINEX M/x.', i);
    end
    if isfield(e, 'name') && ~isempty(e.name), names(i) = string(e.name); end
end
P = size(Nc{1}, 1);
for i = 1:K
    if ~isequal(size(Nc{i}), [P P]) || numel(bc{i}) ~= P
        error('shLowLevel:neqCombine:sizeMismatch', ...
            'Contribution %d: N must be %d x %d and b %d x 1.', i, P, P, P);
    end
end
if ~fixed && (any(isnan(ltpl)) || any(isnan(nobs)))
    error('shLowLevel:neqCombine:needStats', ...
        ['VCE needs ltpl and nobs on every contribution (missing on %d ' ...
         'of %d). Provide them (SINEX +SOLUTION/STATISTICS; readSINEX ' ...
         'parses it since v3.20) or pass fixed Weights.'], ...
        nnz(isnan(ltpl) | isnan(nobs)), K);
end

% ---- iterate
if fixed
    s2 = 1 ./ opts.Weights;
    nIter = 1; converged = true;
    red = NaN(K, 1);
else
    s2 = ones(K, 1);
    converged = false;
    red = NaN(K, 1);
    for nIter = 1:opts.MaxIter
        [N, b] = accumulate(Nc, bc, s2);
        x = N \ b;
        Ninv = inv(N);
        s2new = zeros(K, 1);
        for i = 1:K
            Om = ltpl(i) - 2 * (x.' * bc{i}) + x.' * (Nc{i} * x);
            red(i) = nobs(i) - trace(Ninv * Nc{i}) / s2(i); %#ok<MINV>
            if red(i) <= 0
                error('shLowLevel:neqCombine:nonpositiveRedundancy', ...
                    ['Contribution %d has partial redundancy %.3g <= 0 - ' ...
                     'fewer observations than effective parameters. VCE is ' ...
                     'undefined; fix nobs or pass fixed Weights.'], i, red(i));
            end
            s2new(i) = Om / red(i);
        end
        % invariant: sum(red) == sum(nobs) - P (asserted, cheap insurance)
        assert(abs(sum(red) - (sum(nobs) - P)) < 1e-6 * max(1, sum(nobs)), ...
            'shLowLevel:neqCombine:redundancyInvariant', ...
            'Redundancy invariant violated: numerical breakdown.');
        if max(abs(s2new - s2) ./ s2) < opts.Tol
            s2 = s2new; converged = true; break
        end
        s2 = s2new;
    end
end
[N, b] = accumulate(Nc, bc, s2);
x = N \ b;

out = struct('sigma2', s2, 'w', 1 ./ s2, 'redundancy', red, ...
    'nIter', nIter, 'converged', converged, 'N', N, 'names', names);
if opts.Covariance
    out.Cx = inv(N);
end
if ~fixed && ~converged
    warning('shLowLevel:neqCombine:noConvergence', ...
        'VCE not converged after %d iterations (last rel. change of s^2 above %.1e).', ...
        opts.MaxIter, opts.Tol);
end
end

% ---------------------------------------------------------------- local
function [N, b] = accumulate(Nc, bc, s2)
N = Nc{1} / s2(1); b = bc{1} / s2(1);
for i = 2:numel(Nc)
    N = N + Nc{i} / s2(i);
    b = b + bc{i} / s2(i);
end
end
