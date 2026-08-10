function nse = nashSutcliffe(ref, model)
%NASHSUTCLIFFE Nash-Sutcliffe efficiency of a model against a reference.
%
%   NSE = shx.nashSutcliffe(REF, MODEL) computes the hydrology-standard
%   skill score NSE = 1 - sum((model-ref).^2) / sum((ref-mean(ref)).^2)
%   column-wise. NSE = 1 is perfect, 0 means "no better than the
%   reference mean", negative is worse than the mean. More informative
%   than correlation because bias and amplitude errors are punished.
%   NaN pairs are removed per column.
%
%   Outputs
%     nse        (1,K) double  efficiency per column (K = 1 for vectors)
%
%   Example
%     nse = shx.nashSutcliffe(avg(1,:), avg(2,:));   % e.g. SH vs mascon basin
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5),
%   2026-08-08 (v2.6.0).
arguments
    ref double
    model double
end
if isvector(ref), ref = ref(:); end
if isvector(model), model = model(:); end
if ~isequal(size(ref), size(model))
    error('shx:nashSutcliffe:sizeMismatch', ...
        'REF and MODEL must share the same size.');
end
K = size(ref, 2);
nse = nan(1, K);
for k = 1:K
    use = isfinite(ref(:, k)) & isfinite(model(:, k));
    r = ref(use, k); m = model(use, k);
    if numel(r) < 3, continue, end
    den = sum((r - mean(r)).^2);
    if den > 0, nse(k) = 1 - sum((m - r).^2) / den; end
end
end
