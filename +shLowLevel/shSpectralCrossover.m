function [nCrossover, degreeInterp] = shSpectralCrossover(spec)
%SHSPECTRALCROSSOVER Degree at which the error spectrum first exceeds the
%   signal (degree-amplitude) spectrum -- a standard guide for choosing a
%   truncation degree or Gaussian filter radius.
%
%   [NCROSSOVER, DEGREEINTERP] = SHSPECTRALCROSSOVER(SPEC) requires SPEC
%   to have both '.degAmplitude' and '.errAmplitude' fields (i.e. SPEC was
%   produced by SHDEGREERMS called with 'sigmaC'/'sigmaS'). NCROSSOVER is
%   the first integer degree n where errAmplitude(n) > degAmplitude(n).
%   DEGREEINTERP is a sub-degree estimate via linear interpolation in
%   log-amplitude space between that degree and the one before it (more
%   useful for choosing a smooth filter radius than the integer value).
%
%   Returns NaN for both outputs if the error spectrum never exceeds the
%   signal spectrum over the available degree range.
%
%   Claude (Sonnet 4.6), 2026-07-11; merged into +shLowLevel: Claude (Fable 5), 2026-08-07.
%   Outputs
%     nc         (1 x 1) double   first degree where the error spectrum exceeds the signal (NaN if none)
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

if ~isfield(spec, 'errAmplitude')
    error('shSpectralCrossover:missingField', ...
        'SPEC has no ''errAmplitude'' field -- call shDegreeRMS with sigmaC/sigmaS.');
end

n = spec.degree;
sig = spec.degAmplitude;
err = spec.errAmplitude;

idx = find(err > sig, 1, 'first');
if isempty(idx) || idx == 1
    nCrossover = NaN;
    degreeInterp = NaN;
    return;
end

nCrossover = n(idx);

% log-space linear interpolation between idx-1 (signal above error) and idx
logSig1 = log(sig(idx-1)); logErr1 = log(err(idx-1));
logSig2 = log(sig(idx));   logErr2 = log(err(idx));
% find fraction f in [0,1] where logSig - logErr crosses zero
d1 = logSig1 - logErr1;
d2 = logSig2 - logErr2;
f = d1 / (d1 - d2);
degreeInterp = n(idx-1) + f * (n(idx) - n(idx-1));

end
