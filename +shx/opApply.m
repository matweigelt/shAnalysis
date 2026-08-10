function Yout = opApply(op, Z, t, mode)
%OPAPPLY Apply the month-t filter operator W_t (or W_t') to columns of Z.
%   Yout = shx.opApply(op, Z, t)          computes W_t * Z
%   Yout = shx.opApply(op, Z, t, 'transp') computes W_t' * Z
%
%   W_t = W0_t + S*Ac*(Ac'*S*Ac)^-1 * Ac' * (I - W0_t)  (if constrained)
%   with W0_t = V * diag(g) * Ut, g = lam./(lam + s(t)).
%   O(P^2) per column; no P x P matrix is formed.
%   Outputs
%     xf         (P x cols(x)) double   W_t * x (or W_t' * x in 'transp' mode) without forming the P x P matrix
%
%   Developed by Matthias Weigelt with the help of Claude (Fable 5).

arguments
    op (1,1) struct
    Z double
    t (1,1) double {mustBeInteger, mustBePositive}
    mode {mustBeMember(mode, {'forward','transp'})} = 'forward'
end

g = op.lam ./ (op.lam + op.s(t));

if isfield(op, 'layout') && strcmp(op.layout, 'blocks')
    % block-diagonal operator (no constraints by construction)
    Yout = zeros(size(Z));
    for k = 1:numel(op.blocks)
        b = op.blocks(k);
        if isfield(op, 'sBlocks')
            sk = op.sBlocks(k, t);               % per-order-band VCE (v2.2)
        else
            sk = op.s(t);
        end
        gb = b.lam ./ (b.lam + sk);
        if strcmp(mode, 'forward')
            Yout(b.rows, :) = b.V * (gb .* (b.Ut * Z(b.rows, :)));
        else
            Yout(b.rows, :) = b.Ut' * (gb .* (b.V' * Z(b.rows, :)));
        end
    end
    return
end

if strcmp(mode, 'forward')
    Y0 = op.V * (g .* (op.Ut * Z));
    if isempty(op.Ac)
        Yout = Y0;
    else
        Yout = Y0 + op.SA * (op.M \ (op.Ac' * (Z - Y0)));
    end
else
    if isempty(op.Ac)
        Yout = op.Ut' * (g .* (op.V' * Z));
    else
        % W' = W0' + (I - W0') * Ac * M^-1 * SA'   (M symmetric)
        Q = op.Ac * (op.M \ (op.SA' * Z));
        Zq = Z - Q;
        Yout = op.Ut' * (g .* (op.V' * Zq)) + Q;
    end
end
end
