function [X, Z] = simular_dae(mdl, t, U)
%SIMULAR_DAE Modelo ORIGINAL DAE integrado por IDAS com entradas ZOH.
%  x(k+1) e z(k+1) são obtidos a partir dos valores no instante k.
%  É preciso que t seja uma malha equiespaçada.
    import casadi.*
    t = t(:)';
    dt = diff(t);
    assert(all(abs(dt-dt(1)) < 1e-8*max(1,dt(1))), ...
        'A malha de tempo deve ser equiespaçada.');
    assert(size(U,2)==numel(t) && size(U,1)==mdl.nu, ...
        'U deve ter dimensão [n_entradas_DAE x n_tempos].');

    F = integrator('integrador_dae_bcs', 'idas', mdl.dae, 0, dt(1));
    X = zeros(mdl.nx, numel(t));
    Z = zeros(mdl.nz, numel(t));
    X(:,1) = mdl.x_ss(:);
    Z(:,1) = mdl.z_ss(:);

    for k = 1:numel(t)-1
        res = F('x0', X(:,k), 'z0', Z(:,k), 'p', U(:,k));
        X(:,k+1) = full(res.xf);
        Z(:,k+1) = full(res.zf);
        if ~all(isfinite(X(:,k+1))) || ~all(isfinite(Z(:,k+1)))
            error('BCS:DAE', 'Integração DAE inválida no passo %d.', k);
        end
    end
end
