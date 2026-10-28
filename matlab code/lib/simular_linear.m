function [X, Z] = simular_linear(lin, t, U)
%SIMULAR_LINEAR Exatamente ZOH por exponencial de matriz, em desvio dimensional.
%  Parte do estado estacionário. Mantém a consistência das algébricas de 1a ordem.
    t = t(:)';
    dt = diff(t);
    assert(all(abs(dt-dt(1)) < 1e-8*max(1,dt(1))), ...
        'A malha de tempo deve ser equiespaçada.');
    nx = size(lin.A,1); nu = size(lin.B,2);
    E = expm([lin.A,lin.B;zeros(nu,nx+nu)]*dt(1));
    Ad = E(1:nx,1:nx); Bd = E(1:nx,nx+1:end);
    X_desv = zeros(nx,numel(t));
    dU = U(lin.idx_u,:) - lin.u_ss;
    for k = 1:numel(t)-1
        X_desv(:,k+1) = Ad*X_desv(:,k) + Bd*dU(:,k);
    end
    X = X_desv + lin.x_ss;
    Z = lin.Zx*X_desv + lin.Zu*dU + lin.z_ss;
end
