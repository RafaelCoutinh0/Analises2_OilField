function [Y,X] = simular_passo_siso(A,B,C,D,t)
%SIMULAR_PASSO_SISO Resposta exata a degrau unitário via expm ampliada.
%  Permite matrizes complexas em uma base modal/de autovetores.
    t = t(:)'; dt = t(2)-t(1);
    nx = size(A,1);
    M = expm([A,B;zeros(1,nx+1)]*dt);
    Ad = M(1:nx,1:nx); Bd = M(1:nx,end);
    X = zeros(nx,numel(t),'like',A);
    Y = zeros(1,numel(t),'like',A);
    for k = 1:numel(t)
        Y(k) = C*X(:,k)+D;
        if k<numel(t), X(:,k+1) = Ad*X(:,k)+Bd; end
    end
end
