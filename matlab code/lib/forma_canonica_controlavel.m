function [Ac,Bc,Cc,Dc] = forma_canonica_controlavel(num,den)
%FORMA_CANONICA_CONTROLAVEL Companion no mesmo arranjo do scipy.tf2ss.
%  Para um G(s) próprio: A(1,:)=-[a_(n-1) ... a_0],
%  B=[1;0;...], C=[b_(n-1) ... b_0]-D*[a_(n-1) ... a_0].
%  Coeficientes passados em potências DECRESCENTES de s.
    num = num(:)'; den = den(:)';
    assert(~isempty(den) && den(1)~=0,'Denominador inválido');
    a0 = den(1);
    den = den/a0;num = num/a0;
    ordem = numel(den)-1;
    if numel(num)>ordem+1
        error('BCS:Canonica','Função de transferência imprópria.');
    end
    b = [zeros(1,ordem+1-numel(num)),num];
    Dc = b(1);
    if ordem==0
        Ac = zeros(0); Bc = zeros(0,1); Cc = zeros(1,0);
        return;
    end
    Ac = zeros(ordem);
    Ac(1,:) = -den(2:end);
    if ordem>1, Ac(2:end,1:end-1) = eye(ordem-1); end
    Bc = zeros(ordem,1); Bc(1) = 1;
    Cc = b(2:end)-Dc*den(2:end);
end
