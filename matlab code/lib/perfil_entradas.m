function U = perfil_entradas(t, u0, idx_u, qual, modo)
%PERFIL_ENTRADAS Retorna entradas físicas, linhas: as 6 entradas da DAE.
%  Modo 'pulso': um pulso em uma das 5 entradas manipuladas.
%  Modo 'comparacao': dois patamares em f_BCS_1.
    if nargin < 5, modo = 'pulso'; end
    t = t(:)';
    U = repmat(u0(:), 1, numel(t));
    t_ini = t(1); t_fim = t(end); dur = t_fim - t_ini;
    switch lower(modo)
        case 'pulso'
            delta = [2.0, 2.0, 0.03, 2.0, 0.03];
            mask = t >= t_ini+0.2*dur & t < t_ini+0.55*dur;
            U(idx_u(qual),mask) = u0(idx_u(qual)) + delta(qual);
        case 'comparacao'
            % f_BCS_1 = 50 Hz no nominal, depois 52 Hz e 59 Hz.
            idx = idx_u(2);
            mask_pequeno = t >= t_ini+0.12*dur & t < t_ini+0.33*dur;
            mask_grande  = t >= t_ini+0.33*dur & t < t_ini+0.70*dur;
            U(idx, mask_pequeno) = u0(idx) + 2;
            U(idx, mask_grande)  = u0(idx) + 9;
        otherwise
            error('BCS:Perfil', 'Modo de perturbação desconhecido: %s', modo);
    end
end
