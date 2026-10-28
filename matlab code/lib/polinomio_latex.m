function txt = polinomio_latex(coef)
%POLINOMIO_LATEX Exprime polinômio numérico em potências de s, sem Symbolic.
    coef = coef(:)';
    grau = numel(coef)-1;
    escala = max(abs(coef));
    if escala==0, txt = '0'; return; end
    partes = {};
    for k = 1:numel(coef)
        a = coef(k); g = grau-k+1;
        if abs(a)<escala*1e-13, continue; end
        mag = sprintf('%.4g',abs(a));
        if g==0
            parcela = mag;
        elseif g==1
            parcela = [mag '\,s'];
        else
            parcela = sprintf('%s\\,s^{%d}',mag,g);
        end
        if isempty(partes)
            if a<0
                partes{end+1} = ['-' parcela]; %#ok<AGROW>
            else
                partes{end+1} = parcela; %#ok<AGROW>
            end
        else
            if a<0
                partes{end+1} = [' - ' parcela]; %#ok<AGROW>
            else
                partes{end+1} = [' + ' parcela]; %#ok<AGROW>
            end
        end
    end
    if isempty(partes), txt = '0'; else, txt = [partes{:}]; end
end
