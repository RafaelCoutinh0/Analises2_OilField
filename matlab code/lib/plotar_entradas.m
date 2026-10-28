function plotar_entradas(t, U, nome_arquivo)
%PLOTAR_ENTRADAS Cinquenta por cento similar ao layout dos notebooks: 1 linha.
    n = nomes_modelo();
    fig = figure('Name','Entradas do sistema','Color','w', ...
        'Position',[40 60 1750 370]);
    for k = 1:numel(n.idx_u)
        subplot(1,numel(n.idx_u),k);
        stairs(t, U(n.idx_u(k),:),'LineWidth',1.6);
        xlabel('Tempo / s');
        ylabel(sprintf('%s / %s',n.u_labels{k},n.u_units{k}), ...
            'Interpreter','tex');
        title(n.u{k},'Interpreter','none');grid on;
    end
    salvar_figura(fig,nome_arquivo);
end
