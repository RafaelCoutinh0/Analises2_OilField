function plotar_saidas(t, X, nome_arquivo)
%PLOTAR_SAIDAS 2 linhas x 4 colunas: os oito estados físicos do modelo BCS.
    n = nomes_modelo();
    fig = figure('Name','Saídas do sistema','Color','w', ...
        'Position',[40 50 1700 680]);
    for k = 1:numel(n.x)
        subplot(2,4,k);
        plot(t,X(k,:),'LineWidth',1.6);
        xlabel('Tempo / s');
        ylabel(sprintf('%s / %s',n.x_labels{k},n.x_units{k}), ...
            'Interpreter','tex');
        title(n.x{k},'Interpreter','none');grid on;
    end
    salvar_figura(fig,nome_arquivo);
end
