function plotar_algebricas(t, Z, nome_arquivo)
%PLOTAR_ALGEBRICAS Gráficos adicionais específicos de um modelo DAE.
    n = nomes_modelo();
    fig = figure('Name','Variáveis algébricas','Color','w', ...
        'Position',[80 80 1240 650]);
    for k = 1:numel(n.z)
        subplot(2,2,k);
        plot(t,Z(k,:),'LineWidth',1.6);
        xlabel('Tempo / s');
        ylabel(sprintf('%s / %s',n.z_labels{k},n.z_units{k}), ...
            'Interpreter','tex');
        title(n.z{k},'Interpreter','none');grid on;
    end
    salvar_figura(fig,nome_arquivo);
end
