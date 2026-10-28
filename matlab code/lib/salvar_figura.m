function sucesso = salvar_figura(fig, arquivo)
%SALVAR_FIGURA Exportacao robusta de graficos para fundo branco.
% Nunca perde os resultados numericos nem interrompe a DAE se a camada
% grafica falhar. Para recuperar as figuras, use regenerar_figuras_etapa_1.
% O resultado logico 'sucesso' e opcional (compatibilidade com scripts).
    sucesso = false;
    [pasta,nome,~] = fileparts(arquivo);
    if ~isempty(pasta) && ~isfolder(pasta), mkdir(pasta); end

    if ~isgraphics(fig, 'figure')
        warning('BCS:FiguraInvalida', ...
            'Handle de figura invalido antes de exportar: %s',arquivo);
        return;
    end

    % Garante legibilidade no tema escuro sem alterar as curvas.
    try
        set(fig, 'Color', [1 1 1], 'InvertHardcopy', 'on');
        axes_da_figura = findall(fig, 'Type', 'axes');
        for k = 1:numel(axes_da_figura)
            ax = axes_da_figura(k);
            set(ax, 'Color', [1 1 1], ...
                'XColor', [0 0 0], 'YColor', [0 0 0], ...
                'ZColor', [0 0 0], ...
                'GridColor', [0.72 0.72 0.72], 'GridAlpha', 0.35);
        end
        objetos_texto = findall(fig, 'Type', 'text');
        for k = 1:numel(objetos_texto)
            set(objetos_texto(k), 'Color', [0 0 0]);
        end
        legendas = findall(fig, 'Type', 'legend');
        for k = 1:numel(legendas)
            set(legendas(k),'Color',[1 1 1],'TextColor',[0 0 0]);
        end
        barras_cores = findall(fig, 'Type', 'colorbar');
        for k = 1:numel(barras_cores)
            set(barras_cores(k),'Color',[0 0 0]);
        end
        drawnow;
    catch erro
        warning('BCS:EstiloFigura', ...
            'Ajuste de cores/desenho falhou (%s). Tentarei exportar.',erro.message);
    end

    erros = {};
    % Tentativa 1: exportgraphics em baixa resolucao para consumir menos
    % memoria que as figuras grandes originais.
    if isgraphics(fig, 'figure')
        try
            exportgraphics(fig,arquivo,'Resolution',110, ...
                'BackgroundColor','white');
            sucesso = true;
        catch e
            erros{end+1} = sprintf('exportgraphics: %s',e.message); %#ok<AGROW>
        end
    end

    % Tentativa 2: usa outro renderizador e uma exportacao PNG antiga.
    if ~sucesso && isgraphics(fig,'figure')
        try
            set(fig,'Renderer','painters');
            drawnow;
            print(fig,arquivo,'-dpng','-r100');
            sucesso = true;
        catch e
            erros{end+1} = sprintf('print: %s',e.message); %#ok<AGROW>
        end
    end

    % Tentativa 3: a API saveas pode funcionar em alguns sistemas.
    if ~sucesso && isgraphics(fig,'figure')
        try
            saveas(fig,arquivo);
            sucesso = true;
        catch e
            erros{end+1} = sprintf('saveas: %s',e.message); %#ok<AGROW>
        end
    end

    if ~sucesso
        % Salva a figura MATLAB editavel, mesmo que o PNG falhe.
        caminho_fig = fullfile(pasta,[nome '.fig']);
        backup_ok = false;
        if isgraphics(fig, 'figure')
            try
                savefig(fig,caminho_fig);
                backup_ok = true;
            catch e
                erros{end+1} = sprintf('savefig: %s',e.message); %#ok<AGROW>
            end
        end
        if backup_ok
            aviso = sprintf('Tente abrir o backup %s.',caminho_fig);
        else
            aviso = 'Use os dados .mat e regenerar_figuras_etapa_1.';
        end
        warning('BCS:Exportacao', ...
            'Falha na exportacao do PNG %s. %s Detalhes: %s', ...
            arquivo, aviso, strjoin(erros,' | '));
    end

    % Fecha SEMPRE, inclusive se houve erro de renderizacao.
    if isgraphics(fig,'figure')
        try
            close(fig);
        catch
            % Nao deixe uma falha de interface travar a simulacao.
        end
    end
end
