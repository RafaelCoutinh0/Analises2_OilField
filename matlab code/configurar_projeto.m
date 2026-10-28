function raiz = configurar_projeto()
%CONFIGURAR_PROJETO Configura pastas, CasADi e estilo dos graficos.
    raiz = fileparts(mfilename('fullpath'));
    addpath(raiz);

    % Resolve dependencias apenas na pasta do projeto (ou aninhadas).
    pastas = {'lib', 'modelo'};
    essenciais = {'simular_dae.m', 'sym_dynamics.m'};
    for k = 1:numel(pastas)
        pasta = fullfile(raiz, pastas{k});
        if ~isfolder(pasta)
            alternativa = fullfile(raiz, 'Analise_BCS_MATLAB', pastas{k});
            if isfolder(alternativa), pasta = alternativa; end
        end
        if isfolder(pasta)
            addpath(pasta, '-begin');
        elseif ~isfile(fullfile(raiz, essenciais{k}))
            error('BCS:EstruturaIncompleta', ...
                ['Nao encontrei %s (nem a pasta %s). ', ...
                 'Extraia o ZIP completo preservando as subpastas.\nRaiz: %s'], ...
                essenciais{k}, pastas{k}, raiz);
        end
    end

    % Reconhecimento automatico da instalacao Windows que ja funcionou
    % neste computador. Nao exige import casadi.* no Command Window.
    casadi_ok = existe_casadi();
    if ~casadi_ok
        candidatos = {};
        ambiente = getenv('CASADI_MATLAB_PATH');
        if ~isempty(ambiente), candidatos{end+1} = ambiente; end %#ok<AGROW>
        perfil = getenv('USERPROFILE');
        if ~isempty(perfil)
            pasta_matlab = fullfile(perfil, 'Documents', 'MATLAB');
            % Primeiro tenta a versao especifica informada pelo usuario.
            candidatos{end+1} = fullfile(pasta_matlab, ...
                'casadi-3.8.1-windows64-matlab2018b'); %#ok<AGROW>
            % Outras versoes da distribuicao CasADi para MATLAB.
            instalacoes = dir(fullfile(pasta_matlab, 'casadi-*-windows64-matlab*'));
            for k = 1:numel(instalacoes)
                if instalacoes(k).isdir
                    candidatos{end+1} = fullfile(pasta_matlab, instalacoes(k).name); %#ok<AGROW>
                end
            end
        end
        for k = 1:numel(candidatos)
            if isfolder(candidatos{k})
                addpath(candidatos{k}, '-begin');
                if existe_casadi()
                    fprintf('[configurar_projeto] CasADi: %s\n', candidatos{k});
                    casadi_ok = true;
                    break;
                end
            end
        end
    end
    if ~casadi_ok
        error('BCS:CasADi', [ ...
            'CasADi nao encontrado. Use addpath(''C:\\caminho\\casadi-...''), ', ...
            'ou defina CASADI_MATLAB_PATH antes de executar os scripts.']);
    end

    if exist('ss','file')~=2 || exist('c2d','file')~=2
        error('BCS:ControlToolbox','Requer Control System Toolbox (ss, c2d).');
    end
    if exist('fsolve','file')~=2
        error('BCS:OptimizationToolbox','Requer Optimization Toolbox (fsolve).');
    end

    dirs = {'outputs', 'outputs/simulacoes_nonlineares', ...
        'figures/nonlinear', 'figures/linear', 'figures/comparacao', ...
        'figures/estabilidade', 'figures/formas_canonicas', ...
        'figures/discretizacao'};
    for k = 1:numel(dirs)
        destino = fullfile(raiz, dirs{k});
        if ~isfolder(destino), mkdir(destino); end
    end

    % Grafico claro, mesmo quando o MATLAB usa tema escuro.
    set(groot, 'defaultFigureColor', [1 1 1]);
    set(groot, 'defaultAxesColor', [1 1 1]);
    set(groot, 'defaultAxesXColor', [0 0 0]);
    set(groot, 'defaultAxesYColor', [0 0 0]);
    set(groot, 'defaultAxesZColor', [0 0 0]);
    set(groot, 'defaultTextColor', [0 0 0]);
    set(groot, 'defaultAxesGridColor', [0.72 0.72 0.72]);
end

function ok = existe_casadi()
    ok = false;
    try
        teste = casadi.MX.sym('casadi_teste'); %#ok<NASGU>
        ok = true;
    catch
        % Checara os candidatos de instalacao em seguida.
    end
end
