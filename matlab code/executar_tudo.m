%% EXECUTAR_TUDO - replica a ordem dos cinco notebooks.
% Abra esta pasta como Current Folder no MATLAB e execute este script.
clc;
addpath(fileparts(mfilename('fullpath')));
run(fullfile(fileparts(mfilename('fullpath')),'etapa_1_sistema.m'));
run(fullfile(fileparts(mfilename('fullpath')),'etapa_2_linearizacao.m'));
run(fullfile(fileparts(mfilename('fullpath')),'etapa_3_estabilidade.m'));
run(fullfile(fileparts(mfilename('fullpath')),'etapa_4_formas_canonicas.m'));
run(fullfile(fileparts(mfilename('fullpath')),'etapa_5_discretizacao.m'));
fprintf('\nTODAS AS ETAPAS CONCLUÍDAS. Confira figures/ e outputs/.\n');
