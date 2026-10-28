%% 1. APRESENTAÇÃO DO SISTEMA - equivalente a 1.sistema.ipynb
% Modelo original de dois poços com BCS e manifold. 8 EDOs + 4 algébricas.
% Pulsos independentes nas 5 entradas manipuladas, sem alterar p_topside.
% Gere as figuras de entradas, estados e (extra DAE) variáveis algébricas.
clearvars; close all;
raiz = configurar_projeto();
nomes = nomes_modelo();
mdl = sym_dynamics();

fprintf('\n=========================================================\n');
fprintf(' 1. SISTEMA NÃO LINEAR - MODELO DAE BCS (2 POÇOS)\n');
fprintf('=========================================================\n');
fprintf('Estados diferenciais: %d | Algébricas: %d\n',mdl.nx,mdl.nz);
fprintf('Entradas totais: %d | Manipuladas: %d\n',mdl.nu,numel(nomes.idx_u));
fprintf('Pressão topside fixada em %.3f bar.\n',mdl.u0(2)/1e5);
fprintf('Máximo resíduo do ponto estacionário: %.3e\n', ...
    max(abs(mdl.fsolve_residual)));

% Na referência, t=0..5 h e 1500 pontos. Aqui: segundos, conforme a DAE.
T_final = 600; n_amostras = 451;
t = linspace(0,T_final,n_amostras);
for k = 1:numel(nomes.idx_u)
    U = perfil_entradas(t,mdl.u0,nomes.idx_u,k,'pulso');
    fprintf('\nPulso em %s: entradas e saídas do sistema não linear\n',nomes.u{k});
    [X,Z] = simular_dae(mdl,t,U);
    % IMPORTANTE: salva resultados antes de criar qualquer objeto grafico.
    % Se o MATLAB falhar ao exportar PNG, pode-se recuperar as figuras sem
    % reintegrar o modelo DAE via regenerar_figuras_etapa_1.
    arquivo_dados = fullfile(raiz,'outputs','simulacoes_nonlineares', ...
        sprintf('pulse-%s.mat',nomes.u{k}));
    save(arquivo_dados,'t','U','X','Z','-v7');
    plotar_entradas(t,U,fullfile(raiz,'figures','nonlinear', ...
        sprintf('pulse-%s-inputs.png',nomes.u{k})));
    plotar_saidas(t,X,fullfile(raiz,'figures','nonlinear', ...
        sprintf('pulse-%s-outputs.png',nomes.u{k})));
    plotar_algebricas(t,Z,fullfile(raiz,'figures','nonlinear', ...
        sprintf('pulse-%s-algebraics.png',nomes.u{k})));
    fprintf('  P_man final = %.5f bar; q_tr final = %.5f m^3/h\n', ...
        X(1,end),X(2,end));
end
fprintf('\nFiguras salvas em figures/nonlinear/\n');
