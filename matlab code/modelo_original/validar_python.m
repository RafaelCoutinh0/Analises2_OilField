%% ============================================================
% VALIDAÇÃO MATLAB x PYTHON
% ============================================================

import casadi.*

ref = load('validacao_python.mat');

X_py = ref.X_python;
Z_py = ref.Z_python;
Inputs_py = ref.Inputs_python;

qtd_pts = 25;
tfinal = 1000;

local_grid = linspace(0,tfinal,qtd_pts);

%% Recupera entradas aplicadas por bloco

n_blocos = size(Inputs_py,2) / qtd_pts;

U_blocos = zeros(6,n_blocos);

for k = 1:n_blocos
    idx = (k-1)*qtd_pts + 1;
    U_blocos(:,k) = Inputs_py(:,idx);
end

%% Integrador MATLAB/CasADi

F_val = integrator( ...
    'F_validation', ...
    'idas', ...
    dae, ...
    0, ...
    local_grid);

%% Condições iniciais

x_atual = x_ss(:);
z_atual = z_ss(:);

X_mat = [];
Z_mat = [];

%% Simulação usando EXATAMENTE os u do Python

for k = 1:n_blocos

    u_k = U_blocos(:,k);

    res = F_val( ...
        'x0', x_atual, ...
        'z0', z_atual, ...
        'p', u_k);

    X_bloco = full(res.xf);
    Z_bloco = full(res.zf);

    X_mat = [X_mat X_bloco];
    Z_mat = [Z_mat Z_bloco];

    % Estado final vira condição inicial
    x_atual = X_bloco(:,end);
    z_atual = Z_bloco(:,end);

end

%% ============================================================
% COMPARAÇÃO
% ============================================================

erro_X = X_mat - X_py;
erro_Z = Z_mat - Z_py;

fprintf('\n========================================\n');
fprintf('      MATLAB x PYTHON\n');
fprintf('========================================\n');

fprintf('Erro máximo X  : %.12e\n', ...
    max(abs(erro_X(:))));

fprintf('Erro máximo Z  : %.12e\n', ...
    max(abs(erro_Z(:))));

RMSE_X = sqrt(mean(erro_X(:).^2));
RMSE_Z = sqrt(mean(erro_Z(:).^2));

fprintf('RMSE X         : %.12e\n', RMSE_X);
fprintf('RMSE Z         : %.12e\n', RMSE_Z);

%% Erro por variável

nomes_x = {
    'P_man'
    'q_tr'
    'P_fbhp_1'
    'P_choke_1'
    'q_mean_1'
    'P_fbhp_2'
    'P_choke_2'
    'q_mean_2'
};

nomes_z = {
    'P_intake_1'
    'dP_BCS_1'
    'P_intake_2'
    'dP_BCS_2'
};

fprintf('\nERRO MÁXIMO POR ESTADO\n');

for i = 1:size(X_mat,1)

    erro_i = max(abs(erro_X(i,:)));

    fprintf('%-15s : %.12e\n', ...
        nomes_x{i}, erro_i);

end

fprintf('\nERRO MÁXIMO POR VARIÁVEL ALGÉBRICA\n');

for i = 1:size(Z_mat,1)

    erro_i = max(abs(erro_Z(i,:)));

    fprintf('%-15s : %.12e\n', ...
        nomes_z{i}, erro_i);

end

t = ref.grid_python(:)';

figure

plot(t, X_py(1,:), '-', 'LineWidth', 1.5)
hold on
plot(t, X_mat(1,:), '--', 'LineWidth', 1.5)

xlabel('Tempo [s]')
ylabel('P_{man} [bar]')
legend('Python','MATLAB')
title('Validação - Pressão do Manifold')

set(gca,'XGrid','on','YGrid','on');


figure

plot(t, X_py(2,:), '-', 'LineWidth', 1.5)
hold on
plot(t, X_mat(2,:), '--', 'LineWidth', 1.5)

xlabel('Tempo [s]')
ylabel('q_{tr} [m^3/h]')
legend('Python','MATLAB')
title('Validação - Vazão do Manifold')

set(gca,'XGrid','on','YGrid','on');


figure

plot(t, X_py(3,:), '-', 'LineWidth', 1.5)
hold on
plot(t, X_mat(3,:), '--', 'LineWidth', 1.5)

xlabel('Tempo [s]')
ylabel('P_{FBHP,1} [bar]')
legend('Python','MATLAB')
title('Validação - FBHP Poço 1')

set(gca,'XGrid','on','YGrid','on');

eps_val = 1e-12;

erro_rel_X = abs(X_mat - X_py) ./ ...
    max(abs(X_py),eps_val);

erro_rel_Z = abs(Z_mat - Z_py) ./ ...
    max(abs(Z_py),eps_val);

fprintf('\nErro relativo máximo X = %.6e %%\n', ...
    100*max(erro_rel_X(:)));

fprintf('Erro relativo máximo Z = %.6e %%\n', ...
    100*max(erro_rel_Z(:)));