%% ============================================================
% LINEARIZAÇÃO DO MODELO
% ESPAÇO DE ESTADOS EM DESVIO PERCENTUAL
% FUNÇÕES DE TRANSFERÊNCIA
% POLOS E ZEROS
% ============================================================

clear;
clc;
close all;

import casadi.*


%% ============================================================
% CARREGA/MONTA O MODELO NÃO LINEAR
% ============================================================

mdl = sym_dynamics();


%% ============================================================
% VARIÁVEIS DO MODELO
% ============================================================

x = mdl.x;
z = mdl.z;
u = mdl.u;

f = mdl.f;
g = mdl.g;

x_ss = mdl.x_ss;
z_ss = mdl.z_ss;
u0 = mdl.u0;

nx = mdl.nx;
nz = mdl.nz;
nu = mdl.nu;


%% ============================================================
% ORDEM DAS VARIÁVEIS
% ============================================================

% Estados:
%
% x1 = P_man
% x2 = q_tr
% x3 = P_fbhp_1
% x4 = P_choke_1
% x5 = q_mean_1
% x6 = P_fbhp_2
% x7 = P_choke_2
% x8 = q_mean_2
%
% Algébricas:
%
% z1 = P_intake_1
% z2 = dP_BCS_1
% z3 = P_intake_2
% z4 = dP_BCS_2
%
% Entradas:
%
% u1 = f_booster
% u2 = p_topside
% u3 = f_BCS1
% u4 = alpha1
% u5 = f_BCS2
% u6 = alpha2


%% ============================================================
% JACOBIANAS DO DAE
%
% x_dot = f(x,z,u)
% 0     = g(x,z,u)
% ============================================================

f_x = jacobian(f,x);
f_z = jacobian(f,z);
f_u = jacobian(f,u);

g_x = jacobian(g,x);
g_z = jacobian(g,z);
g_u = jacobian(g,u);


%% ============================================================
% FUNÇÃO CASADI PARA AVALIAR TODAS AS JACOBIANAS
% ============================================================

J_fun = Function( ...
    'J_linearizacao', ...
    {x,z,u}, ...
    {f_x,f_z,f_u,g_x,g_z,g_u});


%% ============================================================
% AVALIAÇÃO NO ESTADO ESTACIONÁRIO
% ============================================================

[F_x,F_z,F_u,G_x,G_z,G_u] = ...
    J_fun(x_ss,z_ss,u0);

F_x = full(F_x);
F_z = full(F_z);
F_u = full(F_u);

G_x = full(G_x);
G_z = full(G_z);
G_u = full(G_u);


%% ============================================================
% VERIFICAÇÃO DAS JACOBIANAS
% ============================================================

fprintf('\n============================================\n');
fprintf(' VERIFICAÇÃO DA LINEARIZAÇÃO DAE\n');
fprintf('============================================\n');

fprintf('max|F_z| = %.12e\n',max(abs(F_z(:))));
fprintf('max|G_u| = %.12e\n',max(abs(G_u(:))));


%% ============================================================
% ELIMINAÇÃO DAS VARIÁVEIS ALGÉBRICAS
%
% dx/dt = f_x dx + f_z dz + f_u du
%
% 0 = g_x dx + g_z dz + g_u du
%
% portanto:
%
% dz = -g_z^-1*g_x dx - g_z^-1*g_u du
%
% Ared = f_x - f_z*g_z^-1*g_x
% Bred = f_u - f_z*g_z^-1*g_u
% ============================================================

A_dim = F_x - F_z * (G_z \ G_x);

B_dim_all = F_u - F_z * (G_z \ G_u);


%% ============================================================
% RELAÇÃO DAS VARIÁVEIS ALGÉBRICAS
% ============================================================

Zx = -(G_z \ G_x);

Zu = -(G_z \ G_u);


%% ============================================================
% TODAS AS POSSÍVEIS SAÍDAS
%
% y_all =
%
% [x]
% [z]
% ============================================================

C_all_dim = [ ...
    eye(nx);
    Zx ...
];

D_all_dim = [ ...
    zeros(nx,nu);
    Zu ...
];


%% ============================================================
% REMOÇÃO DE p_topside
%
% Igual à lógica usada no Python:
%
% entradas manipuladas:
%
% 1 f_booster
% 3 f_BCS1
% 4 alpha1
% 5 f_BCS2
% 6 alpha2
% ============================================================

idx_u = [1 3 4 5 6];

B_dim = B_dim_all(:,idx_u);

u_ss = u0(idx_u);


%% ============================================================
% SELEÇÃO DAS SAÍDAS
%
% [x;z]:
%
% 1  = P_man
% 2  = q_tr
% 3  = P_fbhp_1
% 4  = P_choke_1
% 5  = q_mean_1
% 6  = P_fbhp_2
% 7  = P_choke_2
% 8  = q_mean_2
% 9  = P_intake_1
% 10 = dP_BCS_1
% 11 = P_intake_2
% 12 = dP_BCS_2
%
% Aqui estou usando:
%
% P_man
% P_fbhp_1
% P_fbhp_2
% P_intake_1
% dP_BCS_1
% P_intake_2
% dP_BCS_2
% ============================================================

idx_y = [ ...
    1 ...
    3 ...
    6 ...
    9 ...
    10 ...
    11 ...
    12 ...
];


%% ============================================================
% Caso queira reproduzir os P_choke do código Python antigo,
% use:
%
% idx_y = [1 4 7 9 10 11 12];
% ============================================================


%% ============================================================
% PONTOS ESTACIONÁRIOS DAS SAÍDAS
% ============================================================

xz_ss = [ ...
    x_ss(:);
    z_ss(:) ...
];

y_ss = xz_ss(idx_y);


%% ============================================================
% SELECIONA C E D
% ============================================================

C_dim = C_all_dim(idx_y,:);

D_dim = D_all_dim(idx_y,idx_u);


%% ============================================================
% MODELO DIMENSIONAL
% ============================================================

sys_dimensional = ss( ...
    A_dim, ...
    B_dim, ...
    C_dim, ...
    D_dim ...
);


%% ============================================================
% TRANSFORMAÇÃO PARA DESVIO PERCENTUAL
%
% x_p = (x - x_ss)/x_ss
%
% u_p = (u - u_ss)/u_ss
%
% y_p = (y - y_ss)/y_ss
%
% Se multiplicarmos tudo por 100, o fator 100 cancela.
%
% Portanto:
%
% A_p = inv(Tx) * A * Tx
%
% B_p = inv(Tx) * B * Tu
%
% C_p = inv(Ty) * C * Tx
%
% D_p = inv(Ty) * D * Tu
% ============================================================

Tx = diag(x_ss);

Tu = diag(u_ss);

Ty = diag(y_ss);


%% ============================================================
% EVITA USAR INV()
% ============================================================

A = Tx \ (A_dim * Tx);

B = Tx \ (B_dim * Tu);

C = Ty \ (C_dim * Tx);

D = Ty \ (D_dim * Tu);


%% ============================================================
% MODELO EM ESPAÇO DE ESTADOS PERCENTUAL
% ============================================================

sys_measured = ss( ...
    A, ...
    B, ...
    C, ...
    D ...
);


%% ============================================================
% NOMES DOS ESTADOS
% ============================================================

sys_measured.StateName = { ...
    'P_man', ...
    'q_tr', ...
    'P_fbhp_1', ...
    'P_choke_1', ...
    'q_mean_1', ...
    'P_fbhp_2', ...
    'P_choke_2', ...
    'q_mean_2' ...
};


%% ============================================================
% NOMES DAS ENTRADAS
% ============================================================

sys_measured.InputName = { ...
    'f_booster', ...
    'f_BCS_1', ...
    'alpha_1', ...
    'f_BCS_2', ...
    'alpha_2' ...
};


%% ============================================================
% NOMES DAS SAÍDAS
% ============================================================

sys_measured.OutputName = { ...
    'P_man', ...
    'P_fbhp_1', ...
    'P_fbhp_2', ...
    'P_intake_1', ...
    'dP_BCS_1', ...
    'P_intake_2', ...
    'dP_BCS_2' ...
};


%% ============================================================
% MOSTRA AS MATRIZES
% ============================================================

fprintf('\n============================================\n');
fprintf(' MATRIZ A - DESVIO PERCENTUAL\n');
fprintf('============================================\n');
disp(A);

fprintf('\n============================================\n');
fprintf(' MATRIZ B - DESVIO PERCENTUAL\n');
fprintf('============================================\n');
disp(B);

fprintf('\n============================================\n');
fprintf(' MATRIZ C - DESVIO PERCENTUAL\n');
fprintf('============================================\n');
disp(C);

fprintf('\n============================================\n');
fprintf(' MATRIZ D - DESVIO PERCENTUAL\n');
fprintf('============================================\n');
disp(D);


%% ============================================================
% DIMENSÕES
% ============================================================

fprintf('\nModelo percentual:\n');
fprintf('A = %d x %d\n',size(A,1),size(A,2));
fprintf('B = %d x %d\n',size(B,1),size(B,2));
fprintf('C = %d x %d\n',size(C,1),size(C,2));
fprintf('D = %d x %d\n',size(D,1),size(D,2));


%% ============================================================
% AUTOVALORES / POLOS DO ESPAÇO DE ESTADOS
% ============================================================

autovalores = eig(A);

fprintf('\n============================================\n');
fprintf(' AUTOVALORES DE A\n');
fprintf('============================================\n');

disp(autovalores);


%% ============================================================
% ESTABILIDADE
% ============================================================

if all(real(autovalores) < 0)

    fprintf('\nSistema linearizado ESTÁVEL.\n');

elseif any(real(autovalores) > 0)

    fprintf('\nSistema linearizado INSTÁVEL.\n');

else

    fprintf('\nSistema possui polos sobre o eixo imaginário.\n');

end


%% ============================================================
% FUNÇÃO DE TRANSFERÊNCIA MIMO
% ============================================================

sys_tf = tf(sys_measured);

fprintf('\n============================================\n');
fprintf(' FUNÇÕES DE TRANSFERÊNCIA\n');
fprintf('============================================\n');

disp(sys_tf);


%% ============================================================
% FORMA POLO-ZERO-GANHO
% ============================================================

sys_zpk = zpk(sys_measured);


%% ============================================================
% POLOS GLOBAIS DO MODELO
% ============================================================

polos = pole(sys_measured);

fprintf('\n============================================\n');
fprintf(' POLOS DO SISTEMA\n');
fprintf('============================================\n');

disp(polos);


%% ============================================================
% MAPA DE POLOS DO MODELO
% ============================================================

figure;

pzmap(sys_measured);

grid on;

title('Polos e Zeros - Modelo Linearizado Percentual');


%% ============================================================
% ESCOLHA DE UMA FUNÇÃO DE TRANSFERÊNCIA
%
% Exemplo:
%
% saída 1  = P_man
% entrada 2 = f_BCS_1
%
% portanto:
%
% G = P_man / f_BCS_1
% ============================================================

saida = 1;
entrada = 2;

G = sys_tf(saida,entrada);


%% ============================================================
% MOSTRA A FT ESCOLHIDA
% ============================================================

fprintf('\n============================================\n');
fprintf(' FT SELECIONADA\n');
fprintf('============================================\n');

fprintf('Saída   : %s\n', ...
    sys_measured.OutputName{saida});

fprintf('Entrada : %s\n', ...
    sys_measured.InputName{entrada});

disp(G);


%% ============================================================
% POLOS E ZEROS DA FT
% ============================================================

p = pole(G);
z = zero(G);

fprintf('\nPolos da FT:\n');
disp(p);

fprintf('Zeros da FT:\n');
disp(z);


%% ============================================================
% PZMAP DA FT INDIVIDUAL
% ============================================================

figure;

pzmap(G);

grid on;

title(sprintf( ...
    '%s / %s', ...
    sys_measured.OutputName{saida}, ...
    sys_measured.InputName{entrada} ...
));


%% ============================================================
% PROCURA POLOS E ZEROS PRÓXIMOS
% ============================================================

fprintf('\n============================================\n');
fprintf(' DISTÂNCIA ENTRE POLOS E ZEROS\n');
fprintf('============================================\n');

if isempty(p)

    fprintf('A FT não possui polos.\n');

elseif isempty(z)

    fprintf('A FT não possui zeros.\n');

else

    distancias = abs(p - z.');

    for ip = 1:length(p)

        for iz = 1:length(z)

            fprintf( ...
                'Polo %2d = % .8e %+.8ej | ', ...
                ip, ...
                real(p(ip)), ...
                imag(p(ip)));

            fprintf( ...
                'Zero %2d = % .8e %+.8ej | ', ...
                iz, ...
                real(z(iz)), ...
                imag(z(iz)));

            fprintf( ...
                'dist = %.8e\n', ...
                distancias(ip,iz));

        end

    end

    [dist_min,idx_min] = min(distancias(:));

    [ip_min,iz_min] = ind2sub( ...
        size(distancias), ...
        idx_min);

    fprintf('\nPar mais próximo:\n');

    fprintf('Polo =\n');
    disp(p(ip_min));

    fprintf('Zero =\n');
    disp(z(iz_min));

    fprintf( ...
        'Distância = %.12e\n', ...
        dist_min);

end


%% ============================================================
% NÃO FAZER CANCELAMENTO AUTOMÁTICO AINDA
%
% Quando definirmos uma tolerância justificável:
%
% tol = 1e-5;
% G_reduzida = minreal(G,tol);
%
% figure;
% pzmap(G_reduzida);
% grid on;
%
% ============================================================


%% ============================================================
% TABELA COM TODAS AS 35 FUNÇÕES DE TRANSFERÊNCIA
% ============================================================

ny = size(sys_measured,1);
nmu = size(sys_measured,2);

fprintf('\n============================================\n');
fprintf(' POLOS/ZEROS DE TODAS AS FTs\n');
fprintf('============================================\n');

for iy = 1:ny

    for iu = 1:nmu

        Gij = sys_tf(iy,iu);

        pij = pole(Gij);
        zij = zero(Gij);

        fprintf('\n--------------------------------------------\n');

        fprintf( ...
            '%s  <--  %s\n', ...
            sys_measured.OutputName{iy}, ...
            sys_measured.InputName{iu});

        fprintf('Polos:\n');
        disp(pij);

        fprintf('Zeros:\n');
        disp(zij);

    end

end