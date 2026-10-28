function mdl = sym_dynamics()
% ============================================================
% SYM_DYNAMICS
%
% Modelo não linear DAE do campo de produção com 2 poços + BCS
%
% Saídas:
%   mdl.x          estados diferenciais simbólicos
%   mdl.z          variáveis algébricas simbólicas
%   mdl.u          entradas simbólicas
%   mdl.f          equações diferenciais
%   mdl.g          equações algébricas
%   mdl.dae        estrutura DAE para CasADi
%   mdl.mani_model modelo completo
%   mdl.x_ss       estado estacionário
%   mdl.z_ss       algébricas no estado estacionário
%   mdl.u0         entradas no ponto de linearização
%
% Requer:
%   Pipe.m
%   Pump.m
%   Choke.m
%   Well.m
%   Manifold.m
%   CasADi
% ============================================================

%% ============================================================
% CASADI
% ============================================================

addpath('C:\Users\maxpi\Documents\MATLAB\casadi-3.8.1-windows64-matlab2018b');

import casadi.*


%% ============================================================
% TUBULAÇÃO DO MANIFOLD
% ============================================================

pipe_mani = Pipe( ...
    0.030 * 2, ...
    500, ...
    0, ...
    8.3022e6, ...
    984, ...
    4, ...
    5.752218216772682e6, ...
    3.903249155428134e7 ...
);


%% ============================================================
% FUNÇÕES DA BCS E CHOKE
% ============================================================

f_ca = MX.sym('f');
q_ca = MX.sym('q');
alpha_ca = MX.sym('alpha');

[H_head, eff_fun, pot_fun] = bcs_functions(f_ca,q_ca);

% Curva das BCS dos poços
head_fun = Function( ...
    'head', ...
    {f_ca,q_ca}, ...
    {64 * H_head});

efficiency_fun = Function( ...
    'efficiency', ...
    {f_ca,q_ca}, ...
    {eff_fun});

power_fun = Function( ...
    'power', ...
    {f_ca,q_ca}, ...
    {pot_fun});

% Booster
booster_fun = Function( ...
    'booster', ...
    {f_ca,q_ca}, ...
    {550 * (f_ca/50)^2});

% Choke
valve_fun = Function( ...
    'choke', ...
    {alpha_ca}, ...
    {choke_fun(alpha_ca)});


%% ============================================================
% BOMBAS
% ============================================================

bcs1 = Pump( ...
    head_fun, ...
    efficiency_fun, ...
    power_fun);

bcs2 = Pump( ...
    head_fun, ...
    efficiency_fun, ...
    power_fun);

booster = Pump( ...
    booster_fun, ...
    efficiency_fun, ...
    power_fun);


%% ============================================================
% CHOKES
% ============================================================

K_choke = 212.54 / (6e4 * sqrt(1e5));

choke1 = Choke(K_choke,valve_fun);
choke2 = Choke(K_choke,valve_fun);


%% ============================================================
% POÇO 1
% ============================================================

pipe_sec1_1 = Pipe( ...
    0.090 * 2, ...
    2.80e3, ...
    1000 - 920, ...
    1.5e9, ...
    984, ...
    0.3, ...
    5.752218216772682e6, ...
    3.903249155428134e7 ...
);

pipe_sec2_1 = Pipe( ...
    0.065 * 2, ...
    8.00e2, ...
    920 - 150, ...
    1.5e9, ...
    984, ...
    4, ...
    7.455247950618545e6, ...
    6.264914244217266e7 ...
);


%% ============================================================
% POÇO 2
% ============================================================

pipe_sec1_2 = Pipe( ...
    0.085 * 2, ...
    3.10e3, ...
    1050 - 920, ...
    1.5e9, ...
    984, ...
    0.3, ...
    5.752218216772682e6, ...
    3.903249155428134e7 ...
);

pipe_sec2_2 = Pipe( ...
    0.062 * 2, ...
    9.00e2, ...
    920 - 170, ...
    1.5e9, ...
    984, ...
    4, ...
    7.455247950618545e6, ...
    6.264914244217266e7 ...
);


%% ============================================================
% RESERVATÓRIOS
% ============================================================

PI_1   = 2.0e-9;
PRES_1 = 35.0e6;

PI_2   = 2.5e-9;
PRES_2 = 30.0e6;


%% ============================================================
% POÇOS
% ============================================================

well1 = Well( ...
    pipe_sec1_1, ...
    pipe_sec2_1, ...
    bcs1, ...
    choke1, ...
    PI_1, ...
    PRES_1 ...
);

well2 = Well( ...
    pipe_sec1_2, ...
    pipe_sec2_2, ...
    bcs2, ...
    choke2, ...
    PI_2, ...
    PRES_2 ...
);


%% ============================================================
% MANIFOLD
% ============================================================

mani = Manifold( ...
    pipe_mani, ...
    booster, ...
    0, ...
    0, ...
    {well1,well2} ...
);


%% ============================================================
% VARIÁVEIS SIMBÓLICAS
% ============================================================

% ------------------------------------------------------------
% Entradas
%
% u =
% 1 f_booster
% 2 p_topside
% 3 f_BCS1
% 4 alpha1
% 5 f_BCS2
% 6 alpha2
% ------------------------------------------------------------

f_BP       = MX.sym('f_BP');
p_topside  = MX.sym('p_topside');

f_ESP_1 = MX.sym('f_ESP_1');
alpha_1 = MX.sym('alpha_1');

f_ESP_2 = MX.sym('f_ESP_2');
alpha_2 = MX.sym('alpha_2');

u = [ ...
    f_BP;
    p_topside;
    f_ESP_1;
    alpha_1;
    f_ESP_2;
    alpha_2 ...
];


%% ============================================================
% ESTADOS
%
% x =
% 1 p_man
% 2 q_tr
% 3 P_fbhp_1
% 4 P_choke_1
% 5 q_mean_1
% 6 P_fbhp_2
% 7 P_choke_2
% 8 q_mean_2
% ============================================================

p_man = MX.sym('p_man');
q_tr  = MX.sym('q_tr');

P_fbhp_1  = MX.sym('P_fbhp_1');
P_choke_1 = MX.sym('P_choke_1');
q_mean_1  = MX.sym('q_mean_1');

P_fbhp_2  = MX.sym('P_fbhp_2');
P_choke_2 = MX.sym('P_choke_2');
q_mean_2  = MX.sym('q_mean_2');

x = [ ...
    p_man;
    q_tr;
    P_fbhp_1;
    P_choke_1;
    q_mean_1;
    P_fbhp_2;
    P_choke_2;
    q_mean_2 ...
];


%% ============================================================
% VARIÁVEIS ALGÉBRICAS
%
% z =
% 1 P_intake_1
% 2 dP_bcs_1
% 3 P_intake_2
% 4 dP_bcs_2
% ============================================================

P_intake_1 = MX.sym('P_intake_1');
DP_bcs_1   = MX.sym('dP_bcs_1');

P_intake_2 = MX.sym('P_intake_2');
DP_bcs_2   = MX.sym('dP_bcs_2');

z = [ ...
    P_intake_1;
    DP_bcs_1;
    P_intake_2;
    DP_bcs_2 ...
];


%% ============================================================
% MODELO SIMBÓLICO
% ============================================================

mani_model = mani.model(0,x,z,u);

nx = length(x);
nz = length(z);
nu = length(u);

% 8 primeiras equações = EDO
f = vertcat(mani_model{1:nx});

% 4 últimas equações = algébricas
g = vertcat(mani_model{nx+1:nx+nz});


%% ============================================================
% ENTRADA DO PONTO DE LINEARIZAÇÃO
% ============================================================

u0 = [ ...
    45;
    1e5;
    50;
    0.8;
    50;
    0.8 ...
];


%% ============================================================
% CHUTE INICIAL
% ============================================================

x0 = [ ...
    76.52500;
    4*85;
    64.11666;
    120.91641;
    85;
    64.11666;
    120.91641;
    85 ...
];

z0 = [ ...
    30.03625;
    239.95338 - 30.03625;
    30.03625;
    239.95338 - 30.03625 ...
];


%% ============================================================
% ESTADO ESTACIONÁRIO
% ============================================================

mani_solver = @(y) mani_solver_eval(y,mani,u0);

options = optimoptions( ...
    'fsolve', ...
    'Display','off');

[y_ss,fval,exitflag,output] = fsolve( ...
    mani_solver, ...
    [x0;z0], ...
    options);

x_ss = y_ss(1:nx);
z_ss = y_ss(nx+1:nx+nz);


%% ============================================================
% DAE CASADI
% ============================================================

dae = struct( ...
    'x',x, ...
    'z',z, ...
    'p',u, ...
    'ode',f, ...
    'alg',g ...
);


%% ============================================================
% ORGANIZA RESULTADOS
% ============================================================

mdl = struct;

mdl.mani = mani;

mdl.x = x;
mdl.z = z;
mdl.u = u;

mdl.f = f;
mdl.g = g;

mdl.mani_model = mani_model;
mdl.dae = dae;

mdl.nx = nx;
mdl.nz = nz;
mdl.nu = nu;

mdl.x_ss = x_ss;
mdl.z_ss = z_ss;
mdl.u0 = u0;

mdl.fsolve_residual = fval;
mdl.fsolve_exitflag = exitflag;
mdl.fsolve_output = output;


%% ============================================================
% INFORMAÇÕES
% ============================================================

fprintf('\n============================================\n');
fprintf(' MODELO SIMBÓLICO CRIADO\n');
fprintf('============================================\n');

fprintf('Estados diferenciais : %d\n',nx);
fprintf('Variáveis algébricas : %d\n',nz);
fprintf('Entradas              : %d\n',nu);

fprintf('\nEstado estacionário x_ss:\n');
disp(x_ss);

fprintf('Estado estacionário z_ss:\n');
disp(z_ss);

fprintf('Máximo resíduo estacionário = %.12e\n', ...
    max(abs(fval)));

end


%% ============================================================
% FUNÇÃO DA BCS
% ============================================================

function [H,eff,pot] = bcs_functions(f,q)

f0 = 50;

q0 = q * (f0/f);

H0 = ...
    -29.5615 * (q0/0.0353)^4 ...
    +25.3722 * (q0/0.0353)^3 ...
    -8.7944 * (q0/0.0353)^2 ...
    -8.5706 * (q0/0.0353) ...
    +21.4278;

H = H0 * (f/f0)^2;

eff = 1 + 0*f;
pot = 1 + 0*f;

end


%% ============================================================
% CARACTERÍSTICA DA CHOKE
% ============================================================

function value = choke_fun(alpha)

import casadi.*

zc = [ ...
    0,13,15,17,19,22,25,29,30,32,34,36,39,41,44, ...
    46,49,52,55,58,61,64,67,71,75,78,82,86,91,96,100.01 ...
];

G = [ ...
    0, ...
    0.011052632, ...
    0.024210526, ...
    0.043157895, ...
    0.067894737, ...
    0.097894737, ...
    0.133157895, ...
    0.173684211, ...
    0.195789474, ...
    0.219473684, ...
    0.244736842, ...
    0.271052632, ...
    0.298947368, ...
    0.328421053, ...
    0.358947368, ...
    0.390526316, ...
    0.423684211, ...
    0.458421053, ...
    0.494210526, ...
    0.531578947, ...
    0.570526316, ...
    0.610526316, ...
    0.651578947, ...
    0.694210526, ...
    0.738421053, ...
    0.784210526, ...
    0.830526316, ...
    0.878947368, ...
    0.928421053, ...
    0.979473684, ...
    1 ...
];

fun_choke = interpolant( ...
    'choke_curve', ...
    'bspline', ...
    {zc}, ...
    G);

value = fun_choke(alpha*100);

end


%% ============================================================
% FUNÇÃO PARA O FSOLVE
% ============================================================

function Fval = mani_solver_eval(y,mani,u0)

import casadi.*

nx = 8;

x_num = y(1:nx);
z_num = y(nx+1:end);

model_values = mani.model( ...
    0, ...
    x_num, ...
    z_num, ...
    u0);

Fcas = vertcat(model_values{:});

Fval = full(Fcas);

end