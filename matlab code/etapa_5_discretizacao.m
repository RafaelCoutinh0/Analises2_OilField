%% 5. DISCRETIZAÇÃO - equivalente a 5.discretização.ipynb
% Discretização pelo método de Euler e retentor de ordem zero (ZOH).
% Ts é ESCOLHIDO dos polos de A para que Euler seja numericamente estável;
% copiar Ts=0.02 h do reator não é adequado para esta DAE em segundos.
clearvars; close all;
raiz = configurar_projeto();
arquivo = fullfile(raiz,'outputs','modelo_linear.mat');
if ~isfile(arquivo)
    error('Execute antes etapa_2_linearizacao.m para gerar a matriz A.');
end
s = load(arquivo,'lin');lin = s.lin; n = lin.nomes;
A = lin.A; B = lin.B; C = lin.C; D = lin.D;
nx = size(A,1); nu = size(B,2);
j = 2; % entrada com degrau unitário: f_BCS_1, análogo a FR do notebook
u = zeros(nu,1); u(j) = 1.0;

lambdas = eig(A);
neg = lambdas(real(lambdas)<-1e-10);
if numel(neg) ~= nx
    warning('BCS:Euler', ...
        'Há polos não estáveis: Euler pode divergir independentemente de Ts.');
    Ts = min(0.005,1/max(abs(lambdas)));
else
    T_limite_Euler = min(-2*real(neg)./(abs(neg).^2));
    Ts = 0.80*T_limite_Euler; % 20%% abaixo do limite de estabilidade
end
T_final = 240; % segundos, cobre os modos lentos ~50 segundos
td = 0:Ts:T_final;  % grade uniforme; último instante <= T_final
N = numel(td);
tc = linspace(0,td(end),1200);

% Euler explícito
Ad_euler = eye(nx)+Ts*A;
Bd_euler = Ts*B;
Cd_euler = C; Dd_euler = D;

% ZOH exato: c2d (equivalente ao cont2discrete(...,method='zoh')).
sys_zoh = c2d(ss(A,B,C,D),Ts,'zoh');
Ad_zoh = sys_zoh.A; Bd_zoh = sys_zoh.B;
Cd_zoh = sys_zoh.C; Dd_zoh = sys_zoh.D;

% Contínuo: resposta à entrada constante em desvio desde t=0.
[Yc,tc] = lsim(ss(A,B,C,D),repmat(u.',numel(tc),1),tc); %#ok<ASGLU>
Xc_desv = Yc.';                 % pois C=I e D=0
Xc_abs  = Xc_desv + lin.x_ss;
Z1c = lin.z_ss(1) + lin.Zx(1,:)*Xc_desv + lin.Zu(1,:)*u;
Yc9 = [Xc_abs;Z1c];

% Modelos discretos: x[k+1]=Ad*x[k]+Bd*u; y[k]=Cd*x[k]+Dd*u.
Yd9 = cell(1,2);Xd = cell(1,2);
sistemas = {{'Euler',Ad_euler,Bd_euler,Cd_euler,Dd_euler}, ...
    {'ZOH',Ad_zoh,Bd_zoh,Cd_zoh,Dd_zoh}};
for modo = 1:2
    it = sistemas{modo};
    x = zeros(nx,N);
    for k = 1:N-1
        x(:,k+1) = it{2}*x(:,k) + it{3}*u;
    end
    y = it{4}*x + it{5}*u;       % desvios diferenciais
    z1 = lin.z_ss(1) + lin.Zx(1,:)*x + lin.Zu(1,:)*u;
    Yd9{modo} = [y+lin.x_ss;z1];
    Xd{modo} = x;
end

fprintf('\n=========================================================\n');
fprintf(' 5. DISCRETIZAÇÃO\n');
fprintf('=========================================================\n');
fprintf('Entrada perturbada: %s (+1 %s)\n',n.u{j},n.u_units{j});
fprintf('Ts = %.8f s (Euler com margem de estabilidade)\n',Ts);
fprintf('Horizonte = %.3f s, %d amostras\n',td(end),N);
fprintf('Polos contínuos (lambda):\n');disp(lambdas);
fprintf('Autovalores discretos Euler:\n');disp(eig(Ad_euler));
fprintf('Autovalores discretos ZOH:\n');disp(eig(Ad_zoh));
fprintf('Maior |polo Euler| = %.8f; maior |polo ZOH| = %.8f\n', ...
    max(abs(eig(Ad_euler))),max(abs(eig(Ad_zoh))));
fprintf('A_d (Euler):\n');disp(Ad_euler);
fprintf('B_d (Euler):\n');disp(Bd_euler);
fprintf('A_d (ZOH):\n');disp(Ad_zoh);
fprintf('B_d (ZOH):\n');disp(Bd_zoh);

% Disposição 3 x 3 semelhante ao notebook (9 grandezas representativas).
% O 9º sinal é P_intake_1 (algébrico), pois há 8 estados diferenciais.
ids = [1 3 6; 2 5 8; 4 7 9];
rotulos = {n.x_labels{:},n.z_labels{1}}; %#ok<CCAT>
unidades = {n.x_units{:},n.z_units{1}};
fig = figure('Color','w','Name','Discretização Euler x ZOH', ...
    'Position',[50 50 1650 970]);
% Mostra amostras iniciais densas e decima o restante para não saturar a figura.
sel = unique([1:min(N,ceil(1/Ts)),round(linspace(1,N,min(2000,N)))]);
for r = 1:3
    for c = 1:3
        idx = ids(r,c);
        subplot(3,3,3*(r-1)+c);hold on;
        plot(tc,Yc9(idx,:),'LineWidth',1.5);
        stairs(td(sel),Yd9{1}(idx,sel),'LineWidth',1.1);
        stairs(td(sel),Yd9{2}(idx,sel),'LineWidth',1.1);
        title(rotulos{idx});
        xlabel('Tempo / s');
        ylabel(sprintf('%s / %s',rotulos{idx},unidades{idx}));
        grid on;
        if r==1 && c==1
            legend('Contínuo','Euler','ZOH','Location','best');
        end
    end
end
salvar_figura(fig,fullfile(raiz,'figures','discretizacao', ...
    sprintf('discretizacao-%s.png',n.u{j})));

% Compara ZOH nos instantes de amostragem contra a solução analítica exata.
% x(Ts) = expm([A B;0 0]*Ts) * [0;u] -> Bd*u.
M = expm([A B;zeros(nu,nx+nu)]*Ts);
erro_ZOH = norm(Ad_zoh-M(1:nx,1:nx),inf) + ...
    norm(Bd_zoh-M(1:nx,nx+1:end),inf);
fprintf('Erro da equivalência ZOH vs expm por passo = %.5e\n',erro_ZOH);
save(fullfile(raiz,'outputs','discretizacao.mat'), ...
    'Ts','T_final','Ad_euler','Bd_euler','Cd_euler','Dd_euler', ...
    'Ad_zoh','Bd_zoh','Cd_zoh','Dd_zoh','j','u','td','Xd');
fprintf('Figura salva em figures/discretizacao/.\n');
