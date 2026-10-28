%% 3. ESTABILIDADE - equivalente a 3.estabilidade.ipynb
% Autovalores e mapas de fase 3D de 12 condições iniciais próximas da raiz.
clearvars; close all;
raiz = configurar_projeto();
carregar = fullfile(raiz,'outputs','modelo_linear.mat');
if ~isfile(carregar)
    error('Execute antes etapa_2_linearizacao.m para gerar outputs/modelo_linear.mat');
end
s = load(carregar,'lin'); lin = s.lin;
A = lin.A; n = lin.nomes;
lambdas = eig(A);
tol = 1e-9;
if all(real(lambdas)<-tol)
    classificacao = 'ESTÁVEL';
elseif any(real(lambdas)>tol)
    classificacao = 'INSTÁVEL';
else
    classificacao = 'MARGINAL / REQUER ANÁLISE ADICIONAL';
end
fprintf('\n=========================================================\n');
fprintf(' 3. ESTABILIDADE\n');
fprintf('=========================================================\n');
fprintf('Autovalores de A: \n');disp(lambdas);
fprintf('lambda (partes reais arredondadas): [');
fprintf('%.2f ',real(lambdas));
fprintf('] | %s\n',classificacao);

% Mapa de fase da referência: 12 CIs e 3 painéis 3D.
rng(0);
delta = 0.05*abs(lin.x_ss(:));   % +-5%% do ponto de operação
n_traj = 12;
dX0 = (2*rand(numel(delta),n_traj)-1).*delta;
T_final = 300; tt = linspace(0,T_final,250);
Ad = expm(A*(tt(2)-tt(1)));
traj = zeros(numel(lin.x_ss),numel(tt),n_traj);
for k = 1:n_traj
    traj(:,1,k) = dX0(:,k);
    for p = 2:numel(tt)
        traj(:,p,k) = Ad*traj(:,p-1,k);
    end
end

% 3 projeções no espaço de estados de 8 dimensões.
grupos = {[1 3 6],[2 5 8],[1 4 7]};
titulos = {'Pressão manifold e FBHP','Vazões do campo', ...
    'Pressão manifold e chokes'};
fig = figure('Color','w','Name','Mapa de fase', ...
    'Position',[40 80 1800 650]);
for g = 1:3
    idx = grupos{g};
    subplot(1,3,g);hold on;
    for k = 1:n_traj
        plot3(squeeze(traj(idx(1),:,k)), ...
            squeeze(traj(idx(2),:,k)), ...
            squeeze(traj(idx(3),:,k)),'LineWidth',1.0);
    end
    plot3(dX0(idx(1),:),dX0(idx(2),:),dX0(idx(3),:), ...
        'k.','MarkerSize',10);
    plot3(0,0,0,'ro','MarkerFaceColor','r','MarkerSize',7);
    xlabel(sprintf('\Delta%s / %s',n.x_labels{idx(1)},n.x_units{idx(1)}));
    ylabel(sprintf('\Delta%s / %s',n.x_labels{idx(2)},n.x_units{idx(2)}));
    zlabel(sprintf('\Delta%s / %s',n.x_labels{idx(3)},n.x_units{idx(3)}));
    title(titulos{g});grid on;view(42,25);
end
salvar_figura(fig,fullfile(raiz,'figures','estabilidade','mapa_de_fase.png'));
fprintf('Mapas de fase salvos em figures/estabilidade/mapa_de_fase.png\n');
