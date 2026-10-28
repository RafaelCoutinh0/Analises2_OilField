%% 2. LINEARIZACAO
% Linearizacao do modelo DAE com eliminacao das variaveis algebricas.
% Representacao principal em variaveis de desvio dimensional.

clearvars; close all;

raiz = configurar_projeto();
mdl = sym_dynamics();
lin = linearizar_dae(mdl);

A = lin.A;
B = lin.B;
C = lin.C;
D = lin.D;
n = lin.nomes;

fprintf('\n=========================================================\n');
fprintf(' 2. LINEARIZACAO DO MODELO NAO LINEAR DAE\n');
fprintf('=========================================================\n');
fprintf('Ponto de operacao (regime estacionario):\n');

for k = 1:numel(n.x)
    fprintf('  %-13s = %14.6f %s\n', ...
        n.x{k},lin.x_ss(k),n.x_units{k});
end

fprintf('Entradas no ponto de operacao:\n');
for k = 1:numel(n.u)
    fprintf('  %-13s = %14.6f %s\n', ...
        n.u{k},lin.u_ss(k),n.u_units{k});
end

fprintf('\nMatrizes do modelo linearizado - DESVIO DIMENSIONAL:\n');
fprintf('A = (%d x %d):\n',size(A)); disp(A);
fprintf('B = (%d x %d):\n',size(B)); disp(B);
fprintf('C = (%d x %d):\n',size(C)); disp(C);
fprintf('D = (%d x %d):\n',size(D)); disp(D);
fprintf('A_pct:\n'); disp(lin.A_pct);
fprintf('B_pct:\n'); disp(lin.B_pct);
fprintf('Jacobianas: max|F_z|=%.6e; max|G_u|=%.6e\n', ...
    max(abs(lin.f_z(:))),max(abs(lin.g_u(:))));
fprintf('rcond(G_z)=%.6e\n',rcond(lin.Gz));

%% Exportacao das matrizes
save(fullfile(raiz,'outputs','modelo_linear.mat'), ...
    'lin','A','B','C','D');
save(fullfile(raiz,'outputs','A_B_C_D.mat'), ...
    'A','B','C','D');
writematrix(A,fullfile(raiz,'outputs','A.csv'));
writematrix(B,fullfile(raiz,'outputs','B.csv'));
writematrix(C,fullfile(raiz,'outputs','C.csv'));

%% RESPOSTAS DO SISTEMA LINEARIZADO
% Um unico pulso por entrada, entre 20 e 420 s.
% Todas as figuras usam 20% da amplitude como margem acima e abaixo.

T_final = 800;
t = linspace(0,T_final,801);
t1 = 20;
t2 = 420;
amplitudes = [2, 2, 0.03, 2, 0.03];
percentual_margem = 0.20;

for k = 1:numel(n.idx_u)
    U = repmat(mdl.u0(:),1,numel(t));
    idx = n.idx_u(k);
    u_nominal = mdl.u0(idx);

    U(idx,t >= t1 & t < t2) = u_nominal + amplitudes(k);
    U(idx,t >= t2) = u_nominal;

    [X,Z] = simular_linear(lin,t,U);
    fprintf('\nPerturbacao em %s\n',n.u{k});

    %% Figura das entradas
    fig = figure('Color','w','Name','Entradas', ...
        'Position',[50 100 1700 350]);
    for j = 1:numel(n.idx_u)
        ax = subplot(1,numel(n.idx_u),j);
        idxj = n.idx_u(j);
        valores = U(idxj,:);
        stairs(ax,t,valores,'LineWidth',1.7);
        title(ax,n.u{j},'Interpreter','none');
        xlabel(ax,'Tempo / s');
        ylabel(ax,sprintf('%s / %s',n.u_labels{j},n.u_units{j}), ...
            'Interpreter','tex');
        xlim(ax,[0 T_final]); grid(ax,'on');

        validos = valores(isfinite(valores));
        if ~isempty(validos)
            ymin = min(validos); ymax = max(validos);
            faixa = ymax-ymin;
            if faixa > 0
                margem = percentual_margem*faixa;
            else
                margem = max(0.02*abs(ymax),0.02);
            end
            ylim(ax,[ymin-margem ymax+margem]);
        end
    end
    salvar_figura(fig,fullfile(raiz,'figures','linear', ...
        sprintf('pulse-%s-inputs.png',n.u{k})));

    %% Figura dos 8 estados (2 linhas x 4 colunas)
    fig = figure('Color','w','Name','Estados', ...
        'Position',[40 50 1700 680]);
    for j = 1:numel(n.x)
        ax = subplot(2,4,j);
        valores = X(j,:);
        plot(ax,t,valores,'LineWidth',1.6);
        title(ax,n.x{j},'Interpreter','none');
        xlabel(ax,'Tempo / s');
        ylabel(ax,sprintf('%s / %s',n.x_labels{j},n.x_units{j}), ...
            'Interpreter','tex');
        xlim(ax,[0 T_final]); grid(ax,'on');

        validos = valores(isfinite(valores));
        if ~isempty(validos)
            ymin = min(validos); ymax = max(validos);
            faixa = ymax-ymin;
            if faixa > 0
                margem = percentual_margem*faixa;
            else
                margem = max(0.02*abs(ymax),0.02);
            end
            ylim(ax,[ymin-margem ymax+margem]);
        end
    end
    salvar_figura(fig,fullfile(raiz,'figures','linear', ...
        sprintf('pulse-%s-outputs.png',n.u{k})));

    %% Figura das 4 variaveis algebricas (2 linhas x 2 colunas)
    fig = figure('Color','w','Name','Variaveis algebricas', ...
        'Position',[80 80 1240 650]);
    for j = 1:numel(n.z)
        ax = subplot(2,2,j);
        valores = Z(j,:);
        plot(ax,t,valores,'LineWidth',1.6);
        title(ax,n.z{j},'Interpreter','none');
        xlabel(ax,'Tempo / s');
        ylabel(ax,sprintf('%s / %s',n.z_labels{j},n.z_units{j}), ...
            'Interpreter','tex');
        xlim(ax,[0 T_final]); grid(ax,'on');

        validos = valores(isfinite(valores));
        if ~isempty(validos)
            ymin = min(validos); ymax = max(validos);
            faixa = ymax-ymin;
            if faixa > 0
                margem = percentual_margem*faixa;
            else
                margem = max(0.02*abs(ymax),0.02);
            end
            ylim(ax,[ymin-margem ymax+margem]);
        end
    end
    salvar_figura(fig,fullfile(raiz,'figures','linear', ...
        sprintf('pulse-%s-algebraics.png',n.u{k})));

    fprintf('Figuras salvas para %s\n',n.u{k});
end

fprintf('\n========================================\n');
fprintf(' SIMULACOES CONCLUIDAS\n');
fprintf('========================================\n');
fprintf('Inicio da perturbacao: %.0f s\n',t1);
fprintf('Retorno ao nominal: %.0f s\n',t2);
fprintf('Tempo final: %.0f s\n',T_final);
fprintf('Margem vertical por eixo: %.0f%% da amplitude\n',100*percentual_margem);

%% Comparação linear x não linear
% Saída: incremento de pressão da BCS do poço 1.
% Entrada: frequência da BCS do poço 1.

tc = linspace(0,150,751);

entrada = 2;
idx_entrada = n.idx_u(entrada);

% DeltaP_BCS_1: segunda variável algébrica
idx_dp = 2;

% Entradas inicialmente nos valores nominais
Uc = repmat(mdl.u0(:),1,numel(tc));

f0 = mdl.u0(idx_entrada);

%% Aplicação dos degraus

% Perturbação pequena: 52 Hz
Uc(idx_entrada,tc >= 20 & tc < 70) = 55;

% Perturbação grande: 65 Hz
Uc(idx_entrada,tc >= 70 & tc < 110) = 65;

% Retorno ao ponto nominal
Uc(idx_entrada,tc >= 110) = f0;

%% Simulação dos modelos

[~,Z_nl] = simular_dae(mdl,tc,Uc);
[~,Z_lin] = simular_linear(lin,tc,Uc);

%% Figura comparativa

fig = figure('Color','w', ...
    'Name','Comparação linear x não linear - DeltaP BCS 1', ...
    'Position',[100 80 960 720]);

% Frequência da BCS 1
ax1 = subplot(2,1,1);

stairs(tc,Uc(idx_entrada,:), ...
    'LineWidth',1.8);

ylabel('f_{BCS,1} / Hz');
xlim([0 150]);
grid on;

% Incremento de pressão da BCS 1
ax2 = subplot(2,1,2);

plot(tc,Z_nl(idx_dp,:), ...
    'LineWidth',1.7);
hold on;

plot(tc,Z_lin(idx_dp,:),'--', ...
    'LineWidth',1.7);

ylabel('\Delta p_{BCS,1} / bar');
xlabel('Tempo / s');

legend('Não linear','Linear','Location','best');

xlim([0 150]);
grid on;

%% Ajuste automático dos eixos

ajustar_eixos(fig,0.15);

% Escala da frequência com margem
ylim(ax1,[45 70]);
yticks(ax1,45:5:70);

xlim(ax1,[0 150]);
xlim(ax2,[0 150]);

%% Salvar figura

arquivo_figura = fullfile(raiz, ...
    'figures','comparacao', ...
    'linear_vs_nonlinear_DeltaP_BCS1.png');

salvar_figura(fig,arquivo_figura);

%% Avaliação do erro

erro_dp = abs(Z_nl(idx_dp,:) - Z_lin(idx_dp,:));

idx_local = tc >= 20 & tc < 70;
idx_grande = tc >= 70 & tc < 110;

fprintf('\n========================================\n');
fprintf(' COMPARAÇÃO LINEAR X NÃO LINEAR\n');
fprintf('========================================\n');

fprintf('Erro máximo global: %.6f bar\n', ...
    max(erro_dp));

fprintf('Erro máximo no intervalo de 52 Hz: %.6f bar\n', ...
    max(erro_dp(idx_local)));

fprintf('Erro máximo no intervalo de 65 Hz: %.6f bar\n', ...
    max(erro_dp(idx_grande)));

fprintf('Figura salva em:\n%s\n',arquivo_figura);

