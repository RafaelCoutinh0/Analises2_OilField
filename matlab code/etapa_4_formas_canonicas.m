%% 4. FORMAS CANÔNICAS - equivalente a 4.formas_canônicas.ipynb
% Matriz de transferência MIMO, realizações SISO controlável e observável,
% diagonalização por autovetores (Jordan para o caso diagonalizável).
% Em problemas mal condicionados utiliza-se Schur e relata-se a diferença.
clearvars; close all;
raiz = configurar_projeto();
arquivo = fullfile(raiz,'outputs','modelo_linear.mat');
if ~isfile(arquivo)
    error('Execute antes etapa_2_linearizacao.m para gerar a matriz A.');
end
s = load(arquivo,'lin');lin = s.lin; n = lin.nomes;
sys = ss(lin.A,lin.B,lin.C,lin.D);
G_tf = tf(sys);  % matriz completa: 8 saídas x 5 entradas
ny = size(G_tf,1); nu = size(G_tf,2);

fprintf('\n=========================================================\n');
fprintf(' 4. FORMAS CANÔNICAS\n');
fprintf('=========================================================\n');
fprintf('Matriz de transferência: %d x %d = %d FTs SISO.\n',ny,nu,ny*nu);
save(fullfile(raiz,'outputs','funcoes_transferencia.mat'),'G_tf');

% Como na referência: exportação da matriz de transferência em arquivo TeX.
fid = fopen(fullfile(raiz,'outputs','tf.tex'),'w');
if fid<0, error('Não foi possível criar outputs/tf.tex'); end
for iy = 1:ny
    for iu = 1:nu
        [num0,den0] = tfdata(G_tf(iy,iu),'v');
        fprintf(fid,'G_{%d,%d}(s) = \\frac{%s}{%s}\n\n',iy,iu, ...
            polinomio_latex(num0),polinomio_latex(den0));
    end
end
fclose(fid);

%% Função de transferência escolhida - análoga ao ij da referência.
% Referência original: T3 <- FR, que não existem neste modelo de poços.
% Adaptação: q_mean_1 <- f_BCS_1 (o efeito da bomba sobre a vazão do poço).
i = 5; j = 2;
G = G_tf(i,j);
fprintf('\nG de %s para %s:\n',n.u{j},n.x{i});
disp(G);

% Minimiza cancelamentos de polos/zeros devido ao modelo MIMO reduzido.
% Mantemos o G original e um G_min para obter uma realização SISO robusta.
G_min = minreal(G,1e-5);
fprintf('Ordem antes: %d; após minreal: %d\n', ...
    order(G),order(G_min));
[num,den] = tfdata(G_min,'v');
[Ac,Bc,Cc,Dc] = forma_canonica_controlavel(num,den);
if isempty(Ac)
    error('A FT selecionada é constante; escolha outra para ilustrar estados canônicos.');
end
fprintf('\nForma canônica controlável:\n');
fprintf('A_c =\n');disp(Ac);
fprintf('B_c =\n');disp(Bc);
fprintf('C_c =\n');disp(Cc);
fprintf('D_c =\n');disp(Dc);

%% Forma canônica observável: dualidade, igual ao notebook.
Ao = Ac.'; Bo = Cc.'; Co = Bc.'; Do = Dc;
fprintf('\nForma canônica observável:\n');
fprintf('A_o =\n');disp(Ao);
fprintf('B_o =\n');disp(Bo);
fprintf('C_o =\n');disp(Co);
fprintf('D_o =\n');disp(Do);

%% Jordan quando diagonalizável com condicionamento suficiente.
[V,~] = eig(Ac);
if rcond(V)>1e-10
    Aj = V\(Ac*V); Bj = V\Bc; Cj = Cc*V; Dj = Dc;
    titulo_j = 'Forma de Jordan';
    fprintf('\nJordan: autovetores independentes com condicionamento aceitável.\n');
else
    % Schur unitária evita amplificar erros em autovetores quase dependentes.
    [Q,T] = schur(Ac,'complex');
    Aj = T; Bj = Q'*Bc; Cj = Cc*Q; Dj = Dc;
    titulo_j = 'Forma de Schur (alternativa)';
    fprintf(['\nAVISO: diagonalização por autovetores mal condicionada. ', ...
        'Foi usada a forma de Schur, não a forma de Jordan exata.\n']);
end
fprintf('%s:\n',titulo_j);
fprintf('A_J =\n');disp(Aj);
fprintf('B_J =\n');disp(Bj);
fprintf('C_J =\n');disp(Cj);
fprintf('D_J =\n');disp(Dj);


%% Comparação das 3 realizações em 2 linhas x 3 colunas
% Degrau unitário aplicado em t = 20 s.
% Ajuste automático dos eixos Y com margem de 20%.

polos = eig(Ac);
pestaveis = abs(real(polos(real(polos)<-1e-10)));

if isempty(pestaveis)
    T_final = 200;
else
    T_final = min(500,max(30,5/min(pestaveis)));
end

%% Tempo de aplicação do degrau

t_degrau = 20;

T_final = T_final + t_degrau;

tt = unique([linspace(0,T_final,500),t_degrau]);

%% Definição das três realizações

sistemas = {{'Canônica controlável',Ac,Bc,Cc,Dc}, ...
    {'Canônica observável',Ao,Bo,Co,Do}, ...
    {titulo_j,Aj,Bj,Cj,Dj}};

fig = figure('Color','w','Name','Formas canônicas', ...
    'Position',[50 60 1800 760]);

respostas = cell(1,3);

%% Simulação e representação gráfica

for col = 1:3

    item = sistemas{col};

    % Condições iniciais nulas
    Y = zeros(1,numel(tt));
    X = zeros(size(item{2},1),numel(tt));

    % Aplicar degrau a partir de 20 segundos
    idx = tt >= t_degrau;

    tempo_resposta = tt(idx) - t_degrau;

    [Yp,Xp] = simular_passo_siso( ...
        item{2},item{3},item{4},item{5}, ...
        tempo_resposta);

    % Armazenar respostas após o degrau
    Y(idx) = reshape(Yp,1,[]);
    X(:,idx) = Xp;

    respostas{col} = real(Y);

    %% Gráficos dos estados

    ax1 = subplot(2,3,col);

    plot(tt,real(X.'),'LineWidth',1.15);

    xlabel('Tempo / s');
    ylabel('Amplitude');

    title(['Estados: ' item{1}]);
    grid on;

    if size(X,1)<=12
        legend(arrayfun(@(k)sprintf('x%d',k), ...
            1:size(X,1),'UniformOutput',false), ...
            'Location','best');
    end

    xlim(ax1,[0 T_final]);

    %% Gráficos das saídas

    ax2 = subplot(2,3,3+col);

    plot(tt,real(Y),'LineWidth',1.55);

    xlabel('Tempo / s');

    ylabel(sprintf('\\Delta %s / %s', ...
        n.x_labels{i},n.x_units{i}));

    title(['Saída: ' item{1}]);
    grid on;

    xlim(ax2,[0 T_final]);

end

%% Ajuste automático dos eixos Y
% Margem de 20% acima do máximo e abaixo do mínimo.

eixos = findall(fig,'Type','axes');

for k = 1:numel(eixos)

    ax = eixos(k);

    curvas = findall(ax,'Type','line');

    valores = [];

    % Usar m para não sobrescrever j
    for m = 1:numel(curvas)

        y = curvas(m).YData;

        valores = [valores; real(y(:))];

    end

    valores = valores(isfinite(valores));

    if isempty(valores)
        continue;
    end

    ymin = min(valores);
    ymax = max(valores);

    amplitude = ymax - ymin;

    if amplitude > 0
        margem = 0.20 * amplitude;
    else
        margem = max(0.20 * abs(ymax),0.01);
    end

    ylim(ax,[ymin-margem, ymax+margem]);

    xlim(ax,[0 T_final]);

end

%% Salvar figura com os eixos ajustados

salvar_figura(fig,fullfile(raiz, ...
    'figures','formas_canonicas', ...
    sprintf('formas_canonicas-%s_%s.png', ...
    n.x{i},n.u{j})));

%% Comparação das respostas

erro_co = max(abs(respostas{1}-respostas{2}));

erro_cj = max(abs(respostas{1}-respostas{3}));

fprintf('\nErro máximo das saídas (controlável x observável): %.6e\n', ...
    erro_co);

fprintf('Erro máximo das saídas (controlável x Jordan/Schur): %.6e\n', ...
    erro_cj);

%% Verificação das diferenças numéricas

if max(erro_co,erro_cj)> ...
        1e-3*max(1,max(abs(respostas{1})))

    warning('BCS:CanonicaNumerica', ...
        ['Diferenças acima da tolerância; ', ...
        'verificar condicionamento e cancelamentos.']);

end

%% Salvar resultados

save(fullfile(raiz,'outputs','formas_canonicas.mat'), ...
    'Ac','Bc','Cc','Dc', ...
    'Ao','Bo','Co','Do', ...
    'Aj','Bj','Cj','Dj', ...
    'G','G_min','i','j', ...
    'tt','respostas','titulo_j');

fprintf('Arquivo TeX com 40 FTs em outputs/tf.tex.\n');

fprintf('\n========================================\n');
fprintf(' FORMAS CANÔNICAS CONCLUÍDAS\n');
fprintf('========================================\n');

fprintf('Degrau aplicado em: %.0f s\n',t_degrau);
fprintf('Tempo final: %.2f s\n',T_final);
fprintf('Margem dos eixos Y: 20%%\n');
