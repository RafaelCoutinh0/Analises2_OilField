# Relatório LaTeX — análise do modelo BCS

- `main.tex`: texto completo, estruturado nas 5 etapas do MATLAB.
- `figuras/`: PNG originais selecionados do `figures.zip` enviado.
- `relatorio.pdf`: PDF compilado, se incluído.

## Compilação

No Overleaf, importe este ZIP como projeto e abra `main.tex` com compilador pdfLaTeX. Localmente execute `pdflatex main.tex` duas vezes.

## Notas importantes de rastreabilidade

- Modelo MATLAB utilizado: **2 poços**, 8 estados diferenciais, 4 variáveis algébricas, 6 entradas (5 manipuladas com `p_topside` fixa).
- Os coeficientes numéricos de matrizes e os autovalores foram reconstituídos a partir das equações/constantes do `sym_dynamics.m`, reproduzindo o ponto de equilíbrio do MATLAB (não foram fornecidos arquivos de saída `formas_canonicas.mat` ou `discretizacao.mat`). Conferir com `outputs/formas_canonicas.mat` em caso de modificações locais.
- As figuras não lineares entregues no `figures.zip` usam aproximadamente 600 s, diferentemente do padrão original de 400 s no arquivo `etapa_1_sistema.m` acessível aqui.
- O valor de Ts no relatório vem do código enviado, que usa Ts = 80% do limite de estabilidade de Euler; se Ts foi modificado após a geração do pacote, substituir o valor no texto pelo valor efetivamente executado no MATLAB.
- Textos citam os artigos enviados para o contexto físico e padronização da nomenclatura, mas não atribuem execução de MPC, otimização econômica, treinamento de RNA ou modelos incrementais aumentados a estes scripts MATLAB.
