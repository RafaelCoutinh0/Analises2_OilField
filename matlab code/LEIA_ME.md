# Análise dinâmica de um campo de produção com BCS — MATLAB

## Objetivo

Reprodução **em MATLAB** das cinco etapas dos notebooks Jupyter fornecidos como referência, **substituindo o reator químico da referência pelo modelo DAE de dois poços com BCS fornecido pelo usuário**. O modelo físico é o de `sym_dynamics.m`, que utiliza `Manifold.m`, `Well.m`, `Pipe.m`, `Pump.m` e `Choke.m`.

| Referência Jupyter | Reprodução MATLAB | Análises |
|---|---|---|
| `1.sistema.ipynb` | `etapa_1_sistema.m` | Simulação não linear, uma perturbação por entrada, figuras de entradas e saídas |
| `2.linearização.ipynb` | `etapa_2_linearizacao.m` | Jacobianas, linearização DAE, respostas lineares, comparação com o não linear, exportação das matrizes |
| `3.estabilidade.ipynb` | `etapa_3_estabilidade.m` | Autovalores, classificação de estabilidade, três mapas de fase 3D com 12 condições iniciais |
| `4.formas_canônicas.ipynb` | `etapa_4_formas_canonicas.m` | Matriz completa de FTs, canônicas controlável e observável, Jordan/Schur, resposta ao degrau |
| `5.discretização.ipynb` | `etapa_5_discretizacao.m` | Euler x ZOH x contínuo, figura com 3 linhas e 3 colunas |

### Diferenças necessárias entre os sistemas

- Os notebooks modelam um reator de **9 estados, 3 entradas, tempo em horas**. Seu modelo é uma **DAE com 8 estados diferenciais, 4 variáveis algébricas, 6 entradas totais, tempo em segundos**.
- As 5 entradas analisadas são `f_BP`, `f_BCS_1`, `alpha_1`, `f_BCS_2`, `alpha_2`. `p_topside` é mantida no valor nominal, conforme o `linearizacao_ss_ft.m` original. Assim, a matriz de transferência tem **8 saídas diferenciais × 5 entradas = 40 FTs**.
- O Jupyter original usa `C=I_9`. Aqui a versão equivalente é `C=I_8`; estados algébricos são tratados **separadamente** por `dz=Zx*dx+Zu*du`. Adicionalmente, a discretização mostra `P_intake_1` como nono painel, explicitamente identificado como variável algébrica.
- A referência usa perturbações definidas na pasta `lib` do projeto Python, que **não foi anexada**. Por isso, os pulsos são configurados de forma transparente em `lib/perfil_entradas.m`; não alegamos reproduzir os mesmos pulsos numéricos do reator.
- As matrizes principais `A,B,C,D` usam **desvio dimensional**, compatível com os notebooks: `δx = x − x_ss`, `δu = u − u_ss`. O seu script `linearizacao_ss_ft.m` usa desvio relativo/percentual; para preservá-lo, também são exportadas `A_pct,B_pct,C_pct,D_pct` com `lin.A_pct`, etc.
- **Euler deve usar `Ts` pequeno** devido aos polos rápidos do campo BCS. O script calcula automaticamente a maior amostragem estável pelo espectro de `A` e usa 80% desse limite; não copia literalmente `0.02 h` da referência.
- A realização numérica de Jordan só é usada se a matriz de autovetores for adequadamente condicionada. Caso contrário, o script avisa e usa a forma de Schur (não a apresenta incorretamente como Jordan).

## Dependências

1. MATLAB (preferencialmente versão R2020a ou superior).
2. CasADi para MATLAB (recomendado o mesmo pacote já usado no seu projeto, série **3.8.x**, com suporte a IDAS).
3. **Optimization Toolbox** (`fsolve`), pois o seu `sym_dynamics.m` calcula o ponto estacionário com `fsolve`.
4. **Control System Toolbox** (`ss`, `tf`, `minreal`, `c2d`, `lsim`).

Não é necessário Jupyter, Python ou Symbolic Math Toolbox. As figuras são geradas no MATLAB e exportadas como PNG. Não há dependência de caminho absoluto de outro computador.

## Como executar

1. Extraia o ZIP e abra a pasta `Analise_BCS_MATLAB` no MATLAB, mantendo `lib/` e `modelo/` na mesma estrutura.
2. Configure o CasADi no MATLAB **antes** da execução, usando a sua pasta de instalação, por exemplo:

   ```matlab
   addpath('C:\caminho\para\casadi-3.8.1-windows64-matlab2018b')
   ```

   Alternativamente, defina a variável de ambiente `CASADI_MATLAB_PATH` apontando para a mesma pasta (antes de abrir o MATLAB). A função `configurar_projeto` também procura a pasta nesse ambiente.
3. Para confirmar que o MATLAB encontrou a função de configuração, rode `which configurar_projeto` no Command Window. Ela deve apontar para `Analise_BCS_MATLAB/configurar_projeto.m`, e **não** para um diretório externo.
4. Execute tudo, na ordem correta:

   ```matlab
   executar_tudo
   ```

   Ou execute as etapas individualmente:

   ```matlab
   etapa_1_sistema
   etapa_2_linearizacao
   etapa_3_estabilidade
   etapa_4_formas_canonicas
   etapa_5_discretizacao
   ```

   **Atenção:** as etapas 3, 4 e 5 dependem de `outputs/modelo_linear.mat` gerado na etapa 2. A etapa 1 é independente. As etapas 1 e 2 podem consumir mais tempo, por simularem diretamente uma DAE rígida via IDAS.

## Arquivos gerados

- `figures/nonlinear/`: cinco figuras de entradas, cinco de estados e cinco de algébricas do **modelo não linear**.
- `figures/linear/`: gráficos equivalentes do **modelo linearizado**.
- `figures/comparacao/`: degraus pequenos/grandes em `f_BCS_1` e comparação de `q_mean_1` linear x não linear.
- `figures/estabilidade/mapa_de_fase.png`: três projeções 3D.
- `figures/formas_canonicas/`: figura 2×3 das realizações SISO.
- `figures/discretizacao/`: figura 3×3 Euler, ZOH e contínuo.
- `outputs/modelo_linear.mat`: estrutura `lin` com todas as matrizes do modelo, jacobianas, escalas e pontos nominais.
- `outputs/A_B_C_D.mat`, `A.csv`, `B.csv`, `C.csv`: matrizes principais em desvio dimensional.
- `outputs/funcoes_transferencia.mat`, `tf.tex`, `formas_canonicas.mat`, `discretizacao.mat`: análises derivadas.

## Estrutura do pacote

```text
Analise_BCS_MATLAB/
├── configurar_projeto.m  # prepara os caminhos antes de acessar lib/ e modelo/
├── executar_tudo.m
├── etapa_1_sistema.m
├── etapa_2_linearizacao.m
├── etapa_3_estabilidade.m
├── etapa_4_formas_canonicas.m
├── etapa_5_discretizacao.m
├── lib/                    # auxiliares reutilizados por todos os scripts
├── modelo/                 # classes e sym_dynamics do seu sistema
├── modelo_original/        # cópias INTACTAS dos arquivos MATLAB enviados
├── referencias_ipynb/      # cinco notebooks de referência, sem mudanças
├── figures/                # PNGs criados durante execução no MATLAB
└── outputs/                # matrizes/FTs criadas durante execução no MATLAB
```

A única diferença no `modelo/sym_dynamics.m` em relação à cópia intacta em `modelo_original/sym_dynamics.m` é a **remoção do `addpath` com caminho absoluto**, para que o pacote funcione com a instalação local de CasADi. Os parâmetros e as equações físicas foram mantidos.

## Verificação científica e limitações

- O modelo estacionário é o definido no seu `sym_dynamics.m`, incluindo as constantes, equações e parâmetros de ambos os poços. Um cálculo independente das equações, feito para conferir ordem de grandeza, indicou polos aproximadamente em `[-203.1, -149.0, -2.43, -0.408, -0.385, -0.0182, -0.0676, -0.0553]` **em 1/s** no ponto fornecido em `validacao_python.mat`. Os números calculados *pelo MATLAB* serão impressos ao executar a etapa 3. Esses valores de conferência não foram usados como substituto dos cálculos dos scripts.
- A função `minreal` pode cancelar polos/zeros quase coincidentes na realização SISO; o script informa a ordem antes e depois e mantém a matriz completa `G_tf` sem essa redução.
- Este pacote foi criado e revisado em ambiente sem MATLAB executável. **Não afirmamos que a execução MATLAB/CasADi foi testada de ponta a ponta**. Os PNGs e arquivos `.mat` de `outputs` serão efetivamente gerados ao executar no seu MATLAB.

## Correção da versão de 27/10/2026

- `configurar_projeto.m` foi movido de `lib/` para a **raiz** do projeto: o MATLAB consegue encontrá-lo antes de adicionar `lib/` ao path.
- `executar_tudo.m` adiciona explicitamente a raiz ao path.
- Para execução individual, abra a raiz como **Current Folder**.
- Caso utilize a versão antiga do ZIP, a solução temporária, a partir da raiz, é `addpath(fullfile(pwd,'lib')); executar_tudo`.

## Correcao adicional: objeto grafico invalido no ultimo pulso

Se aparecer `Invalid graphics object` dentro de `saveas`/`salvar_figura`,
utilize a versao atualizada de `lib/salvar_figura.m`. Agora ela:

1. Verifica o handle da figura e mostra a causa original de `exportgraphics`.
2. Usa `saveas` somente como alternativa de compatibilidade.
3. **Fecha a janela APOS salvar o PNG**, para que as 15 figuras da etapa 1,
   e as 16 da etapa 2, nao fiquem todas abertas simultaneamente.

Se os avisos indicarem ausencia de `lib/` ou `modelo/`, extraia o ZIP inteiro,
preservando a estrutura. `configurar_projeto.m` agora reconhece tambem
subpastas acidentalmente aninhadas e arquivos `.m` colocados na raiz.

Para corrigir um projeto ja extraido, pode-se copiar apenas
`salvar_figura.m` atualizado para a pasta que o comando `which salvar_figura -all`
indicar. Em seguida, execute `close all; etapa_1_sistema`.


## Correção: áreas dos gráficos pretas no tema escuro

Se os PNGs tiverem o fundo **preto dentro dos subplots** (mesmo com a borda
branca), substitua `configurar_projeto.m` e `lib/salvar_figura.m` pela versão
atualizada. Ambos forçam fundo branco e texto preto em todos os eixos e
legendas, mantendo valores simulados e cores das curvas originais.

O MATLAB precisa **gerar os PNGs novamente**. Para isso, abra a pasta do
projeto e rode:

```matlab
close all
clear functions
rehash
etapa_1_sistema
```

Para as demais figuras, execute as respectivas etapas. Esta correção não
muda parâmetros nem equações físicas do campo BCS.


## Atualizacao v3 — falhas intermitentes de exportacao PNG e CasADi

Esta versao diferencia **resultado numerico** e **exportacao grafica**:

- Antes de tentar gerar as 3 figuras de cada pulso, `etapa_1_sistema.m`
  grava `outputs/simulacoes_nonlineares/pulse-<entrada>.mat`, contendo
  exatamente `t`, `U`, `X` e `Z`.
- `lib/salvar_figura.m` tenta `exportgraphics` a 110 DPI, depois `print`
  com renderizador `painters`, depois `saveas`. Se todos falharem,
  tenta `.fig` e emite **warning**, sem encerrar a simulacao.
  A exportacao de PNG pode continuar falhando conforme o ambiente MATLAB.
- Ao final da etapa, se houver algum PNG ausente, execute
  `regenerar_figuras_etapa_1` para refazer somente os graficos a partir
  dos `.mat`, sem rodar novamente a DAE.
- `configurar_projeto.m` tenta automaticamente o caminho CasADi que ja
  funcionou no seu Windows, sob `%USERPROFILE%/Documents/MATLAB`,
  procurando tambem outras versoes dessa distribuicao.
- Se o MATLAB imprimir um aviso de `Last Working Folder`, basta abrir
  a pasta principal como **Current Folder**. Essa preferencia nao altera
  o modelo nem e causa de falha na exportacao.

### Depois de substituir os arquivos da correcao rapida

```matlab
cd('C:\Users\maxpi\Desktop\Analise_BCS_MATLAB')
clear functions
rehash
which salvar_figura -all
which configurar_projeto -all
etapa_1_sistema
% Se alguma figura faltar: regenerar_figuras_etapa_1
```

O pacote completo preserva as cinco etapas e o modelo fisico. A versao
v3 nao foi executada em MATLAB neste ambiente; a falha grafica e
intermitente e pode depender do renderizador/driver grafico.
