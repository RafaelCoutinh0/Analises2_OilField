# Análise dinâmica e representação em espaço de estados de um campo de produção de petróleo operado por BCS

**Rafael Costa Coutinho e Guilherme Rocha**  
Universidade Federal da Bahia (UFBA) · Engenharia de Controle e Automação Industrial

[**Ler o relatório completo (PDF)**](Latex/relatorio.pdf) · [**Consultar o arquivo LaTeX**](Latex/main.tex)

## Resumo

Este trabalho apresenta a análise dinâmica de um campo de produção de petróleo composto por **dois poços equipados com Bombeio Centrífugo Submerso (BCS)**, conectados a um *manifold* comum e a uma bomba *booster*. Com base em uma adaptação do modelo fenomenológico de Rønning (2011), foram investigadas as respostas a perturbações operacionais, a linearização em torno de um ponto estacionário, a estabilidade local, as formas canônicas e a discretização da representação linear.

O objetivo é caracterizar a dinâmica do campo e verificar a capacidade do modelo linearizado de reproduzir o comportamento do sistema nas proximidades do ponto de operação.

**Palavras-chave:** Bombeio Centrífugo Submerso; modelo fenomenológico; espaço de estados; estabilidade; discretização.

## 1. Campo de produção

O sistema possui dois poços, cada um equipado com uma BCS e uma válvula *choke*. A produção dos poços converge para um *manifold*, enquanto a bomba *booster* auxilia no transporte do fluido até as instalações de superfície. Por compartilharem a pressão do *manifold*, os poços apresentam acoplamento hidráulico: uma alteração operacional em um poço pode afetar o outro.

![Fluxograma do campo de produção de petróleo](Latex/figuras/20226.2model.drawio.png)

A representação fenomenológica é dada por um sistema de equações diferenciais-algébricas (DAE):

$$
\dot{\mathbf{x}}=f(\mathbf{x},\mathbf{z},\mathbf{u}),
\qquad
\mathbf{0}=g(\mathbf{x},\mathbf{z},\mathbf{u}).
$$

| Conjunto | Dimensão | Variáveis |
|---|---:|---|
| Estados diferenciais $\mathbf{x}$ | 8 | Pressão do *manifold*, vazão de transporte e, para cada poço, pressão de fundo, pressão na *choke* e vazão média |
| Variáveis algébricas $\mathbf{z}$ | 4 | Pressão de entrada e incremento de pressão de cada BCS |
| Entradas $\mathbf{u}$ | 6 | Frequência da *booster*, pressão *topside*, frequências das duas BCS e aberturas das duas *chokes* |

Nas análises de perturbação, a pressão *topside* foi mantida constante; assim, foram consideradas **cinco entradas manipuladas**.

### Ponto estacionário

Para a condição de operação avaliada, foram adotadas frequência da *booster* de **45 Hz**, frequências das BCS de **50 Hz**, aberturas das *chokes* iguais a **0,8** e pressão *topside* de **1 bar**. Entre os valores estacionários obtidos, destacam-se:

| Variável / unidade | Valor estacionário |
|---|---:|
| $p_{\mathrm{man}}/\mathrm{bar}$ | 66,9278 |
| $q_{\mathrm{tr}}/(\mathrm{m}^3/\mathrm{h})$ | 175,9557 |
| $\bar q_1/(\mathrm{m}^3/\mathrm{h})$ | 90,4961 |
| $\bar q_2/(\mathrm{m}^3/\mathrm{h})$ | 85,4595 |
| $\Delta p_{\mathrm{BCS},1}/\mathrm{bar}$ | 76,7618 |
| $\Delta p_{\mathrm{BCS},2}/\mathrm{bar}$ | 82,5215 |

Foram aplicadas perturbações individuais nas entradas manipuladas para observar os transientes, o acoplamento hidráulico e o retorno ao regime estacionário.

![Respostas dos estados a uma perturbação na frequência da BCS do poço 1](Latex/figuras/sistema_saidas.png)

## 2. Linearização em espaço de estados

O sistema DAE foi linearizado em torno do ponto estacionário, por expansão de Taylor de primeira ordem. As variáveis foram expressas em **desvios dimensionais**, definidos por $\hat{\mathbf{x}}=\mathbf{x}-\mathbf{x}_{ss}$ e $\hat{\mathbf{u}}=\mathbf{u}-\mathbf{u}_{ss}$.

Após a eliminação das variáveis algébricas, obteve-se:

$$
\dot{\hat{\mathbf{x}}}=\mathbf{A}\hat{\mathbf{x}}+\mathbf{B}\hat{\mathbf{u}},
\qquad
\hat{\mathbf{y}}=\mathbf{C}\hat{\mathbf{x}},
$$

com

$$
\mathbf{A}=\mathbf{F}_x-\mathbf{F}_z\mathbf{G}_z^{-1}\mathbf{G}_x,
\qquad
\mathbf{B}=\mathbf{F}_u-\mathbf{F}_z\mathbf{G}_z^{-1}\mathbf{G}_u.
$$

Na representação estudada, $\mathbf{A}\in\mathbb{R}^{8\times8}$, $\mathbf{B}\in\mathbb{R}^{8\times5}$ e $\mathbf{C}=\mathbf{I}_8$. As variáveis algébricas foram recuperadas separadamente pela relação linearizada.

A representatividade local foi investigada comparando a resposta do incremento de pressão da BCS do poço 1, $\Delta p_{\mathrm{BCS},1}$, nos modelos linear e não linear, diante de alterações na frequência de operação de **50 para 52 Hz** e, posteriormente, para **65 Hz**.

![Comparação entre as respostas linear e não linear para o incremento de pressão da BCS 1](Latex/figuras/comparacao.png)

A aproximação linear apresentou maior concordância nas proximidades do ponto de operação; diferenças mais evidentes surgiram para perturbações maiores.

## 3. Análise de estabilidade

A estabilidade local foi investigada pelos autovalores da matriz dinâmica $\mathbf{A}$. Todos os oito autovalores apresentaram parte real negativa, satisfazendo a condição de estabilidade assintótica do modelo linear contínuo:

$$
\operatorname{Re}(\lambda_i)<0.
$$

Além disso, foram simuladas **12 condições iniciais** próximas ao equilíbrio, com variações de até $\pm5\%$ dos valores estacionários. As trajetórias do sistema linear homogêneo convergiram para a origem, em concordância com os autovalores.

![Projeções tridimensionais do espaço de estados](Latex/figuras/estabilidade.png)

A conclusão refere-se à **estabilidade local**, não à estabilidade global do sistema não linear.

## 4. Formas canônicas

A partir do modelo linear, foi obtida a matriz de funções de transferência:

$$
\mathbf{G}(s)=\mathbf{C}(s\mathbf{I}-\mathbf{A})^{-1}\mathbf{B}.
$$

A matriz possui **8 saídas × 5 entradas = 40 funções de transferência SISO**. Para a análise das formas canônicas, foi escolhida a relação entre a frequência da BCS do poço 1 e sua vazão média. Após redução por cancelamentos numéricos de polos e zeros, a realização selecionada apresentou **sexta ordem**.

Foram analisadas três representações:

- **Canônica controlável**, com coeficientes do polinômio característico organizados na matriz dinâmica.
- **Canônica observável**, obtida por dualidade da realização controlável.
- **Forma de Jordan**, obtida por transformação de similaridade, com matriz dinâmica aproximadamente diagonal.

![Comparação das respostas das três formas canônicas](Latex/figuras/canonicas.png)

As realizações apresentam diferentes coordenadas internas, mas preservam a mesma relação entrada–saída, dentro da precisão numérica.

## 5. Discretização

A representação contínua foi discretizada pelos métodos de **Euler explícito** e de **retenção de ordem zero (ZOH)**.

| Método | Matriz de estados | Matriz de entradas |
|---|---|---|
| Euler | $\mathbf{A}_{d,E}=\mathbf{I}+T_s\mathbf{A}$ | $\mathbf{B}_{d,E}=T_s\mathbf{B}$ |
| ZOH | $\mathbf{A}_{d,Z}=e^{\mathbf{A}T_s}$ | $\mathbf{B}_{d,Z}=\int_0^{T_s}e^{\mathbf{A}\tau}\mathbf{B}\,d\tau$ |

O período de amostragem foi definido a partir do limite de estabilidade de Euler, adotando-se $T_s\approx0{,}007877\,\mathrm{s}$, equivalente a 80% do limite calculado. A análise comparou o sistema contínuo e os discretizados para um degrau de **+1 Hz** na frequência da BCS do primeiro poço.

![Comparação dos modelos contínuo, Euler e ZOH](Latex/figuras/discretizacao.png)

As respostas discretas permaneceram próximas às do modelo contínuo para o período de amostragem considerado.

## Estrutura da pasta LaTeX

```text
Latex/
├── main.tex                     # Fonte do relatório
├── relatorio.pdf                # Relatório compilado
├── figuras/
│   ├── 20226.2model.drawio.png # Fluxograma do campo
│   ├── sistema_entradas.png     # Entradas perturbadas
│   ├── sistema_saidas.png       # Estados diferenciais
│   ├── sistema_algebricas.png   # Variáveis algébricas
│   ├── comparacao.png           # Modelo linear × não linear
│   ├── estabilidade.png         # Projeções de estabilidade
│   ├── canonicas.png            # Formas canônicas
│   └── discretizacao.png        # Discretização
└── LEIA_ME.md                   # Notas auxiliares do pacote
```

## Como compilar o relatório

É necessário dispor de uma distribuição LaTeX com suporte a **pdfLaTeX** e aos pacotes declarados em `main.tex`.

**Overleaf:** crie um projeto, importe o conteúdo da pasta `Latex/` e compile `main.tex` utilizando pdfLaTeX.

**Localmente:** a partir da raiz do repositório, execute:

```bash
cd Latex
pdflatex main.tex
pdflatex main.tex
```

O resultado será `main.pdf`. O arquivo `relatorio.pdf` disponibilizado no repositório corresponde à versão compilada incluída no pacote.

> **Escopo:** a pasta `Latex/` contém o relatório e suas figuras. As simulações e o processamento numérico que originaram os resultados foram realizados em MATLAB; os arquivos necessários à compilação do relatório não constituem, por si só, uma implementação executável do modelo.

## Considerações finais

As análises evidenciaram o acoplamento entre os poços, a estabilidade assintótica local da representação linearizada, a equivalência entrada–saída das realizações canônicas e a proximidade entre as respostas contínuas e discretas para o período de amostragem considerado. O trabalho fornece uma base para estudos posteriores de análise e controle de campos de produção operados por BCS.

## Referência

RØNNING, R. *Automatic Start-up and Control of Artificially Lifted Wells*. 2011.
