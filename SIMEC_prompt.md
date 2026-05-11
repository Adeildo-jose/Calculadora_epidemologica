# Prompt para Geração da Calculadora Epidemiológica — SIMEC

Atue como um Engenheiro de Software Científico sênior e Especialista em Epidemiologia Computacional. Seu objetivo é desenvolver um código em Python (utilizando programação orientada a objetos) para uma Calculadora Epidemiológica baseada em um modelo compartimental SEIR estendido. O código deve ser modular, escalável, tipado (Type Hinting) e possuir docstrings completas.

A dinâmica clínica deste modelo subdivide a saída do compartimento infeccioso para simular a progressão da doença em alta resolução (casos leves, hospitalizações e óbitos).

---

## 1. O Modelo Matemático (Sistema de EDOs)

As transições de estado (taxas de variação diária) são definidas pelo seguinte conjunto de equações diferenciais ordinárias:

### Dinâmica de Transmissão (SEIR base)

$$\frac{dS}{dt} = -R_t \cdot \frac{1}{T_{inf}} \cdot \frac{I \cdot S}{N}$$

$$\frac{dE}{dt} = R_t \cdot \frac{1}{T_{inf}} \cdot \frac{I \cdot S}{N} - \frac{1}{T_{inc}} \cdot E$$

$$\frac{dI}{dt} = \frac{1}{T_{inc}} \cdot E - \frac{1}{T_{inf}} \cdot I$$

### Dinâmica Clínica e Desfechos (Compartimentos Auxiliares)

> **Nota:** O compartimento $R$ aqui age como um estado de transição pós-infeccioso, derivando para recuperação leve ou hospitalização.

$$\frac{dR}{dt} = \frac{1}{T_{inf}} \cdot I - \frac{1}{T_{rec}} \cdot R - \frac{\tau}{T_{hosp}} \cdot R$$

$$\frac{dH}{dt} = \frac{\tau}{T_{hosp}} \cdot R - \frac{\delta}{T_m} \cdot H - \frac{1}{T_{int}} \cdot H$$

$$\frac{dRec}{dt} = \frac{1}{T_{rec}} \cdot R + \frac{1}{T_{int}} \cdot H$$

$$\frac{dM}{dt} = \frac{\delta}{T_m} \cdot H$$

---

## 2. Parâmetros do Modelo

O código deve aceitar uma estrutura de dados (ex: `dataclass` ou dicionário estruturado) para os seguintes parâmetros, permitindo fácil calibração:

| Parâmetro | Descrição |
|-----------|-----------|
| $N$ | População Total: população inicial suscetível (assumindo $S \approx N$ em $t=0$) |
| $R_t$ | Número reprodutivo efetivo, onde $R_t = (1 - \theta) \cdot R_0$, e $\theta$ é o fator de redução na transmissão (isolamento/vacinação) |
| $T_{inf}$ | Tempo em dias em que o paciente é infeccioso |
| $T_{inc}$ | Tempo de incubação em dias |
| $\delta$ | Taxa de mortalidade: fração dos hospitalizados que evoluem para óbito |
| $\tau$ | Taxa de internamento: fração dos sintomáticos que requerem hospitalização |
| $T_m$ | Tempo decorrido até a morte (para os casos fatais) |
| $T_{int}$ | Tempo de internamento hospitalar médio (para os que sobrevivem ou até o desfecho) |
| $T_{rec}$ | Tempo de recuperação para casos leves |
| $T_{hosp}$ | Tempo decorrido do fim do período infeccioso até a efetiva hospitalização |

---

## 3. Requisitos Técnicos Obrigatórios

### Motor de Integração
Utilize o `scipy.integrate.solve_ivp` para resolver o sistema de EDOs. Garanta precisão na integração definindo métodos robustos (ex: `'LSODA'` ou `'RK45'`).

### Arquitetura
Crie uma classe `EpidemiologicalModel` que inicialize os parâmetros e as condições iniciais ($S_0, E_0, I_0, R_0, H_0, Rec_0, M_0$).

### Métodos da Classe

- `_derivadas(t, y)` — Método interno contendo as equações diferenciais.
- `simulate(days)` — Executa a simulação pelo número de dias estipulado.
- `plot_results()` — Gera um gráfico profissional utilizando a biblioteca `Plotly` (para interatividade) ou `Matplotlib`. O gráfico deve ter eixos nomeados, legenda clara e distinguir bem os compartimentos clínicos (especialmente $H$ — Hospitalizados e $M$ — Óbitos).

### Tratamento de Erros
Adicione verificações para garantir que nenhuma variável populacional fique negativa devido a erros de precisão de ponto flutuante, e garanta que $S + E + I + R + H + Rec + M = N$ em qualquer instante $t$.

### Exemplo de Uso
Ao final do script, inclua um bloco `if __name__ == "__main__":` instanciando o modelo com valores hipotéticos realistas de uma epidemia respiratória para testar a implementação.
