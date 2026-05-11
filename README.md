# 🦠 Calculadora Epidemiológica - Modelo SEIR Estendido

Uma simulação robusta baseada em **modelos compartimentais SEIR estendidos** para prever a dinâmica de propagação de doenças infecciosas em populações. Versão corrigida com validação aprimorada de parâmetros epidemiológicos.

![Python](https://img.shields.io/badge/Python-3.8+-blue)
![License](https://img.shields.io/badge/License-MIT-green)
![Status](https://img.shields.io/badge/Status-Maintenance-yellow)

---

## 📋 Visão Geral

Este projeto implementa um modelo epidemiológico compartimental sofisticado que divide a população em 8 compartimentos distintos:

- **S** → Suscetíveis
- **E** → Expostos (período de incubação)
- **I** → Infectados (período infeccioso)
- **R_mild** → Recuperados por infecção leve
- **R_hosp** → Que necessitam hospitalização
- **H** → Hospitalizados
- **Rec** → Recuperados finais
- **M** → Óbitos

### ✨ Características Principais

✅ **Modelo SEIR Estendido** com compartimentos adicionais para hospitalização e mortalidade  
✅ **Correção de Bugs Críticos** (v2.0):
   - Taxa de hospitalização (tau) agora reflete exatamente o parâmetro declarado
   - Mortalidade hospitalar (delta) com cálculo correto
   - Validação robusta de condições iniciais

✅ **Métricas Avançadas**:
   - Número reprodutivo efetivo (R_eff) em tempo real
   - Incidência diária
   - Tempo de duplicação
   - Conservação de população

✅ **Simulações Precisas** usando integração numérica (scipy.integrate.solve_ivp)

---

## 🚀 Instalação

### Pré-requisitos
- Python 3.8 ou superior
- pip ou conda

### Passos

```bash
# Clone o repositório
git clone https://github.com/seu-usuario/calculadora-epidemiologica.git
cd calculadora-epidemiologica

# Crie um ambiente virtual (recomendado)
python -m venv venv
source venv/bin/activate  # No Windows: venv\Scripts\activate

# Instale as dependências
pip install -r requirements.txt
```

---

## 📦 Dependências

```
numpy>=1.21.0
scipy>=1.7.0
```

---

## 💻 Como Usar

### Exemplo Básico

```python
from Calculadora import EpidemiologicalModel, EpidemiologicalParameters

# Crie os parâmetros padrão
params = EpidemiologicalParameters(
    N=1_000_000,      # População
    R_t=1.5,          # Número reprodutivo efetivo
    T_inf=10.0,       # Dias infeccioso
    T_inc=5.0,        # Dias incubação
    tau=0.15,         # 15% necessitam hospitalização
    delta=0.02,       # 2% mortalidade hospitalar
    I_0=100,          # Infectados iniciais
    E_0=200           # Expostos iniciais
)

# Instancie o modelo
model = EpidemiologicalModel(params)

# Simule 100 dias
resultado = model.simulate(days=100)

# Acesse os resultados
print(f"Pico de infectados: {max(resultado['I']):,.0f} pessoas")
print(f"Total de óbitos: {resultado['M'][-1]:,.0f}")
```

### Parâmetros Principais

| Parâmetro | Descrição | Padrão |
|-----------|-----------|--------|
| `N` | População total | 1.000.000 |
| `R_t` | Número reprodutivo efetivo | 1.5 |
| `T_inf` | Duração período infeccioso (dias) | 10.0 |
| `T_inc` | Duração período incubação (dias) | 5.0 |
| `tau` | Fração que requer hospitalização | 0.15 |
| `delta` | Fração hospitalizados que morrem | 0.02 |
| `T_rec` | Tempo recuperação casos leves (dias) | 14.0 |
| `T_hosp` | Tempo até internamento (dias) | 3.0 |
| `T_int` | Duração internamento (dias) | 10.0 |
| `I_0` | Infectados iniciais | 100 |
| `E_0` | Expostos iniciais | 200 |

---

## 📊 Estrutura de Saída

O método `simulate()` retorna um dicionário com as séries temporais de cada compartimento:

```python
resultado = {
    't': array([0, 1, 2, ...]),      # Vetor de tempo
    'S': array([...]),                # Suscetíveis
    'E': array([...]),                # Expostos
    'I': array([...]),                # Infectados
    'R_mild': array([...]),           # Recuperados leves
    'R_hosp': array([...]),           # Pré-hospitalizados
    'H': array([...]),                # Hospitalizados
    'Rec': array([...]),              # Recuperados finais
    'M': array([...])                 # Óbitos
}
```

---

## 🔬 Fundamentos Matemáticos

### Equações Diferenciais

$$\frac{dS}{dt} = -\beta \cdot \frac{S \cdot I}{N}$$

$$\frac{dE}{dt} = \beta \cdot \frac{S \cdot I}{N} - \sigma \cdot E$$

$$\frac{dI}{dt} = \sigma \cdot E - \gamma \cdot I$$

$$\frac{dR_{mild}}{dt} = (1-\tau) \cdot \gamma \cdot I - \frac{1}{T_{rec}} \cdot R_{mild}$$

$$\frac{dR_{hosp}}{dt} = \tau \cdot \gamma \cdot I - \frac{1}{T_{hosp}} \cdot R_{hosp}$$

$$\frac{dH}{dt} = \frac{1}{T_{hosp}} \cdot R_{hosp} - \frac{1}{T_{int}} \cdot H$$

$$\frac{dRec}{dt} = \frac{1}{T_{rec}} \cdot R_{mild} + \frac{1-\delta}{T_{int}} \cdot H$$

$$\frac{dM}{dt} = \frac{\delta}{T_{int}} \cdot H$$

### Parâmetros Derivados

- **β** (taxa transmissão): $\beta = \frac{R_t}{T_{inf}}$
- **σ** (taxa incubação): $\sigma = \frac{1}{T_{inc}}$
- **γ** (taxa infecciosidade): $\gamma = \frac{1}{T_{inf}}$

---

## 🧪 Testes

Execute os testes incluídos:

```bash
# Teste rápido
python teste.py

# Teste detalhado com verificação de conservação
python teste_detalhado.py
```

O teste detalhado verifica:
- Conservação de população total
- Valores finais de cada compartimento
- Pico de infectados
- Integridade dos dados

---

## ⚠️ Mudanças na v2.0

### Bugs Corrigidos

| Bug | Impacto | Correção |
|-----|---------|----------|
| Taxa hospitalização (tau) | Fração real divergia do parâmetro | Bifurcação explícita em I |
| Mortalidade hospitalar (delta) | Dependia de múltiplos parâmetros | Taxa única em H |
| Validação de S_0 | Aceitava condições negativas | Validação em __post_init__ |

### Melhorias

- Compartimento R dividido em R_mild e R_hosp para clareza
- Adicionadas métricas: R_eff(t), incidência diária, tempo de duplicação
- Função _enforce_constraints reescrita para maior conservatividade

---

## 📚 Referências

- Kermack-McKendrick (1927) - Modelo SIR clássico
- Compartmental models in epidemiology
- Bifurcação explícita para modelagem de múltiplos desfechos

---

## 🤝 Contribuindo

Contribuições são bem-vindas! Para reportar bugs ou sugerir melhorias:

1. Abra uma [Issue](../../issues)
2. Descreva o problema com detalhes
3. Inclua exemplos se possível

### Padrão de Código
- Siga PEP 8
- Adicione docstrings em funções públicas
- Inclua testes para novas funcionalidades

---

## 📄 Licença

Este projeto está licenciado sob a Licença MIT - veja o arquivo [LICENSE](LICENSE) para detalhes.

---

## 👤 Autor

Desenvolvido como parte de estudos em modelagem epidemiológica.

---

## ⚡ Dicas de Uso

### Simulando diferentes cenários

```python
# Cenário 1: Sem intervenção
params_sem = EpidemiologicalParameters(R_t=2.5)
modelo_sem = EpidemiologicalModel(params_sem)
resultado_sem = modelo_sem.simulate(100)

# Cenário 2: Com isolamento (R_t reduzido)
params_com = EpidemiologicalParameters(R_t=0.8)
modelo_com = EpidemiologicalModel(params_com)
resultado_com = modelo_com.simulate(100)

# Comparar picos
print(f"Sem isolamento: {max(resultado_sem['I']):,.0f}")
print(f"Com isolamento: {max(resultado_com['I']):,.0f}")
```

### Validação de Parâmetros

```python
try:
    params = EpidemiologicalParameters(tau=1.5)  # ❌ Erro!
except ValueError as e:
    print(f"Parâmetro inválido: {e}")
```

---

**Desenvolvido com ❤️ para pesquisa em epidemiologia**
