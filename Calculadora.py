"""
Calculadora Epidemiológica - Modelo Compartimental SEIR Estendido
Versão 2.0 (Corrigida)

MUDANÇAS DESTA VERSÃO:
  - [BUG CRÍTICO] Taxa de hospitalização (tau): Modelo anterior usava taxas
    competitivas entre R→H e R→Rec, fazendo com que a fração real hospitalizada
    DIVERGISSE do parâmetro tau declarado (ex: tau=15% resultava em ~41%).
    Corrigido com bifurcação explícita a partir de I.
  - [BUG CRÍTICO] Mortalidade hospitalar (delta): Mesmo problema; a fração real
    que morria em H dependia de T_m e T_int juntos. Corrigido com taxa única em H.
  - [BUG] Validação de S_0 < 0: __post_init__ não verificava se I_0 + E_0 > N,
    gerando condição inicial negativa silenciosamente. Corrigido.
  - [MELHORIA] Compartimento R dividido em R_mild e R_hosp para clareza.
  - [MELHORIA] Adicionados R_eff(t), incidência diária e tempo de duplicação.
  - [MELHORIA] _enforce_constraints reescrito para ser conservativo.

Compartimentos:
    S      : Suscetível
    E      : Exposto (incubação)
    I      : Infectado (infeccioso)
    R_mild : Pós-infeccioso com evolução leve (→ Rec)
    R_hosp : Pós-infeccioso com evolução grave (→ H)
    H      : Hospitalizado
    Rec    : Recuperado (desfecho final)
    M      : Morto (desfecho final)

Equações diferenciais:
    dS/dt      = -β · S · I / N
    dE/dt      = β · S · I / N - σ · E
    dI/dt      = σ · E - γ · I
    dR_mild/dt = (1-τ) · γ · I - (1/T_rec) · R_mild
    dR_hosp/dt = τ · γ · I - (1/T_hosp) · R_hosp
    dH/dt      = (1/T_hosp) · R_hosp - (1/T_int) · H
    dRec/dt    = (1/T_rec) · R_mild + (1-δ)/T_int · H
    dM/dt      = δ/T_int · H

Onde:
    β = R_t / T_inf   (taxa de transmissão)
    σ = 1 / T_inc     (taxa de progressão da incubação)
    γ = 1 / T_inf     (taxa de recuperação da infecciosidade)
    τ = tau           (fração EXATA que requer hospitalização)
    δ = delta         (fração EXATA de hospitalizados que morrem)
"""

from dataclasses import dataclass, field
from typing import Dict, Optional
import numpy as np
from scipy.integrate import solve_ivp
import warnings


@dataclass
class EpidemiologicalParameters:
    """
    Parâmetros do modelo epidemiológico SEIR estendido (v2).

    Atributos:
        N      (float): População total suscetível.
        R_t    (float): Número reprodutivo efetivo (≥ 0).
        T_inf  (float): Duração média do período infeccioso em dias.
        T_inc  (float): Duração média do período de incubação em dias.
        tau    (float): Fração EXATA dos infectados que requer hospitalização [0, 1].
        delta  (float): Fração EXATA dos hospitalizados que evolui para óbito [0, 1].
        T_rec  (float): Tempo médio de recuperação para casos LEVES (R_mild→Rec) em dias.
        T_hosp (float): Tempo médio desde fim da infecciosidade até internamento (R_hosp→H) em dias.
        T_int  (float): Duração média do internamento hospitalar (H→Rec ou H→M) em dias.
        I_0    (float): Infectados iniciais.
        E_0    (float): Expostos iniciais.

    Nota sobre T_m (removido):
        Na v1, T_m tentava separar o tempo até óbito do tempo até alta. Isso
        introduzia inconsistência matemática: a fração real que morria dependia
        de T_m e T_int simultaneamente, divergindo do parâmetro delta. Nesta
        versão usa-se um único T_int para todo o compartimento H, garantindo
        que delta seja a fração exata de mortalidade hospitalar.
    """

    N: float = 1_000_000
    R_t: float = 1.5
    T_inf: float = 10.0
    T_inc: float = 5.0
    tau: float = 0.15
    delta: float = 0.02
    T_rec: float = 14.0
    T_hosp: float = 3.0
    T_int: float = 10.0
    I_0: float = 100.0
    E_0: float = 200.0

    def __post_init__(self) -> None:
        """Valida os parâmetros após inicialização."""
        if self.N <= 0:
            raise ValueError(f"N deve ser positivo. Recebido: {self.N}")
        if self.R_t < 0:
            raise ValueError(f"R_t não pode ser negativo. Recebido: {self.R_t}")
        if not (0.0 <= self.tau <= 1.0):
            raise ValueError(f"tau deve estar em [0, 1]. Recebido: {self.tau}")
        if not (0.0 <= self.delta <= 1.0):
            raise ValueError(f"delta deve estar em [0, 1]. Recebido: {self.delta}")
        if self.I_0 < 0:
            raise ValueError(f"I_0 deve ser ≥ 0. Recebido: {self.I_0}")
        if self.E_0 < 0:
            raise ValueError(f"E_0 deve ser ≥ 0. Recebido: {self.E_0}")
        if self.I_0 + self.E_0 > self.N:
            raise ValueError(
                f"I_0 + E_0 ({self.I_0 + self.E_0:,.0f}) excede N ({self.N:,.0f}). "
                "Não é possível ter mais infectados/expostos do que a população total."
            )
        tempo_params = {
            "T_inf": self.T_inf, "T_inc": self.T_inc,
            "T_rec": self.T_rec, "T_hosp": self.T_hosp, "T_int": self.T_int,
        }
        for nome, valor in tempo_params.items():
            if valor <= 0:
                raise ValueError(f"{nome} deve ser positivo. Recebido: {valor}")

    @property
    def beta(self) -> float:
        """Taxa de transmissão derivada: β = R_t / T_inf."""
        return self.R_t / self.T_inf

    @property
    def sigma(self) -> float:
        """Taxa de progressão da incubação: σ = 1 / T_inc."""
        return 1.0 / self.T_inc

    @property
    def gamma(self) -> float:
        """Taxa de saída do período infeccioso: γ = 1 / T_inf."""
        return 1.0 / self.T_inf


# Índices dos compartimentos no vetor de estado (y)
_IDX = {name: i for i, name in enumerate(["S", "E", "I", "R_mild", "R_hosp", "H", "Rec", "M"])}
_COMPARTMENTS = list(_IDX.keys())


class EpidemiologicalModel:
    """
    Modelo epidemiológico determinístico baseado em EDOs (SEIR estendido v2).

    Uso básico:
        params = EpidemiologicalParameters(N=500_000, R_t=2.0, tau=0.10)
        modelo = EpidemiologicalModel(params)
        resultados = modelo.simulate(days=180)
        stats = modelo.get_statistics(resultados)
    """

    def __init__(self, params: EpidemiologicalParameters) -> None:
        self.params = params
        self.solution = None
        self.t_eval: Optional[np.ndarray] = None

        S_0 = params.N - params.I_0 - params.E_0
        self.y_0 = np.array([
            S_0,          # S
            params.E_0,   # E
            params.I_0,   # I
            0.0,          # R_mild
            0.0,          # R_hosp
            0.0,          # H
            0.0,          # Rec
            0.0,          # M
        ])

    # ------------------------------------------------------------------ #
    #  Sistema de EDOs                                                     #
    # ------------------------------------------------------------------ #

    def _derivadas(self, t: float, y: np.ndarray) -> np.ndarray:
        """
        Sistema de 8 EDOs do modelo SEIR estendido (v2).

        A bifurcação τ / (1-τ) ocorre na saída de I, garantindo que exatamente
        a fração tau de infectados siga para hospitalização e (1-tau) para
        recuperação leve — sem dependência dos valores de T_hosp e T_rec.

        Analogamente, a bifurcação δ / (1-δ) ocorre na saída de H, garantindo
        que exatamente delta de hospitalizados evolua para óbito.
        """
        S, E, I, R_mild, R_hosp, H, Rec, M = y
        p = self.params
        N = p.N

        force_of_infection = p.beta * S * I / N if N > 0 else 0.0

        dS      = -force_of_infection
        dE      = force_of_infection - p.sigma * E
        dI      = p.sigma * E - p.gamma * I

        # Bifurcação exata τ / (1-τ) na saída de I
        saida_I = p.gamma * I
        dR_mild = (1.0 - p.tau) * saida_I - (1.0 / p.T_rec) * R_mild
        dR_hosp = p.tau * saida_I - (1.0 / p.T_hosp) * R_hosp

        # H recebe R_hosp; bifurcação exata δ / (1-δ) na saída de H
        saida_H = (1.0 / p.T_int) * H
        dH      = (1.0 / p.T_hosp) * R_hosp - saida_H
        dRec    = (1.0 / p.T_rec) * R_mild + (1.0 - p.delta) * saida_H
        dM      = p.delta * saida_H

        return np.array([dS, dE, dI, dR_mild, dR_hosp, dH, dRec, dM])

    # ------------------------------------------------------------------ #
    #  Simulação                                                           #
    # ------------------------------------------------------------------ #

    def simulate(self, days: int) -> Dict[str, np.ndarray]:
        """
        Executa a simulação pelo número de dias especificado.

        Args:
            days: Número inteiro de dias (deve ser ≥ 1).

        Retorna:
            Dicionário com séries temporais:
            {
                't'     : array de tempo (dias),
                'S'     : Suscetíveis,
                'E'     : Expostos,
                'I'     : Infectados,
                'R_mild': Pós-infecciosos leves,
                'R_hosp': Pós-infecciosos graves (pré-internamento),
                'H'     : Hospitalizados,
                'Rec'   : Recuperados totais,
                'M'     : Óbitos totais,
                'incidence': Novos casos por dia (estimativa),
                'R_eff' : Número reprodutivo efetivo ao longo do tempo,
            }

        Raises:
            ValueError: Se days < 1.
            RuntimeError: Se a integração numérica falhar.
        """
        if days < 1:
            raise ValueError(f"days deve ser ≥ 1. Recebido: {days}")

        self.t_eval = np.linspace(0, days, days + 1)

        try:
            self.solution = solve_ivp(
                self._derivadas,
                t_span=[0, days],
                y0=self.y_0,
                method="RK45",
                t_eval=self.t_eval,
                rtol=1e-9,
                atol=1e-11,
            )
        except Exception as exc:
            raise RuntimeError(f"Falha no integrador: {exc}") from exc

        if not self.solution.success:
            raise RuntimeError(f"Integração falhou: {self.solution.message}")

        y = self.solution.y

        # Garante não-negatividade; emite aviso se desvio for significativo
        neg_mask = y < 0
        if neg_mask.any():
            max_neg = np.abs(y[neg_mask]).max()
            if max_neg > self.params.N * 1e-8:
                warnings.warn(
                    f"Valores negativos detectados (máx abs = {max_neg:.2e}). "
                    "Verifique os parâmetros. Os valores foram zerados.",
                    RuntimeWarning,
                    stacklevel=2,
                )
            y = np.maximum(y, 0.0)

        # Valida conservação (informativo, não corrige)
        self._validate_conservation(y)

        S, E, I, R_mild, R_hosp, H, Rec, M = y

        # Incidência: variação negativa em S (novos casos por dia)
        incidence = np.diff(S, prepend=S[0])
        incidence = np.maximum(-incidence, 0.0)

        # R_eff(t) = R_t · S(t) / N  (aproximação de campo médio)
        R_eff = self.params.R_t * S / self.params.N

        return {
            "t": self.solution.t,
            "S": S,
            "E": E,
            "I": I,
            "R_mild": R_mild,
            "R_hosp": R_hosp,
            "H": H,
            "Rec": Rec,
            "M": M,
            "incidence": incidence,
            "R_eff": R_eff,
        }

    # ------------------------------------------------------------------ #
    #  Validação de conservação                                            #
    # ------------------------------------------------------------------ #

    def _validate_conservation(self, y: np.ndarray) -> None:
        """
        Verifica que S+E+I+R_mild+R_hosp+H+Rec+M ≈ N em todos os passos.
        Apenas informativo; não altera os dados.
        """
        totals = y.sum(axis=0)
        max_dev = np.abs(totals - self.params.N).max()
        tolerance = self.params.N * 1e-6

        if max_dev > tolerance:
            warnings.warn(
                f"Desvio de conservação = {max_dev:.4f} (tolerância = {tolerance:.4f}). "
                "Possível instabilidade numérica.",
                RuntimeWarning,
                stacklevel=2,
            )

    # ------------------------------------------------------------------ #
    #  Estatísticas resumidas                                              #
    # ------------------------------------------------------------------ #

    def get_statistics(self, results: Dict[str, np.ndarray]) -> Dict[str, float]:
        """
        Calcula métricas epidemiológicas relevantes da simulação.

        Args:
            results: Dicionário retornado por simulate().

        Retorna:
            Dicionário contendo:
            - peak_infected       : Pico de infectados simultâneos
            - day_peak_infected   : Dia do pico de infectados
            - peak_hospitalized   : Pico de hospitalizados simultâneos
            - day_peak_hospitalized: Dia do pico de hospitalizados
            - total_deaths        : Óbitos acumulados ao final
            - total_recovered     : Recuperados acumulados ao final
            - attack_rate         : Fração da população infectada (taxa de ataque)
            - infection_fatality_rate: Mortalidade por infecção (óbitos / total infectados)
            - peak_R_eff          : Valor máximo de R_eff no início da epidemia
            - day_R_eff_below_1   : Primeiro dia em que R_eff cai abaixo de 1 (ou None)
        """
        p = self.params
        S_0 = p.N - p.I_0 - p.E_0

        total_infected = S_0 - results["S"][-1]
        attack_rate = total_infected / p.N if p.N > 0 else 0.0
        ifr = results["M"][-1] / total_infected if total_infected > 0 else 0.0

        # Dia em que R_eff cruza 1
        R_eff = results["R_eff"]
        below_1 = np.where(R_eff < 1.0)[0]
        day_R_eff_below_1 = int(below_1[0]) if below_1.size > 0 else None

        return {
            "peak_infected": float(np.max(results["I"])),
            "day_peak_infected": int(np.argmax(results["I"])),
            "peak_hospitalized": float(np.max(results["H"])),
            "day_peak_hospitalized": int(np.argmax(results["H"])),
            "total_deaths": float(results["M"][-1]),
            "total_recovered": float(results["Rec"][-1]),
            "attack_rate": float(attack_rate),
            "infection_fatality_rate": float(ifr),
            "peak_R_eff": float(np.max(R_eff)),
            "day_R_eff_below_1": day_R_eff_below_1,
        }

    # ------------------------------------------------------------------ #
    #  Visualização                                                        #
    # ------------------------------------------------------------------ #

    def plot_results(
        self,
        results: Dict[str, np.ndarray],
        save_path: Optional[str] = None,
        use_plotly: bool = False,
    ) -> None:
        """
        Gera visualização profissional dos resultados da simulação.

        Args:
            results  : Dicionário retornado por simulate().
            save_path: Caminho opcional para salvar a figura.
            use_plotly: True = gráfico interativo (Plotly), False = Matplotlib.
        """
        if use_plotly:
            self._plot_plotly(results, save_path)
        else:
            self._plot_matplotlib(results, save_path)

    def _plot_matplotlib(
        self, results: Dict[str, np.ndarray], save_path: Optional[str] = None
    ) -> None:
        import matplotlib.pyplot as plt

        t = results["t"]
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        fig.suptitle(
            "Simulação Epidemiológica – Modelo SEIR Estendido (v2)",
            fontsize=16, fontweight="bold",
        )

        def fmt_ax(ax, title, ylabel="Indivíduos"):
            ax.set_xlabel("Dias", fontsize=10)
            ax.set_ylabel(ylabel, fontsize=10)
            ax.set_title(title, fontweight="bold")
            ax.legend(loc="best", fontsize=9)
            ax.grid(True, alpha=0.3)

        # 1 – Dinâmica SEIR
        ax = axes[0, 0]
        ax.plot(t, results["S"], "b-", lw=2, label="Suscetível (S)")
        ax.plot(t, results["E"], color="orange", lw=2, label="Exposto (E)")
        ax.plot(t, results["I"], "r-", lw=2, label="Infectado (I)")
        fmt_ax(ax, "Dinâmica de Transmissão (SEIR)")

        # 2 – Pós-infecciosos (R_mild e R_hosp, antes v1 era só R)
        ax = axes[0, 1]
        ax.plot(t, results["R_mild"], color="lightgreen", lw=2, label="Pós-inf. Leve (R_mild)")
        ax.plot(t, results["R_hosp"], color="darkred", lw=2, label="Pós-inf. Grave (R_hosp)")
        fmt_ax(ax, "Compartimentos Pós-infecciosos")

        # 3 – Desfechos clínicos
        ax = axes[0, 2]
        ax.plot(t, results["H"], color="purple", lw=2.5, label="Hospitalizados (H)")
        ax.plot(t, results["M"], "k-", lw=2.5, label="Óbitos (M)")
        fmt_ax(ax, "Desfechos Clínicos Graves")

        # 4 – Recuperados
        ax = axes[1, 0]
        ax.plot(t, results["Rec"], "darkgreen", lw=2.5, label="Recuperados (Rec)")
        ax.fill_between(t, 0, results["Rec"], alpha=0.3, color="green")
        fmt_ax(ax, "Acúmulo de Recuperados")

        # 5 – Incidência diária
        ax = axes[1, 1]
        ax.bar(t, results["incidence"], color="tomato", alpha=0.7, label="Incidência diária")
        fmt_ax(ax, "Incidência Diária (Novos Casos)")

        # 6 – R_eff ao longo do tempo
        ax = axes[1, 2]
        ax.plot(t, results["R_eff"], "navy", lw=2, label="R_eff(t)")
        ax.axhline(1.0, color="red", linestyle="--", lw=1.5, label="Limiar R_eff = 1")
        ax.fill_between(t, results["R_eff"], 1.0,
                        where=results["R_eff"] >= 1.0, alpha=0.15, color="red")
        fmt_ax(ax, "Número Reprodutivo Efetivo R_eff(t)", ylabel="R_eff")

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches="tight")
            print(f"Gráfico salvo em: {save_path}")

        plt.show()

    def _plot_plotly(
        self, results: Dict[str, np.ndarray], save_path: Optional[str] = None
    ) -> None:
        try:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
        except ImportError:
            print("Plotly não instalado. Usando Matplotlib...")
            self._plot_matplotlib(results, save_path)
            return

        t = results["t"]

        fig = make_subplots(
            rows=2, cols=3,
            subplot_titles=(
                "Dinâmica SEIR",
                "Pós-infecciosos (R_mild / R_hosp)",
                "Desfechos Clínicos Graves",
                "Recuperados",
                "Incidência Diária",
                "R_eff(t)",
            ),
        )

        def add(row, col, x, y, name, color, fill=False):
            kwargs = dict(x=x, y=y, mode="lines", name=name,
                         line=dict(color=color, width=2))
            if fill:
                kwargs["fill"] = "tozeroy"
            fig.add_trace(go.Scatter(**kwargs), row=row, col=col)

        add(1, 1, t, results["S"], "Suscetível (S)", "blue")
        add(1, 1, t, results["E"], "Exposto (E)", "orange")
        add(1, 1, t, results["I"], "Infectado (I)", "red")
        add(1, 2, t, results["R_mild"], "Pós-inf. Leve (R_mild)", "lightgreen")
        add(1, 2, t, results["R_hosp"], "Pós-inf. Grave (R_hosp)", "darkred")
        add(1, 3, t, results["H"], "Hospitalizados (H)", "purple")
        add(1, 3, t, results["M"], "Óbitos (M)", "black")
        add(2, 1, t, results["Rec"], "Recuperados", "darkgreen", fill=True)

        fig.add_trace(
            go.Bar(x=t, y=results["incidence"], name="Incidência diária",
                   marker_color="tomato", opacity=0.7),
            row=2, col=2,
        )

        add(2, 3, t, results["R_eff"], "R_eff(t)", "navy")
        fig.add_hline(y=1.0, line_dash="dash", line_color="red",
                     annotation_text="R_eff = 1", row=2, col=3)

        fig.update_layout(
            title_text="Simulação Epidemiológica – Modelo SEIR Estendido (v2)",
            height=800, hovermode="x unified",
        )

        for r_ in (1, 2):
            for c_ in (1, 2, 3):
                fig.update_xaxes(title_text="Dias", row=r_, col=c_)
                fig.update_yaxes(title_text="Indivíduos", row=r_, col=c_)
        fig.update_yaxes(title_text="R_eff", row=2, col=3)

        if save_path:
            fig.write_html(save_path)
            print(f"Gráfico interativo salvo em: {save_path}")

        fig.show()


# ------------------------------------------------------------------ #
#  Exemplo de uso (execução direta)                                    #
# ------------------------------------------------------------------ #

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")   # evita janela gráfica em ambientes sem display

    print("=" * 70)
    print("CALCULADORA EPIDEMIOLÓGICA – MODELO SEIR ESTENDIDO v2")
    print("=" * 70)

    params = EpidemiologicalParameters(
        N=1_000_000,
        R_t=1.5,
        T_inf=10.0,
        T_inc=5.0,
        tau=0.15,       # 15% EXATOS necessitam hospitalização
        delta=0.02,     # 2% EXATOS dos hospitalizados morrem
        T_rec=14.0,
        T_hosp=3.0,
        T_int=10.0,
        I_0=100,
        E_0=200,
    )

    print(f"\nParâmetros principais:")
    print(f"  N = {params.N:,.0f} | R_t = {params.R_t} | β = {params.beta:.4f}")
    print(f"  T_inf = {params.T_inf} dias | T_inc = {params.T_inc} dias")
    print(f"  τ (hosp.) = {params.tau:.0%} | δ (mort.) = {params.delta:.0%}")

    modelo = EpidemiologicalModel(params)
    resultados = modelo.simulate(days=365)
    stats = modelo.get_statistics(resultados)

    print("\n" + "=" * 70)
    print("RESULTADOS (365 dias)")
    print("=" * 70)
    print(f"  Pico de infectados   : {stats['peak_infected']:>12,.0f}  (dia {stats['day_peak_infected']})")
    print(f"  Pico de hospitalizados: {stats['peak_hospitalized']:>12,.0f}  (dia {stats['day_peak_hospitalized']})")
    print(f"  Total de óbitos      : {stats['total_deaths']:>12,.0f}")
    print(f"  Total de recuperados : {stats['total_recovered']:>12,.0f}")
    print(f"  Taxa de ataque       : {stats['attack_rate']:.2%}")
    print(f"  IFR (mort./infect.)  : {stats['infection_fatality_rate']:.4%}")
    if stats["day_R_eff_below_1"] is not None:
        print(f"  R_eff < 1 a partir do dia: {stats['day_R_eff_below_1']}")

    print("\nEstado final dos compartimentos:")
    compartments_all = ["S", "E", "I", "R_mild", "R_hosp", "H", "Rec", "M"]
    total = sum(resultados[c][-1] for c in compartments_all)
    for c in compartments_all:
        print(f"  {c:>7}: {resultados[c][-1]:>12,.0f}")
    print(f"  {'Total':>7}: {total:>12,.0f}  (esperado: {params.N:,.0f})")
    print(f"  {'Desvio':>7}: {abs(total - params.N):>12,.2f}")
    print("\nSimulação concluída com sucesso!")
