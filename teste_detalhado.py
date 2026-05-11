from Calculadora import EpidemiologicalModel, EpidemiologicalParameters
import numpy as np

# Teste rápido
m = EpidemiologicalModel(EpidemiologicalParameters())
r = m.simulate(100)

print(f'Pico de infectados em 100 dias: {max(r["I"]):,.0f}')
print()

# Verificar conservação
compartments = ['S', 'E', 'I', 'R', 'H', 'Rec', 'M']
totals = np.sum([r[comp] for comp in compartments], axis=0)
print(f"População inicial: {totals[0]:,.0f}")
print(f"População final: {totals[-1]:,.0f}")
print(f"População esperada: {m.params.N:,}")
print(f"Diferença: {totals[-1] - m.params.N:,.0f}")
print()

# Estado final
print("Estado final em dia 100:")
for comp in compartments:
    print(f"  {comp}: {r[comp][-1]:,.0f}")
print(f"  Total: {totals[-1]:,.0f}")
