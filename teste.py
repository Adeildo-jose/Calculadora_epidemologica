from Calculadora import EpidemiologicalModel, EpidemiologicalParameters

# Teste rápido
m = EpidemiologicalModel(EpidemiologicalParameters())
r = m.simulate(100)
print(f'Pico: {max(r["I"]):,.0f}')
