from opensees.lib import backbone
from numpy import linspace

fc = 4e3
Ec = 30e3
epsc = 0.003
strain = linspace(-0.003, 0.003, 100)
with backbone.Popovics(1, fc, epsc, Ec) as b:
    print([b.getStress(e) for e in strain])

