from numpy import linspace, sin
from opensees import uniaxial

fy = 50e3
E  = 29e6

b  = 0.005

R0 = 18
cR1 = 0.915
cR2 = 0.15

# eh

strain = fy/E*sin(linspace(0, 20, 100))

mats = [
    uniaxial.Steel02( fy, E),
    uniaxial.Steel02( fy, E, b),
    uniaxial.Steel02( fy, E, b, R0),
    uniaxial.Steel02( fy, E, b, R0, 0.925, 0.15)
]



kin = [b, R0, cR1, cR2]
mats = [
    *mats,
#   uniaxial.Steel04(None, fy, E), #, iso, ult, sig_init), 
#   uniaxial.Steel04(None, fy, E, kin), #, iso, ult, sig_init), 
#   uniaxial.Steel04(None, fy, E, kin), #, iso, ult, sig_init), 
]

def test(mat):
    with mat as m:
        print("\t", [m.getStress(e, True) for e in strain])

for m in mats: test(m)
