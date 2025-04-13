from math import cos,sin,sqrt,pi
import opensees as ops
#
# MATERIAL PROPERTIES -------------------------------------------------------------------------------
#
# CONCRETE PROPERTIES
fc = 3.5*ksi;                                      # Default (general, if, unspecified) 28-day, concrete, cylinder, compressive, strength   (+Tension, -Compression)
fce = 5.0*ksi;                                                           # Default (general, if, unspecified) Unconfined, compressive, strength (max[5ksi, or, 1.3f'c])   (+Tension, -Compression)
ec0 = 0.002;                                                                         # Unconfined, strain, at, maximum, strength (+Tension, -Compression)
esp = 0.005;                                                                         # Unconfined, crushing (cover, spalling) strain (+Tension, -Compression)
Ec = 57000.0*sqrt(fce*ksi_psi)/ksi_psi;          # Default (general, if, unspecified) Concrete, Modulus, of, Elasticity
Uc = 0.2;                                             # Poisson's, ratio
Gc = Ec/(2.0*(1.0+Uc));                      # Shear, Modulus, of, Elasticity
wconc = 143.96*lb/ft**3;                        # Normal, Concrete, Weight, per, Volume                              
mconc = (143.96*lb/ft**3)/g;                   # Normal, Concrete, Mass, Weight, per, Volume
# REINFORCING STEEL PROPERTIES
Es = 29000.0*ksi;                            # Steel, Tensile, Modulus, of, Elasticity (initial, elastic, tangent)
Esh = 0.02*Es;                                   # Tangent, at, initial, strain, hardening
fy = 68.0*ksi;                       # Yield, Strength
fu = 95.0*ksi;                                   # Ultimate, strength
esh = 0.0075;                                        # Strain, corresponding, to, initial, strain, hardening 
esu = 0.090;                                         # Strain, at, peak, stress
