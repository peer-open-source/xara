from math import cos,sin,sqrt,pi
import opensees as ops
#
# IN-SPAN HINGE SPRINGS -------------------------------------------------------------------------------
#
khlc = 25000.0*kips/in;                       # Pounding, stiffness, is, equal, for, all, hinges', longitudinal, compressive, stiffness
phl = 1.0e+9*kips;                          # Pounding, strength, and, tensile, strength, are, both, a high, number, for, all, hinges', longitudinal, strength
Eb = 0.50*ksi;                                    # Elastic, modulus, for, elastomeric, bearings
gapV = 0.6*in;                                     # The, plastic, portion, of, the, elastomeric, bearing, pad, which, is, flexible
# In span hinge 5 Properties
gaphl5 = 3.50*in;                                    # Longitudinal, gap, width, Hinge, 5  
alr5 = 24.54*in**2;                             # Longitudinal, restrainer, area, Hinge, 5
llr5 = 506.0*in;                                    # Longitudinal, restrainer, length, Hinge, 5
khlt5 = Es*alr5/llr5;                      # Longitudinal, tensile, stiffness, Hinge, 5
kht5 = 41422.0*kips/in;                      # Transverse, stiffness, Hinge, 5
pht5 = 1036.0*kips;                             # Transverse, strength, Hinge, 5
lb5 = 4.5*in;                                    # Hinge, 5 elastomeric, bearing, thickness
ab5 = 16.0*22.0*in**2;                     # Hinge, 5 elastomeric, bearing, area
nb5 = 6.0;                                           # Hinge, 5 number, of, elastomeric, bearings
khv5 = nb5*Eb*ab5/lb5;                      # Hinge, 5 vertical, initial, stiffness (in, kip/in)
# In span hinge 8 Properties
gaphl8 = 2.75*in;                                    # Longitudinal, gap, width, Hinge, 8
alr8 = 24.54*in**2;                             # Longitudinal, restrainer, area, Hinge, 8
llr8 = 458.0*in;                                    # Longitudinal, restrainer, length, Hinge, 8
khlt8 = Es*alr8/llr8;                      # Longitudinal, tensile, stiffness, Hinge, 8
kht8 = 41422.0*kips/in;                      # Transverse, stiffness, Hinge, 8
pht8 = 1036.0*kips;                             # Transverse, strength, Hinge, 8
lb8 = 3.5*in;                                    # Hinge, 8 elastomeric, bearing, thickness
ab8 = 12.0*22.0*in**2;                     # Hinge, 8 elastomeric, bearing, area
nb8 = 6.0;                                           # Hinge, 8 number, of, elastomeric, bearings
khv8 = nb8*Eb*ab8/lb8;                      # Hinge, 8 vertical, initial, stiffness (in, kip/in)
# In span hinge 12 NE Properties
gaphl12NE = 2.25*in;                                    # Longitudinal, gap, width, Hinge, 12, NE
alr12NE = 24.54*in**2;                             # Longitudinal, restrainer, area, Hinge, 12, NE
llr12NE = 326.0*in;                                    # Longitudinal, restrainer, length, Hinge, 12, NE
khlt12NE = Es*alr12NE/llr12NE;               # Longitudinal, tensile, stiffness, Hinge, 12, NE
kht12NE = 41422.0*kips/in;                      # Transverse, stiffness, Hinge, 12, NE
pht12NE = 1036.0*kips;                             # Transverse, strength, Hinge, 12, NE
lb12NE = 3.0*in;                                    # Hinge, 12, NE, elastomeric, bearing, thickness
ab12NE = 14.0*22.0*in**2;                     # Hinge, 12, NE, elastomeric, bearing, area
nb12NE = 6.0;                                           # Hinge, 12, NE, number, of, elastomeric, bearings
khv12NE = nb12NE*Eb*ab12NE/lb12NE; # Hinge, 12, NE, vertical, initial, stiffness (in, kip/in)
# In span hinge 12 NR Properties
gaphl12NR = 1.75*in;                                    # Longitudinal, gap, width, Hinge, 12, NR
alr12NR = 9.82*in**2;                             # Longitudinal, restrainer, area, Hinge, 12, NR
llr12NR = 330.0*in;                                    # Longitudinal, restrainer, length, Hinge, 12, NR
khlt12NR = Es*alr12NR/llr12NR;               # Longitudinal, tensile, stiffness, Hinge, 12, NR
kht12NR = 34518.0*kips/in;                      # Transverse, stiffness, Hinge, 12, NR
pht12NR = 863.0*kips;                             # Transverse, strength, Hinge, 12, NR
lb12NR = 2.5*in;                                    # Hinge, 12, NR, elastomeric, bearing, thickness
ab12NR = 12.0*18.0*in**2;                     # Hinge, 12, NR, elastomeric, bearing, area
nb12NR = 3.0;                                           # Hinge, 12, NR, number, of, elastomeric, bearings
khv12NR = nb12NR*Eb*ab12NR/lb12NR; # Hinge, 12, NR, vertical, initial, stiffness (in, kip/in)
# print("khlc = khlc; phl = phl")
# print("khlt5 = khlt5; khlt8 = khlt8; khlt12NE = khlt12NE; khlt12NR = khlt12NR")
# print("kht5 = kht5; kht8 = kht8; kht12NE = kht12NE, kht12NR = kht12NR")
# print("pht5 = pht5; pht8 = pht8; pht12NE = pht12NE, pht12NR = pht12NR")
# print("khv5 = khv5; khv8 = khv8; khv12NE = khv12NE, khv12NR = khv12NR")
# Hinge 5 Spring Materials
model.uniaxialMaterial('ElasticPPGap',        5011       , (khlc) '[expr', -phl], '[expr', -gaphl5], 1.0e-3, damage);                                                  # LONGITUDINAL, 'Hinge', 5, COMPRESSION
model.uniaxialMaterial('ElasticPPGap',        5012       , (khlt5) '[expr', phl], '[expr', gaphl5], 1.0e-3, damage);                                                  # LONGITUDINAL, 'Hinge', 5, TENSION
model.uniaxialMaterial('Parallel',               501, 5011, 5012);                                                                                                                                      # LONGITUDINAL, 'Hinge', 5
model.uniaxialMaterial('ElasticPP',               502, kht5, pht5/kht5                                                                                                          )# TRANSVERSE, 'Hinge', 5
# uniaxialMaterial Elastic               503       , 1.0e+9*kips/in                                                                                                          # VERTICAL Hinge 5
model.uniaxialMaterial('ENT',                      5031       , 3.5*khv5)
model.uniaxialMaterial('ElasticPPGap',        5032       , (1.0e)+9) '[expr', -1.0e+9] -gapV, 1.0e-3, damage;
model.uniaxialMaterial('Parallel',               503, 5031, 5032);                                                                                                                                      # VERTICAL, 'Hinge', 5
# Hinge 8 Spring Materials
model.uniaxialMaterial('ElasticPPGap',        8011       , (khlc) '[expr', -phl], '[expr', -gaphl8], 1.0e-3, damage);                                                  # Hinge, 8 longitudinal, 'COMPRESSION'
model.uniaxialMaterial('ElasticPPGap',        8012       , (khlt8) '[expr', phl], '[expr', gaphl8], 1.0e-3, damage);                                                  # Hinge, 8 longitudinal, 'TENSION'
model.uniaxialMaterial('Parallel',               801, 8011, 8012);                                                                                                                                      # LONGITUDINAL, 'Hinge', 8
model.uniaxialMaterial('ElasticPP',               802, kht8, pht8/kht8                                                                                                          )# TRANSVERSE, 'Hinge', 8
# uniaxialMaterial Elastic               803       , 1.0e+9*kips/in                                                                                                          # VERTICAL Hinge 8
model.uniaxialMaterial('ENT',                      8031       , 3.5*khv8)
model.uniaxialMaterial('ElasticPPGap',        8032       , (1.0e)+9) '[expr', -1.0e+9] -gapV, 1.0e-3, damage;
model.uniaxialMaterial('Parallel',               803, 8031, 8032);                                                                                                                                      # VERTICAL, 'Hinge', 8
# Hinge 12 NE Spring Materials
model.uniaxialMaterial('ElasticPPGap',        12011       , (khlc) '[expr', -phl], '[expr', -gaphl12NE], 1.0e-3, damage);                                           # Hinge, 12, NE, 'longitudinal', 'COMPRESSION'
model.uniaxialMaterial('ElasticPPGap',        12012       , (khlt12NE) '[expr', phl], '[expr', gaphl12NE], 1.0e-3, damage);                                           # Hinge, 12, NE, 'longitudinal', 'TENSION'
model.uniaxialMaterial('Parallel',               1201, 12011, 12012);                                                                                                                               # LONGITUDINAL, 'Hinge', 12, NE
model.uniaxialMaterial('ElasticPP',               1202, kht12NE, '[expr', pht12NE/kht12NE]);                                                                                            # TRANSVERSE, 'Hinge', 12, NE
# uniaxialMaterial Elastic               1203       , 1.0e+9*kips/in                                                                                                          # VERTICAL Hinge 12 NE
model.uniaxialMaterial('ENT',                      12031       , 3.5*khv12NE)
model.uniaxialMaterial('ElasticPPGap',        12032       , (1.0e)+9) '[expr', -1.0e+9] -gapV, 1.0e-3, damage;
model.uniaxialMaterial('Parallel',               1203, 12031, 12032);                                                                                                                               # VERTICAL, 'Hinge', 12, NE
# Hinge 12 NR Spring Materials
model.uniaxialMaterial('ElasticPPGap',        12041       , (khlc) '[expr', -phl], '[expr', -gaphl12NR], 1.0e-3, damage);                                           # Hinge, 12, NR, 'longitudinal', 'COMPRESSION'
model.uniaxialMaterial('ElasticPPGap',        12042       , (khlt12NR) '[expr', phl], '[expr', gaphl12NR], 1.0e-3, damage);                                           # Hinge, 12, NR, 'longitudinal', 'TENSION'
model.uniaxialMaterial('Parallel',               1204, 12011, 12012);                                                                                                                               # LONGITUDINAL, 'Hinge', 12, NR
model.uniaxialMaterial('ElasticPP',               1205, kht12NR, '[expr', pht12NR/kht12NR]);                                                                                            # TRANSVERSE, 'Hinge', 12, NR
# uniaxialMaterial Elastic               1206       , 1.0e+9*kips/in                                                                                                          # VERTICAL Hinge 12 NR
model.uniaxialMaterial('ENT',                      12061       , 3.5*khv12NR)
model.uniaxialMaterial('ElasticPPGap',        12062       , (1.0e)+9) '[expr', -1.0e+9] -gapV, 1.0e-3, damage;
model.uniaxialMaterial('Parallel',               1206, 12061, 12062);                                                                                                                               # VERTICAL, 'Hinge', 12, NR
# Spring Material for all hinges
model.uniaxialMaterial('Elastic',               504       , 0.0*kips/in                                                                                                                 )# ROTATIONAL, 'stiffness', for, 'all', hinges, 'set', to, 'zero'
# ZeroLength elements defined for each hinge
model.element('zeroLength',                             500000, 50001, 500010 -mat, 501, 502, 503, 504, 504, 504 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag -orient, 456.63 -162.15534, 0 0, 1 0); # IN-SPAN, 'HINGE', BENT, 5
model.element('zeroLength',                             800000, 80004, 800040 -mat, 801, 802, 803, 504, 504, 504 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag -orient, 702.52 -66.80943, 0 0, 1 0); # IN-SPAN, 'HINGE', BENT, 8
model.element('zeroLength',                             1200000, 120001, 1200010 -mat, 1201, 1202, 1203, 504, 504, 504 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag -orient, 353.11, 41.662, 0 0, 1 0); # IN-SPAN, 'HINGE', BENT, 12, NE
model.element('zeroLength',                             1200001, 120005, 1200050 -mat, 1204, 1205, 1206, 504, 504, 504 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag -orient, 195.72 -6.42, 0 0, 1 0); # IN-SPAN, 'HINGE', BENT, 12, NR

