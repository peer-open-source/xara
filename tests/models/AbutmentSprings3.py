from math import cos,sin,sqrt,pi
import opensees as ops
#
# ABUTMENT SPRINGS -------------------------------------------------------------------------------
#
CW = 4.0/3.0;               # Wall, participation, coefficient (for, wingwall, contribution, to, transverse, stiffness)
CL = 2.0/3.0;               # Wall, effectiveness (for, wingwall, contribution, to, transverse, stiffness)
CWT = 1.0/3.0;                # Multiplier, for, assumed, effective, wingwall, width = 1/3, of, backwall, width
Eb = 0.50*ksi;               # Elastic, modulus, for, elastomeric, bearings
gapV = 0.6*in;                # The, plastic, portion, of, the, elastomeric, bearing, pad, which, is, flexible
# Abutment 1 Properties
wabut1 = 54.0;                                                                              # Abutment, 1 backwall, width (in, feet)
habut1 = 7.0;                                                                              # Abutment, 1 backwall, height (in, feet)
theta1 = 0.25;                                                                              # Abutment, 1 skew, angle (in, degrees)
Rsk1 = exp(-theta1/45);                                                         # Abutment, 1 skew, reduction, factor
kabutl1 = wabut1*(5.5*habut1+20.0)*Rsk1*kips/in;               # Abutment, 1 longitudinal, stiffness (in, kip/in)
kabutt1 = CW*CL*CWT*kabutl1*kips/in;                                    # Abutment, 1 transverse, stiffness (in, kip/in)
pabutl1 = wabut1*(5.5*habut1**2.5)*Rsk1/(1+2.37*habut1);# Abutment, 1 longitudinal, resistance, force (in, kips)
pabutt1 = CW*CL*CWT*pabutl1*kips;                                    # Abutment, 1 transverse, resistance, force (in, kips)
gap1 = 1.75*in;                                                                        # Abutment, 1 longitudinal, gap (in, inches)
lb1 = 2.0*in;                                                                       # Abutment, 1 elastomeric, bearing, thickness
ab1 = 12.0*18.0*in**2;                                                        # Abutment, 1 elastomeric, bearing, area
nb1 = 6.0;                                                                              # Abutment, 1 number, of, elastomeric, bearings
kabutv1 = nb1*Eb*ab1/lb1;                                                         # Abutment, 1 vertical, initial, stiffness (in, kip/in)
# Abutment 15 NE Properties
wabut15NE = 75.5;                                                                              # Abutment, 15, NE, width (in, feet)
habut15NE = 5.5;                                                                              # Abutment, 15, NE, backwall, height (in, feet)
theta15NE = 44.34;                                                                              # Abutment, 15, NE, skew, angle (in, degrees)
Rsk15NE = exp(-theta15NE/45);                                                  # Abutment, 15, NE, skew, reduction, factor
kabutl15NE = wabut15NE*(5.5*habut15NE+20.0)*Rsk15NE*kips/in;# Abutment, 15, NE, longitudinal, stiffness (in, kip/in)
kabutt15NE = CW*CL*CWT*kabutl15NE*kips/in;                             # Abutment, 15, NE, transverse, stiffness (in, kip/in)
pabutl15NE = wabut15NE*(5.5*habut15NE**2.5)*Rsk15NE/(1+2.37*habut15NE);# Abutment, 15, NE, longitudinal, resistance, force (in, kips)
pabutt15NE = CW*CL*CWT*pabutl15NE*kips;                                    # Abutment, 15, NE, transverse, resistance, force (in, kips)
gap15NE = 1.0*in;                                                                        # Abutment, 15, NE, longitudinal, gap (in, inches)
lb15NE = 1.50*in;                                                                       # Abutment, 15, NE, elastomeric, bearing, thickness
ab15NE = 16.0*20.0*in**2;                                                        # Abutment, 15, NE, elastomeric, bearing, area
nb15NE = 6.0;                                                                              # Abutment, 15, NE, number, of, elastomeric, bearings
kabutv15NE = nb15NE*Eb*ab15NE/lb15NE;                                    # Abutment, 15, NE, vertical, initial, stiffness (in, kip/in)
# Abutment 15 NR Properties
wabut15NR = 32.1;                                                                              # Abutment, 15, NR, width (in, feet)
habut15NR = 5.5;                                                                              # Abutment, 15, NR, backwall, height (in, feet)
theta15NR = 35.896;                                                                              # Abutment, 15, NR, skew, angle (in, degrees)
Rsk15NR = exp(-theta15NR/45);                                                  # Abutment, 15, NR, skew, reduction, factor
kabutl15NR = wabut15NR*(5.5*habut15NR+20.0)*Rsk15NR;               # Abutment, 15, NR, longitudinal, stiffness (in, kip/in)
kabutt15NR = CW*CL*CWT*kabutl15NR*kips/in;                             # Abutment, 15, NR, transverse, stiffness (in, kip/in)
pabutl15NR = wabut15NR*(5.5*habut15NR**2.5)*Rsk15NR/(1+2.37*habut15NR);# Abutment, 15, NR, longitudinal, resistance, force (in, kips)
pabutt15NR = CW*CL*CWT*pabutl15NR*kips;                                    # Abutment, 15, NR, transverse, resistance, force (in, kips)
gap15NR = 1.0*in;                                                                        # Abutment, 15, NR, longitudinal, gap (in, inches)
lb15NR = 1.50*in;                                                                       # Abutment, 15, NR, elastomeric, bearing, thickness
ab15NR = 16.0*20.0*in**2;                                                        # Abutment, 15, NR, elastomeric, bearing, area
nb15NR = 3.0;                                                                              # Abutment, 15, NR, number, of, elastomeric, bearings
kabutv15NR = nb15NR*Eb*ab15NR/lb15NR;                                    # Abutment, 15, NR, vertical, initial, stiffness (in, kip/in)
# print("kabutl1 = kabutl1; kabutl15NE = kabutl15NE; kabutl15NR = kabutl15NR")
# print("kabutt1 = kabutt1; kabutt15NE = kabutt15NE; kabutt15NR = kabutt15NR")
# print("kabutv1 = kabutv1; kabutv15NE = kabutv15NE; kabutv15NR = kabutv15NR")
# print("pabutl1 = pabutl1; pabutl15NE = pabutl15NE, pabutl15NR = pabutl15NR")
# print("pabutt1 = pabutt1, pabutt15NE = pabutt15NE, pabutt15NR = pabutt15NR")
# Longitudinal (local x) stiffness is working only in compression
# Transverse (local y) stiffness working in tension and compression, because they work parallel with half of SDC stiffness assigned for each
# Vertical (local z) stiffness has two portions constructed by two parallel uniaxialMaterials 
cfactor =                 0.5;
# Abutment 1 Spring Materials
model.uniaxialMaterial('Elastic',                201       , kabutl1*cfactor        )# LONGITUDINAL, 'Abutment', 1
model.uniaxialMaterial('Elastic',               202       , kabutt1*cfactor        )# TRANSVERSE, 'Abutment', 1
model.uniaxialMaterial('Elastic',               203       , 3.5*kabutv1               )# VERTICAL, 'Abutment', 1
# uniaxialMaterial ElasticPPGap        201       , (kabutl1*cfactor), (-pabutl1*cfactor), (-gap1) 1.0e-3 damage;        # LONGITUDINAL Abutment 1
# uniaxialMaterial ElasticPP               202       , (kabutt1*cfactor), pabutt1/kabutt1                                                                # TRANSVERSE Abutment 1
# uniaxialMaterial ENT                      2031       , 3.5*kabutv1
# uniaxialMaterial ElasticPPGap        2032       , (1.0e+9) [expr -1.0e+9] -gapV 1.0e-3 damage;
# uniaxialMaterial Parallel               203        2031 2032;                                                                                                                                      # VERTICAL Abutment 1
# Abutment 15 NE Spring Materials
model.uniaxialMaterial('Elastic',                215       , kabutl15NE*cfactor        )# LONGITUDINAL, 'Abutment', 15, NE
model.uniaxialMaterial('Elastic',               216       , kabutt15NE*cfactor        )# TRANSVERSE, 'Abutment', 15, NE
model.uniaxialMaterial('Elastic',               217       , 3.5*kabutv15NE               )# VERTICAL, 'Abutment', 15, NE
# uniaxialMaterial ElasticPPGap        215       , (kabutl15NE*cfactor), (-pabutl15NE*cfactor), (-gap15NE) 1.0e-3 damage;        # LONGITUDINAL Abutment 15 NE
# uniaxialMaterial ElasticPP               216       , (kabutt15NE*cfactor), pabutt15NE/kabutt15NE                                                         # TRANSVERSE Abutment 15 NE
# uniaxialMaterial ENT                      2171       , 3.5*kabutv15NE
# uniaxialMaterial ElasticPPGap        2172       , (1.0e+9) [expr -1.0e+9] -gapV 1.0e-3 damage;
# uniaxialMaterial Parallel               217        2171 2172;                                                                                                                                                    # VERTICAL Abutment 15 NE
# Abutment 15 NR Spring Materials
model.uniaxialMaterial('Elastic',                218       , kabutl15NR*cfactor        )# LONGITUDINAL, 'Abutment', 15, NR
model.uniaxialMaterial('Elastic',               219       , kabutt15NR*cfactor        )# TRANSVERSE, 'Abutment', 15, NR
model.uniaxialMaterial('Elastic',               220       , 3.5*kabutv15NR               )# VERTICAL, 'Abutment', 15, NR
# uniaxialMaterial ElasticPPGap        218       , (kabutl15NR*cfactor), (-pabutl15NR*cfactor), (-gap15NR) 1.0e-3 damage;        # LONGITUDINAL Abutment 15 NR
# uniaxialMaterial ElasticPP               219       , (kabutt15NR*cfactor), pabutt15NR/kabutt15NR                                                         # TRANSVERSE Abutment 15 NR
# uniaxialMaterial ENT                      2201       , 3.5*kabutv15NR
# uniaxialMaterial ElasticPPGap        2202       , (1.0e+9) [expr -1.0e+9] -gapV 1.0e-3 damage;
# uniaxialMaterial Parallel               220        2201 2202;                                                                                                                                                    # VERTICAL Abutment 15 NR
# Spring Material for all abutments
model.uniaxialMaterial('Elastic',               2000       , 0.0*kips/in                                                  )# ROTATIONAL
# zerolength element defined from fixed end to the free end of the spring and the z axis defined by cross product of "x" and "yp"
# to have the rigid motion of the abutment the rotational degree of freedom of the zerolength element is released
# Note that the local x axis is parallel to the backwall (abutment element), which corresponds to the transverse direction. The local y axis is perpendicular to the backwall and corresponds to the longitudinal direction.
model.element('zeroLength', 120000, 1021, 1020, mat=202,  201, 203, 2000, 2000, 2000, dir=[1, 2 3, 4 5, 6 ], doRayleigh=rFlag,  orient=True, -156.24 -202.12, 0 176.95 -135.55, 0); 
model.element('zeroLength', 130000, 1031, 1030, mat=202,  201, 203, 2000, 2000, 2000, dir=[1, 2 3, 4 5, 6 ], doRayleigh=rFlag,  orient=156.24,  202.12, 0 -176.95, 135.55, 0); 
model.element('zeroLength', 1510000, 15021, 15020, mat=216,  215, 217, 2000, 2000, 2000, dir=[1, 2 3, 4 5, 6 ], doRayleigh=rFlag,  orient=True, -223.99 -290.74, 0 329.51, 38.87, 0); 
model.element('zeroLength', 1520000, 15031, 15030, mat=216,  215, 217, 2000, 2000, 2000, dir=[1, 2 3, 4 5, 6 ], doRayleigh=rFlag,  orient=223.99,  290.74, 0 -329.51 -38.87, 0); 
model.element('zeroLength', 1550000, 15051, 15050, mat=219,  218, 220, 2000, 2000, 2000, dir=[1, 2 3, 4 5, 6 ], doRayleigh=rFlag,  orient=True, -77.97 -101.20, 0 348.96 -10.50, 0); 
model.element('zeroLength', 1560000, 15061, 15060, mat=219,  218, 220, 2000, 2000, 2000, dir=[1, 2 3, 4 5, 6 ], doRayleigh=rFlag,  orient=77.97,  101.20, 0 -348.96, 10.50, 0);

