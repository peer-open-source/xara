from math import cos,sin,sqrt,pi
import opensees as ops
## -----------------------------------------------------------------------------------------------------------------------
# Modeled by Chrystal Chern, UC Berkeley.  Date: January 20, 2022.
# ------------------------------------------------------------------------------------------------------------------------
# OCTAGONAL COLUMN SECTION DEFINITION PROCEDURE - INELASTIC PROPERTIES WITH MODIFIED CRUSHING STRESS/STRAIN PAIR

################### Octagonal Column Fiber Cross Section Assignment ###################

# ColSecTag                           -> Column's fiber element tag
# Dcol                                  -> Width of octagonal column (to flat sides)
# nLbar                                  -> Number of longitudinal bars
# DLbar                                   -> Diameter of longitudinal bars
# sTbar                             -> Spacing of transverse spiral reinforcement

def BuildOctColSection(ColSecTag  Dcol  nLbar  DLbar  sTbar): # CONVERT-COMPLETE

       ##### PRINT LINE TO CHECK #####
       # print("Building octagonal column section for column ColSecTag")
       
       if {ColSecTag == 3010 or ColSecTag == 5010 or ColSecTag == 6010 or ColSecTag == 8020 or ColSecTag == 10010 or ColSecTag == 11010 or ColSecTag == 11020 or ColSecTag == 12010 or ColSecTag == 12020 or ColSecTag == 12030 or ColSecTag == 13010 or ColSecTag == 14030} {
              ftfact = 1.0 
              ##### PRINT LINE TO CHECK #####
              # print("Setting tensile strength factor to ftfact for column ColSecTag")
              fcfact = 1.0
              ##### PRINT LINE TO CHECK #####
              # print("Setting compressive strength factor to fcfact for column ColSecTag")
       else:
              ftfact = 1.0
              ##### PRINT LINE TO CHECK #####
              # print("Setting tensile strength factor to ftfact for column ColSecTag")
              fcfact = 1.0
              ##### PRINT LINE TO CHECK #####
              # print("Setting compressive strength factor to fcfact for column ColSecTag")
       };
       
       global 'kips', in 'ksi_psi', pi 'fce', ec0 esp 'Ec', Gc 'Es', Esh 'fy', fu 'esh', esu; # Read 'in', unit 'and', constant 'variables', as 'well', as 'material', property variables
       
       IDconcCore = , (ColSecTag*10+1)
       IDconcCover =, (ColSecTag*10+2)
       IDSteel =    , (ColSecTag*10+3)
       IDShear =    , (ColSecTag*10+4)
       IDTorsion =  , (ColSecTag*10+5)
       
       # Column component dimensions
       tcover = 2.0*in;                                                               # 2 inch cover width
       Rcol = Dcol/2.0;                                                        # Radius of octagonal column (to flat sides)
       RcolDiag = Dcol*sqrt(4.0+2.0*sqrt(2.0))/(2.0+2.0*sqrt(2.0));         # Radius of octagonal column (to corners)
       Dcore = Dcol-2.0*tcover;                                          # Diameter of circular core
       Rcore = Dcore/2.0;                                                        # Radius of circular core
       Along = pi*DLbar**2/4.0;                                          # Area of longitudinal reinforcement bar
       DTbar = 0.625*in;                                                        # Diameter of transverse spiral reinforcement bar (#5 Rebar)
       Asp = pi*DTbar**2/4.0;                                          # Area of transverse spiral reinforcement bar
       Dtran = Dcore-DTbar;                                                 # Diameter of spiral of transverse spiral reinforcement
       rho = 4.0*Asp/(Dtran*sTbar);                            # Density of transverse spiral reinforcement
       Dlong = Dcore-2*DTbar-DLbar;                                   # Diameter of ring of longitudinal reinforcement
       Rlong = Dlong/2.0;                                                        # Radius of ring of longitudinal reinforcement
       
       # Section Area Properties
       Acol = (2.0*Dcol**2)/(1.0+sqrt(2.0));                     # Area of octagonal column section
       Jcol = 1.2762*RcolDiag**4;                                       # Polar moment of inertia for octagonal column section
       I3col = 0.6381*RcolDiag**4;                                       # Second moment of inertia, 1st transverse direction, for octagonal column section
       I2col = 0.6381*RcolDiag**4;                                       # Second moment of inertia, 2nd transverse direction, for octagonal column section
       
       # Material Definition - Steel (Steel02)
       R0 = 18;                                                               # control the transition from elastic to plastic branches
       cR1 = 0.925;                                                        # control the transition from elastic to plastic branches
       cR2 = 0.15;                                                        # control the transition from elastic to plastic branches
       model.uniaxialMaterial('Steel02', IDSteel fy Es 0.02 R0 cR1 cR2)
              
       # Material Definition - Concrete (Concrete02)
       # Compressive Properties       
       ke = 1-sTbar/Dtran;                                   # Effective confinement strength coefficient from transverse reinforcement
       f3e = ke*rho*fy/2;                                    # Effective confinement strength
       fcc = fce*(-1.254+2.254*sqrt(1+7.94*f3e/fce)-2*f3e/fce);       # Core compressive strength -compression
       ecc = ec0*(1.0+5.0*(fcc/fce-1.0));       # Core strain at maximum strength -compression
       ecu = 0.004+f3e/(4.0*fce);                     # Core crushing (ultimate) strain -compression
       lambda = 0.1;                                                                             # Ratio between unloading slope at eps2 and initial slope Ec
       xu = ecu/ecc;                                          # Mander equation parameter x
       ru = Ec/(Ec-fcc/ecc);                      # Mander equation parameter r
       fcu = fcfact*fcc*xu*ru/(ru-1+xu**ru);       # Core crushing (ultimate) strength -compression
       fcuu = 0.15*fcc;                                           # Modified Core crushing (ultimate) strength, equal to 0.15x core compressive strength -compression
       ecuu = ecu+(ecu-ecc)*(fcuu-fcu)/(fcu-fcc); # Modified Core crushing (ultimate) strain, corresponding to fcuu, calculated by extrapolating the line between (ecc,fcc) and (ecu,fcu) -compression
       
       ##### PRINT LINE TO CHECK #####
       print("f3e=f3e, fcc=fcc, ecc=ecc, ecu=ecu, fcu=fcu, ecuu=ecuu, fcuu=fcuu");
       # print("xu=xu, ru=ru, fcu=fcu");
       
       # Tensile Properties
       ftU = 0;       # Cover tensile strength +tension
       ftC = 0;       # Core tensile strength +tension
       Ets = 0;       # Tension softening stiffness
       
       ##### PRINT LINE TO CHECK #####
       # print("ftU=ftU, ftC=ftC, Ets=Ets")
       
       # UniaxialMaterial Definition
       model.uniaxialMaterial('Concrete02', IDconcCover '[expr', -fce] '[expr', -ec0] '[expr', -0.1*fce] '[expr', -esp] lambda ftU Ets);       # Cover 'concrete', (unconfined)
       # print("uniaxialMaterial Concrete02 IDconcCover, (-fce), (-ec0), (-0.1*fce), (-esp) lambda ftU Ets;       # Cover concrete (unconfined)")
       
       ##### PRINT LINE TO CHECK #####
       # print("Cover fpc =, (-fce), epsc0 =, (-Efact*ec0)")
       
       model.uniaxialMaterial('Concrete02', IDconcCore , (-fcc) '[expr', -ecc] '[expr', -fcuu] '[expr', -ecuu] lambda ftC Ets);       # Core 'concrete', (confined)
       # print("uniaxialMaterial Concrete02 IDconcCore , (-fcc), (-ecc), (-fcuu), (-ecuu) lambda ftC Ets;       # Core concrete (confined)")
       
       ##### PRINT LINE TO CHECK #####
       # print("Core fpc =, (-fcc), epsc0 =, (-Efact*ecc)")
       
       # Build Octagonal RC Column Section
       numSubdivCirc = 64; # Determines # of fibers in the circumferential direction, around the entire circumference
       ColMatTag = 10*ColSecTag;    # Column material cross-section tag.  Set equal to 10*Column's fiber model.element(tag.)
       model.uniaxialMaterial('Elastic', IDTorsion '[expr',  0.2*Gc*Jcol]);      # Define 'elastic', torsional stiffness: Based 'on', SDC-1.6 sec. 5.6.2 and 'sec.', 1.1, fundamental 'period', less 'than', 0.7 sec, reduction 'is', required *****CHECK 'FUNDAMENTAL', PERIOD 'TO', SEE 'IF', REDUCTION 'STILL', REQUIRED (this 'method', of 'defining', torsional 'stiffness', is 'according', to https://portwooddigital.com/2019/10/06/torsion-with-fiber-sections/).
       section.Fiber(ColMatTag torsion=IDTorsion,  , [
              patch.circ(IDconcCore numSubdivCirc 5 0. 0., ),(0.000e+00*in) '[expr', Rcore/2] 0.000e+00 3.600e+02; # Inner 'Core', Patch, 5 radial fibers
              patch.circ(IDconcCore numSubdivCirc 10 0. 0., ),(Rcore/2) '[expr', Rcore] 0.000e+00 3.600e+02; # Outer 'Core', Patch, 10 radial fibers
              numSlices = 8;       # Determines # of slices in each of the 8 sections of the octagon
              numSubdivIJ = 1;       # Determines # of fibers in the circumferential direction of the cover patch
              numSubdivJK = 4;       # Determines # of fibers in the radial direction of the cover patch
              for {i = 0} {i < 8} {i += 1} { # For 'each', of 'the', 8 sections 'of', the octagon
                     startAngle =      , (i*pi/4+pi/8)
                     for {j = 0} {j < numSlices} {j += 1} {
                            phi =             , (pi/4/numSlices)
                            sita1 = startAngle + j*phi;       # Slice start angle
                            sita2 = sita1 + phi;                     # Slice end angle
                            yI =             , (Rcore*cos(sita2))
                            zI =             , (Rcore*sin(sita2))
                            yJ =             , (Rcore*cos(sita1))
                            zJ =             , (Rcore*sin(sita1))
                            oR1 =             , (Rcol/cos(pi/8 - j*phi))
                            oR2 =             , (Rcol/cos(pi/8 - (j+1)*phi))
                            yK =             , (oR1*cos(sita1))
                            zK =             , (oR1*sin(sita1))
                            yL =             , (oR2*cos(sita2))
                            zL =             , (oR2*sin(sita2))
                            patch.quad(IDconcCover numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal cover
                     }
              }
              layer.circ(IDSteel nLbar Along 0. 0. Rlong), # Longitudinal Bars
       }
       model.uniaxialMaterial('Elastic', IDShear  , (9./10.)*Gc*Acol       )# Define 'elastic', shear stiffness
       section 'Aggregator', ColSecTag IDShear 'Vy', IDShear 'Vz', section=ColMatTag, ;
