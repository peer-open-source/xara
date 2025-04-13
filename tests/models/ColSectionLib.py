
def BuildOctColSectionElastic(ColSecTag  Dcol  nLbar  DLbar  sTbar):
  # OCTAGONAL COLUMN SECTION DEFINITION PROCEDURE

  ################### Octagonal Column Fiber Cross Section Assignment ###################

      # ColSecTag        -> Column's fiber element tag
      # Dcol             -> Width of octagonal column (to flat sides)
      # nLbar            -> Number of longitudinal bars
      # DLbar            -> Diameter of longitudinal bars
      # sTbar            -> Spacing of transverse spiral reinforcement


       ##### PRINT LINE TO CHECK #####
       # print("Building octagonal column section for column ColSecTag")

       ftfact = 1.0;
       fcfact = 1.0;

       global 'kips', in 'ksi_psi', pi 'fce', ec0 esp 'Ec', Gc 'Es', Esh 'fy', fu 'esh', esu; # Read 'in', unit 'and', constant 'variables', as 'well', as 'material', property 'variables'

       IDconcCore =  (ColSecTag*10+1)
       IDconcCover = (ColSecTag*10+2)
       IDSteel =     (ColSecTag*10+3)
       IDShear =     (ColSecTag*10+4)
       IDTorsion =   (ColSecTag*10+5)

       # Column component dimensions
       tcover = 2.0*inch                                                 # 2 inch cover width
       Rcol = Dcol/2.0;                                                  # Radius of octagonal column (to flat sides)
       RcolDiag = Dcol*sqrt(4.0+2.0*sqrt(2.0))/(2.0+2.0*sqrt(2.0));      # Radius of octagonal column (to corners)
       Dcore = Dcol-2.0*tcover;                                          # Diameter of circular core
       Rcore = Dcore/2.0;                                                        # Radius of circular core
       Along = pi*DLbar**2/4.0;                                          # Area of longitudinal reinforcement bar
       DTbar = 0.625*inch;                                                      # Diameter of transverse spiral reinforcement bar (#5 Rebar)
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

       # Material Definition - Steel (Elastic)
       model.uniaxialMaterial('Elastic', IDSteel Es);

       # Material Definition - Concrete (Concrete02)
       # Compressive Properties
       ke = 1-sTbar/Dtran;                                   # Effective confinement strength coefficient from transverse reinforcement
       f3e = ke*rho*fy/2;                                    # Effective confinement strength
       fcc = fce*(-1.254+2.254*sqrt(1+7.94*f3e/fce)-2*f3e/fce);       # Core compressive strength -compression
       ecc = ec0*(1.0+5.0*(fcc/fce-1.0));       # Core strain at maximum strength -compression
       ecu = 0.004+f3e/(4.0*fce);                     # Core crushing (ultimate) strain -compression
       lamda = 0.1                                                       # Ratio between unloading slope at eps2 and initial slope Ec
       xu = ecu/ecc;                                          # Mander equation parameter x
       ru = Ec/(Ec-fcc/ecc);                      # Mander equation parameter r
       fcu = fcfact*fcc*xu*ru/(ru-1+xu**ru);       # Core crushing (ultimate) strength -compression

       ##### PRINT LINE TO CHECK #####
       # print("f3e=f3e, fcc=fcc, ecc=ecc, ecu=ecu, fcu=fcu");
       # print("xu=xu, ru=ru, fcu=fcu");

       # Tensile Properties
       ftU = 7.5*sqrt(fce*ksi_psi)/ksi_psi;              # Cover tensile strength +tension
       ftC = ftfact*7.5*sqrt(fcc*ksi_psi)/ksi_psi;              # Core tensile strength +tension
       Ets = Ec/5.0                                                       # Tension softening stiffness *Check if causes numerical issues (Divide by constant to make flatter if needed)

       ##### PRINT LINE TO CHECK #####
       # print("ftU=ftU, ftC=ftC, Ets=Ets")

       # UniaxialMaterial Definition
       # uniaxialMaterial Concrete02 IDconcCover, (-fce) [-expr ec0], (-0.1*fce), (-esp) lamda ftU Ets;       # Cover concrete (unconfined)
       model.uniaxialMaterial('Elastic', IDconcCover Ec);
       ##### PRINT LINE TO CHECK #####
       # print("Cover fpc = (-fce), epsc0 =, (-Efact*ec0)")
       # uniaxialMaterial Concrete02 IDconcCore, (-fcc), (-ecc), (-fcu), (-ecu) lamda ftC Ets;       # Core concrete (confined)
       model.uniaxialMaterial('Elastic', IDconcCore Ec);
       ##### PRINT LINE TO CHECK #####
       # print("Core fpc = (-fcc), epsc0 =, (-Efact*ecc)")

       # Build Octagonal RC Column Section
       numSubdivCirc = 64; # Determines # of fibers in the circumferential direction, around the entire circumference
       ColMatTag = 10*ColSecTag;    # Column material cross-section tag.  Set equal to 10*Column's fiber model.element(tag.)
       # Define elastic torsional stiffness: Based on SDC-1.6 sec. 5.6.2 and sec. 1.1,
       # fundamental period less than 0.7 sec, reduction is required
       # *****CHECK FUNDAMENTAL PERIOD TO SEE IF REDUCTION STILL REQUIRED
       model.uniaxialMaterial('Elastic', IDTorsion '[expr',  0.2*Gc*Jcol]);
       section.Fiber(ColMatTag torsion=IDTorsion, , [
              patch.circ(IDconcCore numSubdivCirc 5 0. 0., ),(0.000e+00*in) '[expr', Rcore/2] 0.000e+00 3.600e+02; # Inner 'Core', Patch, 5 radial 'fibers'
              patch.circ(IDconcCore numSubdivCirc 10 0. 0., ),(Rcore/2) '[expr', Rcore] 0.000e+00 3.600e+02; # Outer 'Core', Patch, 10 radial 'fibers'
              numSlices = 8;       # Determines # of slices in each of the 8 sections of the octagon
              numSubdivIJ = 1;       # Determines # of fibers in the circumferential direction of the cover patch
              numSubdivJK = 4;       # Determines # of fibers in the radial direction of the cover patch
              for {i = 0} {i < 8} {i += 1} { # For 'each', of 'the', 8 sections 'of', the 'octagon'
                     startAngle =       (i*pi/4+pi/8)
                     for {j = 0} {j < numSlices} {j += 1} {
                            phi =              (pi/4/numSlices)
                            sita1 = startAngle + j*phi;       # Slice start angle
                            sita2 = sita1 + phi;                     # Slice end angle
                            yI =              (Rcore*cos(sita2))
                            zI =              (Rcore*sin(sita2))
                            yJ =              (Rcore*cos(sita1))
                            zJ =              (Rcore*sin(sita1))
                            oR1 =              (Rcol/cos(pi/8 - j*phi))
                            oR2 =              (Rcol/cos(pi/8 - (j+1)*phi))
                            yK =              (oR1*cos(sita1))
                            zK =              (oR1*sin(sita1))
                            yL =              (oR2*cos(sita2))
                            zL =              (oR2*sin(sita2))
                            patch.quad(IDconcCover numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal 'cover'


              layer.circ(IDSteel nLbar Along 0. 0. Rlong), # Longitudinal 'Bars'

       model.uniaxialMaterial('Elastic', IDShear , (9./10.)*Gc*Acol       )# Define 'elastic', shear 'stiffness'
       section 'Aggregator', ColSecTag IDShear 'Vy', IDShear 'Vz', section=ColMatTag, ;


def BuildWideOctColSectionElastic(ColSecTag  Dcol  Wcol  nLbar  nLbar2  DLbar  DLbar2  sTbar):
  # OCTAGONAL WIDE COLUMN SECTION DEFINITION PROCEDURE - INELASTIC PROPERTIES WITH MODIFIED CRUSHING STRESS/STRAIN PAIR

  ################### Octagonal WIDE (Interlocking) Column Fiber Cross Section Assignment ###################

  # ColSecTag                           -> Column's fiber element tag
  # Dcol                                  -> Shorter Width of octagonal column (to flat sides)
  # Wcol                                  -> Longer Width of octagonal column (to flat sides)
  # nLbar                                  -> Number of main (outer) longitudinal bars
  # nLbar2                           -> Number of secondary (inner) longitudinal bars
  # DLbar                                   -> Diameter of main (outer) longitudinal bars
  # DLbar2                             -> Diameter of secondary (inner) longitudinal bars
  # sTbar                             -> Spacing of transverse spiral reinforcement
  # Local Y axis is the horizontal axis.  Local Z axis is the vertical axis.  The longer width of the column is in the local Y direction.

       ##### PRINT LINE TO CHECK #####
       # print("Building WIDE octagonal column section for column ColSecTag")

       ftfact = 1.0;
       fcfact = 1.0;

       global 'kips', in 'ksi_psi', pi 'fce', ec0 esp 'Ec', Gc 'Es', Esh 'fy', fu 'esh', esu; # Read 'in', unit 'and', constant 'variables', as 'well', as 'material', property 'variables'

       IDconcCore =  (ColSecTag*10+1)
       IDconcCover = (ColSecTag*10+2)
       IDSteel =     (ColSecTag*10+3)
       IDShear =     (ColSecTag*10+4)
       IDTorsion =   (ColSecTag*10+5)

       # Column component dimensions
       tcover = 2.0*inch                                                     # 2 inch cover width
       Rcol = Dcol/2.0;                                                        # Radius of octagonal column (to flat sides)
       RcolDiag = Dcol*sqrt(4.0+2.0*sqrt(2.0))/(2.0+2.0*sqrt(2.0));         # Radius of octagonal column (to corners)
       Dcore = Dcol-2.0*tcover;                                          # Diameter of circular core
       Dcore2 = Wcol-2.0*tcover;                                          # Longer width of core
       Rcore = Dcore/2.0;                                                        # Radius of circular core
       Along = pi*DLbar**2/4.0;                                          # Area of OUTER longitudinal reinforcement bar
       Along2 = pi*DLbar2**2/4.0;                                          # Area of INNER longitudinal reinforcement bar
       DTbar = 0.625*inch;                                                      # Diameter of transverse spiral reinforcement bar (#5 Rebar)
       Asp = pi*DTbar**2/4.0;                                          # Area of transverse spiral reinforcement bar
       Dtran = Dcore-DTbar;                                                 # Diameter of spiral of transverse spiral reinforcement
       Dtran2 = Dcore2-DTbar;                                                  # Longer width of transverse spiral reinforcement
       rho = 4.0*Asp/(Dtran*sTbar);                            # Density of transverse spiral reinforcement in shorter direction
       rho2 = 4.0*Asp/(Dtran2*sTbar);                            # Density of transverse spiral reinforcement in longer direction
       Dlong = Dcore-2*DTbar-DLbar;                                   # Diameter of ring of longitudinal reinforcement
       Rlong = Dlong/2.0;                                                        # Radius of ring of longitudinal reinforcement

       # Section Area Properties
       spO = (Wcol-Dcol)/2.0;                                                         # Offof = octagonal sections from centroid (along horizontal axis)
       Acol = (2*Dcol**2)/(1+sqrt(2)) + (Wcol-Dcol)*Dcol; # Area of WIDE octagonal column section
       I2col = 0.6381*RcolDiag**4;                                       # Second moment of inertia, horizontal axis, for WIDE octagonal column section
       I3col = I2col+Acol*spO**2+(Wcol-Dcol)*Dcol**3/12; # Second moment of inertia, vertical axis, for WIDE octagonal column section
       Jcol = I2col+I3col;                                              # Polar moment of inertia for WIDE octagonal column section

       # Material Definition - Steel (Elastic)
       R0 = 18                                                       # control the transition from elastic to plastic branches
       cR1 = 0.925;                                                        # control the transition from elastic to plastic branches
       cR2 = 0.15;                                                        # control the transition from elastic to plastic branches
       # uniaxialMaterial Steel02 IDSteel fy Es 0.02 R0 cR1 cR2
       model.uniaxialMaterial('Elastic', IDSteel Es);

       # Material Definition - Concrete (Concrete02)
       # Compressive Properties
       ke = 1-sTbar/Dtran;                                   # Effective confinement strength coefficient from transverse reinforcement, shorter direction
       ke2 = 1-sTbar/Dtran2;                            # Effective confinement strength coefficient from transverse reinforcement, longer direction
       f2e = ke*rho*fy/2;                                    # Effective confinement strength in the shorter direction
       f3e = ke2*rho2*fy/2;                             # Effective confinement strength in the longer direction
       fcc2 = fce*(-1.254+2.254*sqrt(1+7.94*f2e/fce)-2*f2e/fce);       # Reference core compressive strength for f2e f2e -compression
       fcc3 = fce*(-1.254+2.254*sqrt(1+7.94*f3e/fce)-2*f3e/fce);       # Reference core compressive strength for f3e f3e -compression
       fcc = 9.05                                                       # Actual core compressive strength for f2e f3e

       ##### PRINT LINE TO CHECK #####
       # print("fce=fce, f2e=f2e, f3e=f3e, f2e/fce=[expr f2e/fce], f3e/fce=[expr f3e/fce], fcc2=fcc2, fcc3=fcc3, fcc=fcc")

       fcc = fce*(-1.254+2.254*sqrt(1+7.94*f3e/fce)-2*f3e/fce);       # Core compressive strength -compression
       ecc = ec0*(1.0+5.0*(fcc/fce-1.0));       # Core strain at maximum strength -compression
       ecu = 0.004+f3e/(4.0*fce);                     # Core crushing (ultimate) strain -compression
       lamda = 0.1                                                       # Ratio between unloading slope at eps2 and initial slope Ec
       xu = ecu/ecc;                                          # Mander equation parameter x
       ru = Ec/(Ec-fcc/ecc);                      # Mander equation parameter r
       fcu = fcfact*fcc*xu*ru/(ru-1+xu**ru);       # Core crushing (ultimate) strength -compression
       fcuu = 0.15*fcc;                                           # Modified Core crushing (ultimate) strength, equal to 0.15x core compressive strength -compression
       ecuu = ecu+(ecu-ecc)*(fcuu-fcu)/(fcu-fcc); # Modified Core crushing (ultimate) strain, corresponding to fcuu, calculated by extrapolating the line between (ecc,fcc) and (ecu,fcu) -compression

       ##### PRINT LINE TO CHECK #####
       # print("f3e=f3e, fcc=fcc, ecc=ecc, ecu=ecu, fcu=fcu, ecuu=ecuu, fcuu=fcuu");
       # print("xu=xu, ru=ru, fcu=fcu");

       # Tensile Properties
       ftU = 0;       # Cover tensile strength +tension
       ftC = 0;       # Core tensile strength +tension
       Ets = 0;       # Tension softening stiffness

       ##### PRINT LINE TO CHECK #####
       # print("ftU=ftU, ftC=ftC, Ets=Ets")

       # UniaxialMaterial Definition
       # uniaxialMaterial Concrete02 IDconcCover, (-fce), (-ec0), (-0.1*fce), (-esp) lamda ftU Ets;       # Cover concrete (unconfined)
       # print("uniaxialMaterial Concrete02 IDconcCover, (-fce), (-ec0), (-0.1*fce), (-esp) lamda ftU Ets;       # Cover concrete (unconfined)")
       model.uniaxialMaterial('Elastic', IDconcCover Ec);
       ##### PRINT LINE TO CHECK #####
       # print("Cover fpc = (-fce), epsc0 =, (-Efact*ec0)")
       # uniaxialMaterial Concrete02 IDconcCore, (-fcc), (-ecc), (-fcuu), (-ecuu) lamda ftC Ets;       # Core concrete (confined)
       # print("uniaxialMaterial Concrete02 IDconcCore, (-fcc), (-ecc), (-fcuu), (-ecuu) lamda ftC Ets;       # Core concrete (confined)")
       model.uniaxialMaterial('Elastic', IDconcCore Ec);
       ##### PRINT LINE TO CHECK #####
       # print("Core fpc = (-fcc), epsc0 =, (-Efact*ecc)")

       # Build Octagonal RC Column Section
       numSubdivCirc = 64; # Determines # of fibers in the circumferential direction, around the entire circumference
       ColMatTag = 10*ColSecTag;    # Column material cross-section tag.  Set equal to 10*Column's fiber model.element(tag.)
       # Define elastic torsional stiffness: Based on SDC-1.6 sec. 5.6.2 and sec. 1.1,
       # fundamental period less than 0.7 sec, reduction is required
       # *****CHECK FUNDAMENTAL PERIOD TO SEE IF REDUCTION STILL REQUIRED
       model.uniaxialMaterial('Elastic', IDTorsion '[expr',  0.2*Gc*Jcol]);
       section.Fiber(ColMatTag torsion=IDTorsion, , [
              patch.circ(IDconcCore numSubdivCirc 5, ),(-spO)  0., (0.000e+00*in) '[expr', Rcore/2] 90.0 270.0; # LEFT 'Inner', circular 'Core', Patch, 5 radial 'fibers'
              patch.circ(IDconcCore numSubdivCirc 5, ),(spO) 0., (0.000e+00*in) '[expr', Rcore/2] 270.0 450.0; # RIGHT 'Inner', circular 'Core', Patch, 5 radial 'fibers'
              patch.circ(IDconcCore numSubdivCirc 10, ),(-spO) 0., (Rcore/2) '[expr', Rcore] 90.0 270.0; # LEFT 'Outer', circular 'Core', Patch, 10 radial 'fibers'
              patch.circ(IDconcCore numSubdivCirc 10, ),(spO) 0., (Rcore/2) '[expr', Rcore] 270.0 450.0; # RIGHT 'Outer', circular 'Core', Patch, 10 radial 'fibers'
              patch.quad(IDconcCore 5 5, ),(-spO) '[expr', -Rcore/2] '[expr', spO] '[expr', -Rcore/2] '[expr', spO] '[expr', Rcore/2] '[expr', -spO] '[expr', Rcore/2]; # Inner 'rectangular', Core Patch, 5 fibers 'in', y & z
              patch.quad(IDconcCore 10 10, ),(-spO) '[expr', -Rcore] '[expr', spO] '[expr', -Rcore] '[expr', spO] '[expr', -Rcore/2] '[expr', -spO] '[expr', -Rcore/2]; # BOTTOM 'Outer', rectangular 'Core', Patch, 10 fibers 'in', y & z
              patch.quad(IDconcCore 10 10, ),(-spO) '[expr', Rcore/2] '[expr', spO] '[expr', Rcore/2] '[expr', spO] '[expr', Rcore] '[expr', -spO] '[expr', Rcore]; # TOP 'Outer', rectangular 'Core', Patch, 10 fibers 'in', y & z
              patch.quad(IDconcCover 10 4, ),(-spO) '[expr', -Rcol] '[expr', spO] '[expr', -Rcol] '[expr', spO] '[expr', -Rcore] '[expr', -spO] '[expr', -Rcore]; # BOTTOM 'Outer', rectangular 'COVER', Patch, 10 fibers 'in', y, 2 fibers 'in', z
              patch.quad(IDconcCover 10 4, ),(-spO) '[expr', Rcore] '[expr', spO] '[expr', Rcore] '[expr', spO] '[expr', Rcol] '[expr', -spO] '[expr', Rcol]; # TOP 'Outer', rectangular 'COVER', Patch, 10 fibers 'in', y, 2 fibers 'in', z
              numSlices = 8;       # Determines # of slices in each of the 8 sections of the octagon
              numSubdivIJ = 1;       # Determines # of fibers in the circumferential direction of the cover patch
              numSubdivJK = 4;       # Determines # of fibers in the radial direction of the cover patch
              for {i = 0} {i < 2} {i += 1} { # For 'the', first 2 sections 'of', the 'octagon', (TOP RIGHT; offset 'all', y 'coordinates', by +spO)
                     startAngle =       (i*pi/4+pi/8)
                     for {j = 0} {j < numSlices} {j += 1} {
                            phi =              (pi/4/numSlices)
                            sita1 = startAngle + j*phi;       # Slice start angle
                            sita2 = sita1 + phi;                     # Slice end angle
                            yI =              (Rcore*cos(sita2)+spO)
                            zI =              (Rcore*sin(sita2))
                            yJ =              (Rcore*cos(sita1)+spO)
                            zJ =              (Rcore*sin(sita1))
                            oR1 =              (Rcol/cos(pi/8 - j*phi))
                            oR2 =              (Rcol/cos(pi/8 - (j+1)*phi))
                            yK =              (oR1*cos(sita1)+spO)
                            zK =              (oR1*sin(sita1))
                            yL =              (oR2*cos(sita2)+spO)
                            zL =              (oR2*sin(sita2))
                            patch.quad(IDconcCover numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal 'cover'


              for {i = 2} {i < 6} {i += 1} { # For 'the', 3rd 'through', 6th 'sections', of 'the', octagon (TOP & BOTTOM LEFT; offset 'all', y 'coordinates', by -spO)
                     startAngle =       (i*pi/4+pi/8)
                     for {j = 0} {j < numSlices} {j += 1} {
                            phi =              (pi/4/numSlices)
                            sita1 = startAngle + j*phi;       # Slice start angle
                            sita2 = sita1 + phi;                     # Slice end angle
                            yI =              (Rcore*cos(sita2)-spO)
                            zI =              (Rcore*sin(sita2))
                            yJ =              (Rcore*cos(sita1)-spO)
                            zJ =              (Rcore*sin(sita1))
                            oR1 =              (Rcol/cos(pi/8 - j*phi))
                            oR2 =              (Rcol/cos(pi/8 - (j+1)*phi))
                            yK =              (oR1*cos(sita1)-spO)
                            zK =              (oR1*sin(sita1))
                            yL =              (oR2*cos(sita2)-spO)
                            zL =              (oR2*sin(sita2))
                            patch.quad(IDconcCover numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal 'cover'


              for {i = 6} {i < 8} {i += 1} { # For 'the', last 2 sections 'of', the 'octagon', (BOTTOM RIGHT; offset 'all', y 'coordinates', by +spO)
                     startAngle =       (i*pi/4+pi/8)
                     for {j = 0} {j < numSlices} {j += 1} {
                            phi =              (pi/4/numSlices)
                            sita1 = startAngle + j*phi;       # Slice start angle
                            sita2 = sita1 + phi;                     # Slice end angle
                            yI =              (Rcore*cos(sita2)+spO)
                            zI =              (Rcore*sin(sita2))
                            yJ =              (Rcore*cos(sita1)+spO)
                            zJ =              (Rcore*sin(sita1))
                            oR1 =              (Rcol/cos(pi/8 - j*phi))
                            oR2 =              (Rcol/cos(pi/8 - (j+1)*phi))
                            yK =              (oR1*cos(sita1)+spO)
                            zK =              (oR1*sin(sita1))
                            yL =              (oR2*cos(sita2)+spO)
                            zL =              (oR2*sin(sita2))
                            patch.quad(IDconcCover numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal 'cover'


              layer.circ(IDSteel '[expr', nLbar),/2] Along '[expr', -spO] 0. Rlong 45.0 315.0; # LEFT 'outer', Longitudinal 'Bars'
              layer.circ(IDSteel '[expr', nLbar),/2] Along '[expr', spO] 0. Rlong 225.0 495.0; # RIGHT 'outer', Longitudinal 'Bars'
              layer.circ(IDSteel '[expr', nLbar2),/2] Along2, (-spO) 0. Rlong 315.0 405.0; # RIGHT 'inner', Longitudinal 'Bars'
              layer.circ(IDSteel '[expr', nLbar2),/2] Along2, (spO) 0. Rlong 135.0 225.0; # LEFT 'inner', Longitudinal 'Bars'

       model.uniaxialMaterial('Elastic', IDShear , (9./10.)*Gc*Acol       )# Define 'elastic', shear 'stiffness'
       section 'Aggregator', ColSecTag IDShear 'Vy', IDShear 'Vz', section=ColMatTag, ;


#--------------------------------------------------------------------------
# From ColSectionOctIE
#--------------------------------------------------------------------------
def BuildOctColSectionInel(ColSecTag  Dcol  nLbar  DLbar  sTbar):
    # OCTAGONAL COLUMN SECTION DEFINITION PROCEDURE - INELASTIC PROPERTIES

    ################### Octagonal Column Fiber Cross Section Assignment ###################

    # ColSecTag                           -> Column's fiber element tag
    # Dcol                                  -> Width of octagonal column (to flat sides)
    # nLbar                                  -> Number of longitudinal bars
    # DLbar                                   -> Diameter of longitudinal bars
    # sTbar                             -> Spacing of transverse spiral reinforcement


       ##### PRINT LINE TO CHECK #####
       # print("Building octagonal column section for column ColSecTag")

       ftfact = 1.0;
       fcfact = 1.0;

       global 'kips', in 'ksi_psi', pi 'fce', ec0 esp 'Ec', Gc 'Es', Esh 'fy', fu 'esh', esu; # Read 'in', unit 'and', constant 'variables', as 'well', as 'material', property 'variables'

       IDconcCore =  (ColSecTag*10+1)
       IDconcCover = (ColSecTag*10+2)
       IDSteel =     (ColSecTag*10+3)
       IDShear =     (ColSecTag*10+4)
       IDTorsion =   (ColSecTag*10+5)

       # Column component dimensions
       tcover = 2.0*inch                                                     # 2 inch cover width
       Rcol = Dcol/2.0;                                                        # Radius of octagonal column (to flat sides)
       RcolDiag = Dcol*sqrt(4.0+2.0*sqrt(2.0))/(2.0+2.0*sqrt(2.0));         # Radius of octagonal column (to corners)
       Dcore = Dcol-2.0*tcover;                                          # Diameter of circular core
       Rcore = Dcore/2.0;                                                        # Radius of circular core
       Along = pi*DLbar**2/4.0;                                          # Area of longitudinal reinforcement bar
       DTbar = 0.625*inch;                                                      # Diameter of transverse spiral reinforcement bar (#5 Rebar)
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
       R0 = 18                                                       # control the transition from elastic to plastic branches
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
       lamda = 0.1                                                       # Ratio between unloading slope at eps2 and initial slope Ec
       xu = ecu/ecc;                                          # Mander equation parameter x
       ru = Ec/(Ec-fcc/ecc);                      # Mander equation parameter r
       fcu = fcfact*fcc*xu*ru/(ru-1+xu**ru);       # Core crushing (ultimate) strength -compression

       ##### PRINT LINE TO CHECK #####
       # print("f3e=f3e, fcc=fcc, ecc=ecc, ecu=ecu, fcu=fcu");
       # print("xu=xu, ru=ru, fcu=fcu");

       # Tensile Properties
       ftU = 0;       # Cover tensile strength +tension
       ftC = 0;       # Core tensile strength +tension
       Ets = 0;       # Tension softening stiffness

       ##### PRINT LINE TO CHECK #####
       # print("ftU=ftU, ftC=ftC, Ets=Ets")

       # UniaxialMaterial Definition
       model.uniaxialMaterial('Concrete02', IDconcCover '[expr', -fce] '[expr', -ec0] '[expr', -0.1*fce] '[expr', -esp] lamda ftU Ets);       # Cover 'concrete', (unconfined)
       # print("uniaxialMaterial Concrete02 IDconcCover, (-fce), (-ec0), (-0.1*fce), (-esp) lamda ftU Ets;       # Cover concrete (unconfined)")

       ##### PRINT LINE TO CHECK #####
       # print("Cover fpc = (-fce), epsc0 =, (-Efact*ec0)")

       model.uniaxialMaterial('Concrete02', IDconcCore, (-fcc) '[expr', -ecc] '[expr', -fcu] '[expr', -ecu] lamda ftC Ets);       # Core 'concrete', (confined)
       # print("uniaxialMaterial Concrete02 IDconcCore, (-fcc), (-ecc), (-fcu), (-ecu) lamda ftC Ets;       # Core concrete (confined)")

       # Build Octagonal RC Column Section
       numSubdivCirc = 64; # Determines # of fibers in the circumferential direction, around the entire circumference
       ColMatTag = 10*ColSecTag;    # Column material cross-section tag.  Set equal to 10*Column's fiber model.element(tag.)
       # Define elastic torsional stiffness: Based on SDC-1.6 sec. 5.6.2 and sec. 1.1,
       # fundamental period less than 0.7 sec, reduction is required
       # *****CHECK FUNDAMENTAL PERIOD TO SEE IF REDUCTION STILL REQUIRED
       model.uniaxialMaterial('Elastic', IDTorsion '[expr',  0.2*Gc*Jcol]);
       section.Fiber(ColMatTag torsion=IDTorsion, , [
              patch.circ(IDconcCore numSubdivCirc 5 0. 0., ),(0.000e+00*in) '[expr', Rcore/2] 0.000e+00 3.600e+02; # Inner 'Core', Patch, 5 radial 'fibers'
              patch.circ(IDconcCore numSubdivCirc 10 0. 0., ),(Rcore/2) '[expr', Rcore] 0.000e+00 3.600e+02; # Outer 'Core', Patch, 10 radial 'fibers'
              numSlices = 8;       # Determines # of slices in each of the 8 sections of the octagon
              numSubdivIJ = 1;       # Determines # of fibers in the circumferential direction of the cover patch
              numSubdivJK = 4;       # Determines # of fibers in the radial direction of the cover patch
              for {i = 0} {i < 8} {i += 1} { # For 'each', of 'the', 8 sections 'of', the 'octagon'
                     startAngle =       (i*pi/4+pi/8)
                     for {j = 0} {j < numSlices} {j += 1} {
                            phi =              (pi/4/numSlices)
                            sita1 = startAngle + j*phi;       # Slice start angle
                            sita2 = sita1 + phi;                     # Slice end angle
                            yI =              (Rcore*cos(sita2))
                            zI =              (Rcore*sin(sita2))
                            yJ =              (Rcore*cos(sita1))
                            zJ =              (Rcore*sin(sita1))
                            oR1 =              (Rcol/cos(pi/8 - j*phi))
                            oR2 =              (Rcol/cos(pi/8 - (j+1)*phi))
                            yK =              (oR1*cos(sita1))
                            zK =              (oR1*sin(sita1))
                            yL =              (oR2*cos(sita2))
                            zL =              (oR2*sin(sita2))
                            patch.quad(IDconcCover numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal 'cover'


              layer.circ(IDSteel nLbar Along 0. 0. Rlong), # Longitudinal 'Bars'

       model.uniaxialMaterial('Elastic', IDShear , (9./10.)*Gc*Acol       )# Define 'elastic', shear 'stiffness'
       section 'Aggregator', ColSecTag IDShear 'Vy', IDShear 'Vz', section=ColMatTag, ;

#--------------------------------------------------------------------------
# From ColSectionOctWideIE.tcl
#--------------------------------------------------------------------------
def BuildWideOctColSectionInel(ColSecTag  Dcol  Wcol  nLbar  nLbar2  DLbar  DLbar2  sTbar):
  # OCTAGONAL WIDE COLUMN SECTION DEFINITION PROCEDURE - INELASTIC PROPERTIES WITH MODIFIED CRUSHING STRESS/STRAIN PAIR

  ################### Octagonal WIDE (Interlocking) Column Fiber Cross Section Assignment ###################

  # ColSecTag                           -> Column's fiber element tag
  # Dcol                                  -> Shorter Width of octagonal column (to flat sides)
  # Wcol                                  -> Longer Width of octagonal column (to flat sides)
  # nLbar                                  -> Number of main (outer) longitudinal bars
  # nLbar2                           -> Number of secondary (inner) longitudinal bars
  # DLbar                                   -> Diameter of main (outer) longitudinal bars
  # DLbar2                             -> Diameter of secondary (inner) longitudinal bars
  # sTbar                             -> Spacing of transverse spiral reinforcement
  # Local Y axis is the horizontal axis.  Local Z axis is the vertical axis.  The longer width of the column is in the local Y direction.
       ##### PRINT LINE TO CHECK #####
       # print("Building WIDE octagonal column section for column ColSecTag")

       ftfact = 1.0;
       fcfact = 1.0;

       global 'kips', in 'ksi_psi', pi 'fce', ec0 esp 'Ec', Gc 'Es', Esh 'fy', fu 'esh', esu; # Read 'in', unit 'and', constant 'variables', as 'well', as 'material', property 'variables'

       IDconcCore =  (ColSecTag*10+1)
       IDconcCover = (ColSecTag*10+2)
       IDSteel =     (ColSecTag*10+3)
       IDShear =     (ColSecTag*10+4)
       IDTorsion =   (ColSecTag*10+5)

       # Column component dimensions
       tcover = 2.0*inch                                                     # 2 inch cover width
       Rcol = Dcol/2.0;                                                        # Radius of octagonal column (to flat sides)
       RcolDiag = Dcol*sqrt(4.0+2.0*sqrt(2.0))/(2.0+2.0*sqrt(2.0));         # Radius of octagonal column (to corners)
       Dcore = Dcol-2.0*tcover;                                          # Diameter of circular core
       Dcore2 = Wcol-2.0*tcover;                                          # Longer width of core
       Rcore = Dcore/2.0;                                                        # Radius of circular core
       Along = pi*DLbar**2/4.0;                                          # Area of OUTER longitudinal reinforcement bar
       Along2 = pi*DLbar2**2/4.0;                                          # Area of INNER longitudinal reinforcement bar
       DTbar = 0.625*inch;                                                      # Diameter of transverse spiral reinforcement bar (#5 Rebar)
       Asp = pi*DTbar**2/4.0;                                          # Area of transverse spiral reinforcement bar
       Dtran = Dcore-DTbar;                                                 # Diameter of spiral of transverse spiral reinforcement
       Dtran2 = Dcore2-DTbar;                                                  # Longer width of transverse spiral reinforcement
       rho = 4.0*Asp/(Dtran*sTbar);                            # Density of transverse spiral reinforcement in shorter direction
       rho2 = 4.0*Asp/(Dtran2*sTbar);                            # Density of transverse spiral reinforcement in longer direction
       Dlong = Dcore-2*DTbar-DLbar;                                   # Diameter of ring of longitudinal reinforcement
       Rlong = Dlong/2.0;                                                        # Radius of ring of longitudinal reinforcement

       # Section Area Properties
       spO = (Wcol-Dcol)/2.0;                                                         # Offof = octagonal sections from centroid (along horizontal axis)
       Acol = (2*Dcol**2)/(1+sqrt(2)) + (Wcol-Dcol)*Dcol; # Area of WIDE octagonal column section
       I2col = 0.6381*RcolDiag**4;                                       # Second moment of inertia, horizontal axis, for WIDE octagonal column section
       I3col = I2col+Acol*spO**2+(Wcol-Dcol)*Dcol**3/12; # Second moment of inertia, vertical axis, for WIDE octagonal column section
       Jcol = I2col+I3col;                                              # Polar moment of inertia for WIDE octagonal column section

       # Material Definition - Steel (Steel02)
       R0 = 18                                                       # control the transition from elastic to plastic branches
       cR1 = 0.925;                                                        # control the transition from elastic to plastic branches
       cR2 = 0.15;                                                        # control the transition from elastic to plastic branches
       model.uniaxialMaterial('Steel02', IDSteel fy Es 0.02 R0 cR1 cR2)

       # Material Definition - Concrete (Concrete02)
       # Compressive Properties
       ke = 1-sTbar/Dtran;                                   # Effective confinement strength coefficient from transverse reinforcement, shorter direction
       ke2 = 1-sTbar/Dtran2;                            # Effective confinement strength coefficient from transverse reinforcement, longer direction
       f2e = ke*rho*fy/2;                                    # Effective confinement strength in the shorter direction
       f3e = ke2*rho2*fy/2;                             # Effective confinement strength in the longer direction
       fcc2 = fce*(-1.254+2.254*sqrt(1+7.94*f2e/fce)-2*f2e/fce);       # Reference core compressive strength for f2e f2e -compression
       fcc3 = fce*(-1.254+2.254*sqrt(1+7.94*f3e/fce)-2*f3e/fce);       # Reference core compressive strength for f3e f3e -compression
       fcc = 9.05                                                       # Actual core compressive strength for f2e f3e

       ##### PRINT LINE TO CHECK #####
       # print("fce=fce, f2e=f2e, f3e=f3e, f2e/fce=[expr f2e/fce], f3e/fce=[expr f3e/fce], fcc2=fcc2, fcc3=fcc3, fcc=fcc")

       fcc = fce*(-1.254+2.254*sqrt(1+7.94*f3e/fce)-2*f3e/fce);       # Core compressive strength -compression
       ecc = ec0*(1.0+5.0*(fcc/fce-1.0));       # Core strain at maximum strength -compression
       ecu = 0.004+f3e/(4.0*fce);                     # Core crushing (ultimate) strain -compression
       lamda = 0.1                                                       # Ratio between unloading slope at eps2 and initial slope Ec
       xu = ecu/ecc;                                          # Mander equation parameter x
       ru = Ec/(Ec-fcc/ecc);                      # Mander equation parameter r
       fcu = fcfact*fcc*xu*ru/(ru-1+xu**ru);       # Core crushing (ultimate) strength -compression
       fcuu = 0.15*fcc;                                           # Modified Core crushing (ultimate) strength, equal to 0.15x core compressive strength -compression
       ecuu = ecu+(ecu-ecc)*(fcuu-fcu)/(fcu-fcc); # Modified Core crushing (ultimate) strain, corresponding to fcuu, calculated by extrapolating the line between (ecc,fcc) and (ecu,fcu) -compression

       ##### PRINT LINE TO CHECK #####
       # print("f3e=f3e, fcc=fcc, ecc=ecc, ecu=ecu, fcu=fcu, ecuu=ecuu, fcuu=fcuu");
       # print("xu=xu, ru=ru, fcu=fcu");

       # Tensile Properties
       ftU = 0;       # Cover tensile strength +tension
       ftC = 0;       # Core tensile strength +tension
       Ets = 0;       # Tension softening stiffness

       ##### PRINT LINE TO CHECK #####
       # print("ftU=ftU, ftC=ftC, Ets=Ets")

       # UniaxialMaterial Definition
       model.uniaxialMaterial('Concrete02', IDconcCover '[expr', -fce] '[expr', -ec0] '[expr', -0.1*fce] '[expr', -esp] lamda ftU Ets);       # Cover 'concrete', (unconfined)
       # print("uniaxialMaterial Concrete02 IDconcCover, (-fce), (-ec0), (-0.1*fce), (-esp) lamda ftU Ets;       # Cover concrete (unconfined)")
       ##### PRINT LINE TO CHECK #####
       # print("Cover fpc = (-fce), epsc0 =, (-Efact*ec0)")

       model.uniaxialMaterial('Concrete02', IDconcCore, (-fcc) '[expr', -ecc] '[expr', -fcu] '[expr', -ecu] lamda ftC Ets);       # Core 'concrete', (confined)
       # print("uniaxialMaterial Concrete02 IDconcCore, (-fcc), (-ecc), (-fcu), (-ecu) lamda ftC Ets;       # Core concrete (confined)")
       ##### PRINT LINE TO CHECK #####
       # print("Core fpc = (-fcc), epsc0 =, (-Efact*ecc)")

       # Build Octagonal RC Column Section
       numSubdivCirc = 64; # Determines # of fibers in the circumferential direction, around the entire circumference
       ColMatTag = 10*ColSecTag;    # Column material cross-section tag.  Set equal to 10*Column's fiber model.element(tag.)
       # Define elastic torsional stiffness: Based on SDC-1.6 sec. 5.6.2 and sec. 1.1,
       # fundamental period less than 0.7 sec, reduction is required
       # *****CHECK FUNDAMENTAL PERIOD TO SEE IF REDUCTION STILL REQUIRED
       model.uniaxialMaterial('Elastic', IDTorsion '[expr',  0.2*Gc*Jcol]);
       section.Fiber(ColMatTag torsion=IDTorsion, , [
              patch.circ(IDconcCore numSubdivCirc 5, ),(-spO)  0., (0.000e+00*in) '[expr', Rcore/2] 90.0 270.0; # LEFT 'Inner', circular 'Core', Patch, 5 radial 'fibers'
              patch.circ(IDconcCore numSubdivCirc 5, ),(spO) 0., (0.000e+00*in) '[expr', Rcore/2] 270.0 450.0; # RIGHT 'Inner', circular 'Core', Patch, 5 radial 'fibers'
              patch.circ(IDconcCore numSubdivCirc 10, ),(-spO) 0., (Rcore/2) '[expr', Rcore] 90.0 270.0; # LEFT 'Outer', circular 'Core', Patch, 10 radial 'fibers'
              patch.circ(IDconcCore numSubdivCirc 10, ),(spO) 0., (Rcore/2) '[expr', Rcore] 270.0 450.0; # RIGHT 'Outer', circular 'Core', Patch, 10 radial 'fibers'
              patch.quad(IDconcCore 5 5, ),(-spO) '[expr', -Rcore/2] '[expr', spO] '[expr', -Rcore/2] '[expr', spO] '[expr', Rcore/2] '[expr', -spO] '[expr', Rcore/2]; # Inner 'rectangular', Core Patch, 5 fibers 'in', y & z
              patch.quad(IDconcCore 10 10, ),(-spO) '[expr', -Rcore] '[expr', spO] '[expr', -Rcore] '[expr', spO] '[expr', -Rcore/2] '[expr', -spO] '[expr', -Rcore/2]; # BOTTOM 'Outer', rectangular 'Core', Patch, 10 fibers 'in', y & z
              patch.quad(IDconcCore 10 10, ),(-spO) '[expr', Rcore/2] '[expr', spO] '[expr', Rcore/2] '[expr', spO] '[expr', Rcore] '[expr', -spO] '[expr', Rcore]; # TOP 'Outer', rectangular 'Core', Patch, 10 fibers 'in', y & z
              patch.quad(IDconcCover 10 4, ),(-spO) '[expr', -Rcol] '[expr', spO] '[expr', -Rcol] '[expr', spO] '[expr', -Rcore] '[expr', -spO] '[expr', -Rcore]; # BOTTOM 'Outer', rectangular 'COVER', Patch, 10 fibers 'in', y, 2 fibers 'in', z
              patch.quad(IDconcCover 10 4, ),(-spO) '[expr', Rcore] '[expr', spO] '[expr', Rcore] '[expr', spO] '[expr', Rcol] '[expr', -spO] '[expr', Rcol]; # TOP 'Outer', rectangular 'COVER', Patch, 10 fibers 'in', y, 2 fibers 'in', z
              numSlices = 8;       # Determines # of slices in each of the 8 sections of the octagon
              numSubdivIJ = 1;       # Determines # of fibers in the circumferential direction of the cover patch
              numSubdivJK = 4;       # Determines # of fibers in the radial direction of the cover patch
              for {i = 0} {i < 2} {i += 1} { # For 'the', first 2 sections 'of', the 'octagon', (TOP RIGHT; offset 'all', y 'coordinates', by +spO)
                     startAngle =       (i*pi/4+pi/8)
                     for {j = 0} {j < numSlices} {j += 1} {
                            phi =              (pi/4/numSlices)
                            sita1 = startAngle + j*phi;       # Slice start angle
                            sita2 = sita1 + phi;                     # Slice end angle
                            yI =              (Rcore*cos(sita2)+spO)
                            zI =              (Rcore*sin(sita2))
                            yJ =              (Rcore*cos(sita1)+spO)
                            zJ =              (Rcore*sin(sita1))
                            oR1 =              (Rcol/cos(pi/8 - j*phi))
                            oR2 =              (Rcol/cos(pi/8 - (j+1)*phi))
                            yK =              (oR1*cos(sita1)+spO)
                            zK =              (oR1*sin(sita1))
                            yL =              (oR2*cos(sita2)+spO)
                            zL =              (oR2*sin(sita2))
                            patch.quad(IDconcCover numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal 'cover'


              for {i = 2} {i < 6} {i += 1} { # For 'the', 3rd 'through', 6th 'sections', of 'the', octagon (TOP & BOTTOM LEFT; offset 'all', y 'coordinates', by -spO)
                     startAngle =       (i*pi/4+pi/8)
                     for {j = 0} {j < numSlices} {j += 1} {
                            phi =              (pi/4/numSlices)
                            sita1 = startAngle + j*phi;       # Slice start angle
                            sita2 = sita1 + phi;                     # Slice end angle
                            yI =              (Rcore*cos(sita2)-spO)
                            zI =              (Rcore*sin(sita2))
                            yJ =              (Rcore*cos(sita1)-spO)
                            zJ =              (Rcore*sin(sita1))
                            oR1 =              (Rcol/cos(pi/8 - j*phi))
                            oR2 =              (Rcol/cos(pi/8 - (j+1)*phi))
                            yK =              (oR1*cos(sita1)-spO)
                            zK =              (oR1*sin(sita1))
                            yL =              (oR2*cos(sita2)-spO)
                            zL =              (oR2*sin(sita2))
                            patch.quad(IDconcCover numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal 'cover'


              for {i = 6} {i < 8} {i += 1} { # For 'the', last 2 sections 'of', the 'octagon', (BOTTOM RIGHT; offset 'all', y 'coordinates', by +spO)
                     startAngle =       (i*pi/4+pi/8)
                     for {j = 0} {j < numSlices} {j += 1} {
                            phi =              (pi/4/numSlices)
                            sita1 = startAngle + j*phi;       # Slice start angle
                            sita2 = sita1 + phi;                     # Slice end angle
                            yI =              (Rcore*cos(sita2)+spO)
                            zI =              (Rcore*sin(sita2))
                            yJ =              (Rcore*cos(sita1)+spO)
                            zJ =              (Rcore*sin(sita1))
                            oR1 =              (Rcol/cos(pi/8 - j*phi))
                            oR2 =              (Rcol/cos(pi/8 - (j+1)*phi))
                            yK =              (oR1*cos(sita1)+spO)
                            zK =              (oR1*sin(sita1))
                            yL =              (oR2*cos(sita2)+spO)
                            zL =              (oR2*sin(sita2))
                            patch.quad(IDconcCover numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal 'cover'


              layer.circ(IDSteel '[expr', nLbar),/2] Along '[expr', -spO] 0. Rlong 45.0 315.0; # LEFT 'outer', Longitudinal 'Bars'
              layer.circ(IDSteel '[expr', nLbar),/2] Along '[expr', spO] 0. Rlong 225.0 495.0; # RIGHT 'outer', Longitudinal 'Bars'
              layer.circ(IDSteel '[expr', nLbar2),/2] Along2, (-spO) 0. Rlong 315.0 405.0; # RIGHT 'inner', Longitudinal 'Bars'
              layer.circ(IDSteel '[expr', nLbar2),/2] Along2, (spO) 0. Rlong 135.0 225.0; # LEFT 'inner', Longitudinal 'Bars'

       model.uniaxialMaterial('Elastic', IDShear , (9./10.)*Gc*Acol       )# Define 'elastic', shear 'stiffness'
       section 'Aggregator', ColSecTag IDShear 'Vy', IDShear 'Vz', section=ColMatTag, ;


#---------------------------------------------------------------------
# From ColSectionOctIEu.tcl
#---------------------------------------------------------------------
def BuildOctColSectionInelMod(ColSecTag  Dcol  nLbar  DLbar  sTbar):
  # OCTAGONAL COLUMN SECTION DEFINITION PROCEDURE - INELASTIC PROPERTIES WITH MODIFIED CRUSHING STRESS/STRAIN PAIR

  ################### Octagonal Column Fiber Cross Section Assignment ###################

  # ColSecTag                           -> Column's fiber element tag
  # Dcol                                  -> Width of octagonal column (to flat sides)
  # nLbar                                  -> Number of longitudinal bars
  # DLbar                                   -> Diameter of longitudinal bars
  # sTbar                             -> Spacing of transverse spiral reinforcement
  # materials
  #        core
  #        cover
  #   long
  #   tran
  #   shear

       ##### PRINT LINE TO CHECK #####
       # print("Building octagonal column section for column ColSecTag")

       ftfact = 1.0;
       fcfact = 1.0;

       global 'kips', in 'ksi_psi', pi 'fce', ec0 esp 'Ec', Gc 'Es', Esh 'fy', fu 'esh', esu; # Read 'in', unit 'and', constant 'variables', as 'well', as 'material', property 'variables'

       IDconcCore =  (ColSecTag*10+1)
       IDconcCover = (ColSecTag*10+2)
       IDSteel =     (ColSecTag*10+3)
       IDShear =     (ColSecTag*10+4)
       IDTorsion =   (ColSecTag*10+5)

       # Column component dimensions
       tcover = 2.0*inch                                                     # 2 inch cover width
       Rcol = Dcol/2.0;                                                        # Radius of octagonal column (to flat sides)
       RcolDiag = Dcol*sqrt(4.0+2.0*sqrt(2.0))/(2.0+2.0*sqrt(2.0));         # Radius of octagonal column (to corners)
       Dcore = Dcol-2.0*tcover;                                          # Diameter of circular core
       Rcore = Dcore/2.0;                                                        # Radius of circular core
       Along = pi*DLbar**2/4.0;                                          # Area of longitudinal reinforcement bar
       DTbar = 0.625*inch;                                                      # Diameter of transverse spiral reinforcement bar (#5 Rebar)
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
       R0 = 18                                                       # control the transition from elastic to plastic branches
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
       lamda = 0.1                                                       # Ratio between unloading slope at eps2 and initial slope Ec
       xu = ecu/ecc;                                          # Mander equation parameter x
       ru = Ec/(Ec-fcc/ecc);                      # Mander equation parameter r
       fcu = fcfact*fcc*xu*ru/(ru-1+xu**ru);       # Core crushing (ultimate) strength -compression
       fcuu = 0.15*fcc;                                           # Modified Core crushing (ultimate) strength, equal to 0.15x core compressive strength -compression
       ecuu = ecu+(ecu-ecc)*(fcuu-fcu)/(fcu-fcc); # Modified Core crushing (ultimate) strain, corresponding to fcuu, calculated by extrapolating the line between (ecc,fcc) and (ecu,fcu) -compression

       ##### PRINT LINE TO CHECK #####
       # print("f3e=f3e, fcc=fcc, ecc=ecc, ecu=ecu, fcu=fcu, ecuu=ecuu, fcuu=fcuu");
       # print("xu=xu, ru=ru, fcu=fcu");

       # Tensile Properties
       ftU = 0;       # Cover tensile strength +tension
       ftC = 0;       # Core tensile strength +tension
       Ets = 0;       # Tension softening stiffness

       ##### PRINT LINE TO CHECK #####
       # print("ftU=ftU, ftC=ftC, Ets=Ets")

       # UniaxialMaterial Definition
       model.uniaxialMaterial('Concrete02', IDconcCover '[expr', -fce] '[expr', -ec0] '[expr', -0.1*fce] '[expr', -esp] lamda ftU Ets);       # Cover 'concrete', (unconfined)
       # print("uniaxialMaterial Concrete02 IDconcCover, (-fce), (-ec0), (-0.1*fce), (-esp) lamda ftU Ets;       # Cover concrete (unconfined)")

       ##### PRINT LINE TO CHECK #####
       # print("Cover fpc = (-fce), epsc0 =, (-Efact*ec0)")

       model.uniaxialMaterial('Concrete02', IDconcCore, (-fcc) '[expr', -ecc] '[expr', -fcuu] '[expr', -ecuu] lamda ftC Ets);       # Core 'concrete', (confined)
       # print("uniaxialMaterial Concrete02 IDconcCore, (-fcc), (-ecc), (-fcuu), (-ecuu) lamda ftC Ets;       # Core concrete (confined)")

       # Build Octagonal RC Column Section
       numSubdivCirc = 64; # Determines # of fibers in the circumferential direction, around the entire circumference
       ColMatTag = 10*ColSecTag;    # Column material cross-section tag.  Set equal to 10*Column's fiber model.element(tag.)
       # Define elastic torsional stiffness: Based on SDC-1.6 sec. 5.6.2 and sec. 1.1,
       # fundamental period less than 0.7 sec, reduction is required
       # *****CHECK FUNDAMENTAL PERIOD TO SEE IF REDUCTION STILL REQUIRED
       model.uniaxialMaterial('Elastic', IDTorsion '[expr',  0.2*Gc*Jcol]);
       section.Fiber(ColMatTag torsion=IDTorsion, , [
              patch.circ(IDconcCore numSubdivCirc 5 0. 0., ),(0.000e+00*in) '[expr', Rcore/2] 0.000e+00 3.600e+02; # Inner 'Core', Patch, 5 radial 'fibers'
              patch.circ(IDconcCore numSubdivCirc 10 0. 0., ),(Rcore/2) '[expr', Rcore] 0.000e+00 3.600e+02; # Outer 'Core', Patch, 10 radial 'fibers'
              numSlices = 8;       # Determines # of slices in each of the 8 sections of the octagon
              numSubdivIJ = 1;       # Determines # of fibers in the circumferential direction of the cover patch
              numSubdivJK = 4;       # Determines # of fibers in the radial direction of the cover patch
              for {i = 0} {i < 8} {i += 1} { # For 'each', of 'the', 8 sections 'of', the 'octagon'
                     startAngle =       (i*pi/4+pi/8)
                     for {j = 0} {j < numSlices} {j += 1} {
                            phi =              (pi/4/numSlices)
                            sita1 = startAngle + j*phi;       # Slice start angle
                            sita2 = sita1 + phi;                     # Slice end angle
                            yI =              (Rcore*cos(sita2))
                            zI =              (Rcore*sin(sita2))
                            yJ =              (Rcore*cos(sita1))
                            zJ =              (Rcore*sin(sita1))
                            oR1 =              (Rcol/cos(pi/8 - j*phi))
                            oR2 =              (Rcol/cos(pi/8 - (j+1)*phi))
                            yK =              (oR1*cos(sita1))
                            zK =              (oR1*sin(sita1))
                            yL =              (oR2*cos(sita2))
                            zL =              (oR2*sin(sita2))
                            patch.quad(IDconcCover numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal 'cover'


              layer.circ(IDSteel nLbar Along 0. 0. Rlong), # Longitudinal 'Bars'

       model.uniaxialMaterial('Elastic', IDShear , (9./10.)*Gc*Acol       )# Define 'elastic', shear 'stiffness'
       section 'Aggregator', ColSecTag IDShear 'Vy', IDShear 'Vz', section=ColMatTag, ;


#---------------------------------------------------------------------
# From ColSectionOctWideIEu.tcl
#---------------------------------------------------------------------
def BuildWideOctColSectionInelMod(ColSecTag  Dcol  Wcol  nLbar  nLbar2  DLbar  DLbar2  sTbar):
  # OCTAGONAL WIDE COLUMN SECTION DEFINITION PROCEDURE - INELASTIC PROPERTIES WITH MODIFIED CRUSHING STRESS/STRAIN PAIR

  ################### Octagonal WIDE (Interlocking) Column Fiber Cross Section Assignment ###################

  # ColSecTag                           -> Column's fiber element tag
  # Dcol                                  -> Shorter Width of octagonal column (to flat sides)
  # Wcol                                  -> Longer Width of octagonal column (to flat sides)
  # nLbar                                  -> Number of main (outer) longitudinal bars
  # nLbar2                           -> Number of secondary (inner) longitudinal bars
  # DLbar                                   -> Diameter of main (outer) longitudinal bars
  # DLbar2                             -> Diameter of secondary (inner) longitudinal bars
  # sTbar                             -> Spacing of transverse spiral reinforcement
  # Local Y axis is the horizontal axis.  Local Z axis is the vertical axis.  The longer width of the column is in the local Y direction.
       ##### PRINT LINE TO CHECK #####
       # print("Building WIDE octagonal column section for column ColSecTag")

       ftfact = 1.0;
       fcfact = 1.0;

       global 'kips', in 'ksi_psi', pi 'fce', ec0 esp 'Ec', Gc 'Es', Esh 'fy', fu 'esh', esu; # Read 'in', unit 'and', constant 'variables', as 'well', as 'material', property 'variables'

       IDconcCore =  (ColSecTag*10+1)
       IDconcCover = (ColSecTag*10+2)
       IDSteel =     (ColSecTag*10+3)
       IDShear =     (ColSecTag*10+4)
       IDTorsion =   (ColSecTag*10+5)

       # Column component dimensions
       tcover = 2.0*inch                                                     # 2 inch cover width
       Rcol = Dcol/2.0;                                                        # Radius of octagonal column (to flat sides)
       RcolDiag = Dcol*sqrt(4.0+2.0*sqrt(2.0))/(2.0+2.0*sqrt(2.0));         # Radius of octagonal column (to corners)
       Dcore = Dcol-2.0*tcover;                                          # Diameter of circular core
       Dcore2 = Wcol-2.0*tcover;                                          # Longer width of core
       Rcore = Dcore/2.0;                                                        # Radius of circular core
       Along = pi*DLbar**2/4.0;                                          # Area of OUTER longitudinal reinforcement bar
       Along2 = pi*DLbar2**2/4.0;                                          # Area of INNER longitudinal reinforcement bar
       DTbar = 0.625*inch;                                                      # Diameter of transverse spiral reinforcement bar (#5 Rebar)
       Asp = pi*DTbar**2/4.0;                                          # Area of transverse spiral reinforcement bar
       Dtran = Dcore-DTbar;                                                 # Diameter of spiral of transverse spiral reinforcement
       Dtran2 = Dcore2-DTbar;                                                  # Longer width of transverse spiral reinforcement
       rho = 4.0*Asp/(Dtran*sTbar);                            # Density of transverse spiral reinforcement in shorter direction
       rho2 = 4.0*Asp/(Dtran2*sTbar);                            # Density of transverse spiral reinforcement in longer direction
       Dlong = Dcore-2*DTbar-DLbar;                                   # Diameter of ring of longitudinal reinforcement
       Rlong = Dlong/2.0;                                                        # Radius of ring of longitudinal reinforcement

       # Section Area Properties
       spO = (Wcol-Dcol)/2.0;                                                         # Offof = octagonal sections from centroid (along horizontal axis)
       Acol = (2*Dcol**2)/(1+sqrt(2)) + (Wcol-Dcol)*Dcol; # Area of WIDE octagonal column section
       I2col = 0.6381*RcolDiag**4;                                       # Second moment of inertia, horizontal axis, for WIDE octagonal column section
       I3col = I2col+Acol*spO**2+(Wcol-Dcol)*Dcol**3/12; # Second moment of inertia, vertical axis, for WIDE octagonal column section
       Jcol = I2col+I3col;                                              # Polar moment of inertia for WIDE octagonal column section

       # Material Definition - Steel (Steel02)
       R0 = 18                                                       # control the transition from elastic to plastic branches
       cR1 = 0.925;                                                        # control the transition from elastic to plastic branches
       cR2 = 0.15;                                                        # control the transition from elastic to plastic branches
       model.uniaxialMaterial('Steel02', IDSteel fy Es 0.02 R0 cR1 cR2)

       # Material Definition - Concrete (Concrete02)
       # Compressive Properties
       ke = 1-sTbar/Dtran;                                   # Effective confinement strength coefficient from transverse reinforcement, shorter direction
       ke2 = 1-sTbar/Dtran2;                            # Effective confinement strength coefficient from transverse reinforcement, longer direction
       f2e = ke*rho*fy/2;                                    # Effective confinement strength in the shorter direction
       f3e = ke2*rho2*fy/2;                             # Effective confinement strength in the longer direction
       fcc2 = fce*(-1.254+2.254*sqrt(1+7.94*f2e/fce)-2*f2e/fce);       # Reference core compressive strength for f2e f2e -compression
       fcc3 = fce*(-1.254+2.254*sqrt(1+7.94*f3e/fce)-2*f3e/fce);       # Reference core compressive strength for f3e f3e -compression
       fcc = 9.05                                                       # Actual core compressive strength for f2e f3e

       ##### PRINT LINE TO CHECK #####
       # print("fce=fce, f2e=f2e, f3e=f3e, f2e/fce=[expr f2e/fce], f3e/fce=[expr f3e/fce], fcc2=fcc2, fcc3=fcc3, fcc=fcc")

       fcc = fce*(-1.254+2.254*sqrt(1+7.94*f3e/fce)-2*f3e/fce);       # Core compressive strength -compression
       ecc = ec0*(1.0+5.0*(fcc/fce-1.0));       # Core strain at maximum strength -compression
       ecu = 0.004+f3e/(4.0*fce);                     # Core crushing (ultimate) strain -compression
       lamda = 0.1                                                       # Ratio between unloading slope at eps2 and initial slope Ec
       xu = ecu/ecc;                                          # Mander equation parameter x
       ru = Ec/(Ec-fcc/ecc);                      # Mander equation parameter r
       fcu = fcfact*fcc*xu*ru/(ru-1+xu**ru);       # Core crushing (ultimate) strength -compression
       fcuu = 0.15*fcc;                                           # Modified Core crushing (ultimate) strength, equal to 0.15x core compressive strength -compression
       ecuu = ecu+(ecu-ecc)*(fcuu-fcu)/(fcu-fcc); # Modified Core crushing (ultimate) strain, corresponding to fcuu, calculated by extrapolating the line between (ecc,fcc) and (ecu,fcu) -compression

       ##### PRINT LINE TO CHECK #####
       # print("f3e=f3e, fcc=fcc, ecc=ecc, ecu=ecu, fcu=fcu, ecuu=ecuu, fcuu=fcuu");
       # print("xu=xu, ru=ru, fcu=fcu");

       # Tensile Properties
       ftU = 0;       # Cover tensile strength +tension
       ftC = 0;       # Core tensile strength +tension
       Ets = 0;       # Tension softening stiffness

       ##### PRINT LINE TO CHECK #####
       # print("ftU=ftU, ftC=ftC, Ets=Ets")

       # UniaxialMaterial Definition
       model.uniaxialMaterial('Concrete02', IDconcCover '[expr', -fce] '[expr', -ec0] '[expr', -0.1*fce] '[expr', -esp] lamda ftU Ets);       # Cover 'concrete', (unconfined)
       # print("uniaxialMaterial Concrete02 IDconcCover, (-fce), (-ec0), (-0.1*fce), (-esp) lamda ftU Ets;       # Cover concrete (unconfined)")
       ##### PRINT LINE TO CHECK #####
       # print("Cover fpc = (-fce), epsc0 =, (-Efact*ec0)")

       model.uniaxialMaterial('Concrete02', IDconcCore, (-fcc) '[expr', -ecc] '[expr', -fcuu] '[expr', -ecuu] lamda ftC Ets);       # Core 'concrete', (confined)
       # print("uniaxialMaterial Concrete02 IDconcCore, (-fcc), (-ecc), (-fcuu), (-ecuu) lamda ftC Ets;       # Core concrete (confined)")
       ##### PRINT LINE TO CHECK #####
       # print("Core fpc = (-fcc), epsc0 =, (-Efact*ecc)")

       # Build Octagonal RC Column Section
       numSubdivCirc = 64; # Determines # of fibers in the circumferential direction, around the entire circumference
       ColMatTag = 10*ColSecTag;    # Column material cross-section tag.  Set equal to 10*Column's fiber model.element(tag.)
       # Define elastic torsional stiffness: Based on SDC-1.6 sec. 5.6.2 and sec. 1.1,
       # fundamental period less than 0.7 sec, reduction is required
       # *****CHECK FUNDAMENTAL PERIOD TO SEE IF REDUCTION STILL REQUIRED
       model.uniaxialMaterial('Elastic', IDTorsion '[expr',  0.2*Gc*Jcol]);
       section.Fiber(ColMatTag torsion=IDTorsion, , [
              patch.circ(IDconcCore numSubdivCirc 5, ),(-spO)  0., (0.000e+00*in) '[expr', Rcore/2] 90.0 270.0; # LEFT 'Inner', circular 'Core', Patch, 5 radial 'fibers'
              patch.circ(IDconcCore numSubdivCirc 5, ),(spO) 0., (0.000e+00*in) '[expr', Rcore/2] 270.0 450.0; # RIGHT 'Inner', circular 'Core', Patch, 5 radial 'fibers'
              patch.circ(IDconcCore numSubdivCirc 10, ),(-spO) 0., (Rcore/2) '[expr', Rcore] 90.0 270.0; # LEFT 'Outer', circular 'Core', Patch, 10 radial 'fibers'
              patch.circ(IDconcCore numSubdivCirc 10, ),(spO) 0., (Rcore/2) '[expr', Rcore] 270.0 450.0; # RIGHT 'Outer', circular 'Core', Patch, 10 radial 'fibers'
              patch.quad(IDconcCore 5 5, ),(-spO) '[expr', -Rcore/2] '[expr', spO] '[expr', -Rcore/2] '[expr', spO] '[expr', Rcore/2] '[expr', -spO] '[expr', Rcore/2]; # Inner 'rectangular', Core Patch, 5 fibers 'in', y & z
              patch.quad(IDconcCore 10 10, ),(-spO) '[expr', -Rcore] '[expr', spO] '[expr', -Rcore] '[expr', spO] '[expr', -Rcore/2] '[expr', -spO] '[expr', -Rcore/2]; # BOTTOM 'Outer', rectangular 'Core', Patch, 10 fibers 'in', y & z
              patch.quad(IDconcCore 10 10, ),(-spO) '[expr', Rcore/2] '[expr', spO] '[expr', Rcore/2] '[expr', spO] '[expr', Rcore] '[expr', -spO] '[expr', Rcore]; # TOP 'Outer', rectangular 'Core', Patch, 10 fibers 'in', y & z
              patch.quad(IDconcCover 10 4, ),(-spO) '[expr', -Rcol] '[expr', spO] '[expr', -Rcol] '[expr', spO] '[expr', -Rcore] '[expr', -spO] '[expr', -Rcore]; # BOTTOM 'Outer', rectangular 'COVER', Patch, 10 fibers 'in', y, 2 fibers 'in', z
              patch.quad(IDconcCover 10 4, ),(-spO) '[expr', Rcore] '[expr', spO] '[expr', Rcore] '[expr', spO] '[expr', Rcol] '[expr', -spO] '[expr', Rcol]; # TOP 'Outer', rectangular 'COVER', Patch, 10 fibers 'in', y, 2 fibers 'in', z
              numSlices = 8;       # Determines # of slices in each of the 8 sections of the octagon
              numSubdivIJ = 1;       # Determines # of fibers in the circumferential direction of the cover patch
              numSubdivJK = 4;       # Determines # of fibers in the radial direction of the cover patch
              for {i = 0} {i < 2} {i += 1} { # For 'the', first 2 sections 'of', the 'octagon', (TOP RIGHT; offset 'all', y 'coordinates', by +spO)
                     startAngle =       (i*pi/4+pi/8)
                     for {j = 0} {j < numSlices} {j += 1} {
                            phi =              (pi/4/numSlices)
                            sita1 = startAngle + j*phi;       # Slice start angle
                            sita2 = sita1 + phi;                     # Slice end angle
                            yI =              (Rcore*cos(sita2)+spO)
                            zI =              (Rcore*sin(sita2))
                            yJ =              (Rcore*cos(sita1)+spO)
                            zJ =              (Rcore*sin(sita1))
                            oR1 =              (Rcol/cos(pi/8 - j*phi))
                            oR2 =              (Rcol/cos(pi/8 - (j+1)*phi))
                            yK =              (oR1*cos(sita1)+spO)
                            zK =              (oR1*sin(sita1))
                            yL =              (oR2*cos(sita2)+spO)
                            zL =              (oR2*sin(sita2))
                            patch.quad(IDconcCover numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal 'cover'


              for {i = 2} {i < 6} {i += 1} { # For 'the', 3rd 'through', 6th 'sections', of 'the', octagon (TOP & BOTTOM LEFT; offset 'all', y 'coordinates', by -spO)
                     startAngle =       (i*pi/4+pi/8)
                     for {j = 0} {j < numSlices} {j += 1} {
                            phi =              (pi/4/numSlices)
                            sita1 = startAngle + j*phi;       # Slice start angle
                            sita2 = sita1 + phi;                     # Slice end angle
                            yI =              (Rcore*cos(sita2)-spO)
                            zI =              (Rcore*sin(sita2))
                            yJ =              (Rcore*cos(sita1)-spO)
                            zJ =              (Rcore*sin(sita1))
                            oR1 =              (Rcol/cos(pi/8 - j*phi))
                            oR2 =              (Rcol/cos(pi/8 - (j+1)*phi))
                            yK =              (oR1*cos(sita1)-spO)
                            zK =              (oR1*sin(sita1))
                            yL =              (oR2*cos(sita2)-spO)
                            zL =              (oR2*sin(sita2))
                            patch.quad(IDconcCover numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal 'cover'


              for {i = 6} {i < 8} {i += 1} { # For 'the', last 2 sections 'of', the 'octagon', (BOTTOM RIGHT; offset 'all', y 'coordinates', by +spO)
                     startAngle =       (i*pi/4+pi/8)
                     for {j = 0} {j < numSlices} {j += 1} {
                            phi =              (pi/4/numSlices)
                            sita1 = startAngle + j*phi;       # Slice start angle
                            sita2 = sita1 + phi;                     # Slice end angle
                            yI =              (Rcore*cos(sita2)+spO)
                            zI =              (Rcore*sin(sita2))
                            yJ =              (Rcore*cos(sita1)+spO)
                            zJ =              (Rcore*sin(sita1))
                            oR1 =              (Rcol/cos(pi/8 - j*phi))
                            oR2 =              (Rcol/cos(pi/8 - (j+1)*phi))
                            yK =              (oR1*cos(sita1)+spO)
                            zK =              (oR1*sin(sita1))
                            yL =              (oR2*cos(sita2)+spO)
                            zL =              (oR2*sin(sita2))
                            patch.quad(IDconcCover numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal 'cover'


              layer.circ(IDSteel '[expr', nLbar),/2] Along '[expr', -spO] 0. Rlong 45.0 315.0; # LEFT 'outer', Longitudinal 'Bars'
              layer.circ(IDSteel '[expr', nLbar),/2] Along '[expr', spO] 0. Rlong 225.0 495.0; # RIGHT 'outer', Longitudinal 'Bars'
              layer.circ(IDSteel '[expr', nLbar2),/2] Along2, (-spO) 0. Rlong 315.0 405.0; # RIGHT 'inner', Longitudinal 'Bars'
              layer.circ(IDSteel '[expr', nLbar2),/2] Along2, (spO) 0. Rlong 135.0 225.0; # LEFT 'inner', Longitudinal 'Bars'

       model.uniaxialMaterial('Elastic', IDShear , (9./10.)*Gc*Acol       )# Define 'elastic', shear 'stiffness'
       section 'Aggregator', ColSecTag IDShear 'Vy', IDShear 'Vz', section=ColMatTag, ;


#---------------------------------------------------------------------
# From ColSectionCircIEu.tcl
#---------------------------------------------------------------------
def BuildCircColSectionInelMod(ColSecTag  Dcol  nLbar  DLbar  sTbar):
    # CIRCULAR COLUMN SECTION DEFINITION PROCEDURE - INELASTIC PROPERTIES WITH MODIFIED CRUSHING STRESS/STRAIN PAIR

    ################### Circular Column Fiber Cross Section Assignment ###################

    # ColSecTag                           -> Column's fiber element tag
    # Dcol                                  -> Column diameter
    # nLbar                                  -> Number of longitudinal bars
    # DLbar                                   -> Diameter of longitudinal bars
    # sTbar                             -> Spacing of transverse spiral reinforcement

       ##### PRINT LINE TO CHECK #####
       # print("Building circular column section for column ColSecTag")

       ftfact = 1.0;
       fcfact = 1.0;

       global 'kips', in 'ksi_psi', pi 'fce', ec0 esp 'Ec', Gc 'Es', Esh 'fy', fu 'esh', esu; # Read 'in', unit 'and', constant 'variables', as 'well', as 'material', property 'variables'

       IDconcCore =  (ColSecTag*10+1)
       IDconcCover = (ColSecTag*10+2)
       IDSteel =     (ColSecTag*10+3)
       IDShear =     (ColSecTag*10+4)
       IDTorsion =   (ColSecTag*10+5)

       # Column component dimensions
       tcover = 2.0*inch                                 # 2 inch cover width
       Rcol = Dcol/2.0;                                  # Radius of column
       Dcore = Dcol-2.0*tcover;                          # Diameter of core
       Rcore = Dcore/2.0;                                # Radius of core
       Along = pi*DLbar**2/4.0;                          # Area of longitudinal reinforcement bar
       DTbar = 0.625*inch;                               # Diameter of transverse spiral reinforcement bar (#5 Rebar)
       Asp = pi*DTbar**2/4.0;                            # Area of transverse spiral reinforcement bar
       Dtran = Dcore-DTbar;                              # Diameter of spiral of transverse spiral reinforcement
       rho = 4.0*Asp/(Dtran*sTbar);                      # Density of transverse spiral reinforcement
       Dlong = Dcore-2*DTbar-DLbar;                      # Diameter of ring of longitudinal reinforcement
       Rlong = Dlong/2.0;                                # Radius of ring of longitudinal reinforcement

       # Section Area Properties
       Acol = (2.0*Dcol**2)/(1.0+sqrt(2.0));                     # Area of column section
       Jcol = pi*Rcol**4/2.0;                                    # Polar moment of inertia for column section
       I3col = pi*Rcol**4/4.0;                                   # Second moment of inertia, 1st transverse direction, for column section
       I2col = pi*Rcol**4/4.0;                                   # Second moment of inertia, 2nd transverse direction, for column section

       # Material Definition - Steel (Steel02)
       R0 = 18                                                       # control the transition from elastic to plastic branches
       cR1 = 0.925;                                                  # control the transition from elastic to plastic branches
       cR2 = 0.15;                                                   # control the transition from elastic to plastic branches
       model.uniaxialMaterial('Steel02', IDSteel fy Es 0.02 R0 cR1 cR2)

       # Material Definition - Concrete (Concrete02)
       # Compressive Properties
       ke = 1-sTbar/Dtran;                                   # Effective confinement strength coefficient from transverse reinforcement
       f3e = ke*rho*fy/2;                                    # Effective confinement strength
       fcc = fce*(-1.254+2.254*sqrt(1+7.94*f3e/fce)-2*f3e/fce);       # Core compressive strength -compression
       ecc = ec0*(1.0+5.0*(fcc/fce-1.0));       # Core strain at maximum strength -compression
       ecu = 0.004+f3e/(4.0*fce);                     # Core crushing (ultimate) strain -compression
       lamda = 0.1                                                       # Ratio between unloading slope at eps2 and initial slope Ec
       xu = ecu/ecc;                                          # Mander equation parameter x
       ru = Ec/(Ec-fcc/ecc);                      # Mander equation parameter r
       fcu = fcfact*fcc*xu*ru/(ru-1+xu**ru);       # Core crushing (ultimate) strength -compression
       fcuu = 0.15*fcc;                                           # Modified Core crushing (ultimate) strength, equal to 0.15x core compressive strength -compression
       ecuu = ecu+(ecu-ecc)*(fcuu-fcu)/(fcu-fcc); # Modified Core crushing (ultimate) strain, corresponding to fcuu, calculated by extrapolating the line between (ecc,fcc) and (ecu,fcu) -compression

       ##### PRINT LINE TO CHECK #####
       # print("f3e=f3e, fcc=fcc, ecc=ecc, ecu=ecu, fcu=fcu, ecuu=ecuu, fcuu=fcuu");
       # print("xu=xu, ru=ru, fcu=fcu");

       # Tensile Properties
       ftU = 0;       # Cover tensile strength +tension
       ftC = 0;       # Core tensile strength +tension
       Ets = 0;       # Tension softening stiffness

       ##### PRINT LINE TO CHECK #####
       # print("ftU=ftU, ftC=ftC, Ets=Ets")

       # UniaxialMaterial Definition
       model.uniaxialMaterial('Concrete02', IDconcCover '[expr', -fce] '[expr', -ec0] '[expr', -0.1*fce] '[expr', -esp] lamda ftU Ets);       # Cover 'concrete', (unconfined)
       # print("uniaxialMaterial Concrete02 IDconcCover, (-fce), (-ec0), (-0.1*fce), (-esp) lamda ftU Ets;       # Cover concrete (unconfined)")

       ##### PRINT LINE TO CHECK #####
       # print("Cover fpc = (-fce), epsc0 =, (-Efact*ec0)")

       model.uniaxialMaterial('Concrete02', IDconcCore, (-fcc) '[expr', -ecc] '[expr', -fcuu] '[expr', -ecuu] lamda ftC Ets);       # Core 'concrete', (confined)
       # print("uniaxialMaterial Concrete02 IDconcCore, (-fcc), (-ecc), (-fcuu), (-ecuu) lamda ftC Ets;       # Core concrete (confined)")

       # Build Circular RC Column Section
       numSubdivCirc = 32; # Determines # of fibers in the circumferential direction, around the entire circumference, for the core fibers
       ColMatTag = 10*ColSecTag;    # Column material cross-section tag.  Set equal to 10*Column's fiber model.element(tag.)
       # Define elastic torsional stiffness: Based on SDC-1.6 sec. 5.6.2 and sec. 1.1,
       # fundamental period less than 0.7 sec, reduction is required
       # *****CHECK FUNDAMENTAL PERIOD TO SEE IF REDUCTION STILL REQUIRED
       model.uniaxialMaterial('Elastic', IDTorsion '[expr',  0.2*Gc*Jcol]);
       section.Fiber(ColMatTag torsion=IDTorsion, , [
              patch.circ(IDconcCore numSubdivCirc 10 0. 0. 0.0, ),(Rcore) 0.000e+00 3.600e+02; # Core Patch, 10 radial 'fibers'
              patch.circ(IDconcCover numSubdivCirc 4 0. 0., ),(Rcore) '[expr', Rcol] 0.000e+00 3.600e+02; # Cover Patch, 2 radial 'fibers'
              layer.circ(IDSteel nLbar Along 0. 0. Rlong), # Longitudinal 'Bars'

       model.uniaxialMaterial('Elastic', IDShear , (9./10.)*Gc*Acol       )# Define 'elastic', shear 'stiffness'
       section 'Aggregator', ColSecTag IDShear 'Vy', IDShear 'Vz', section=ColMatTag, ;


#---------------------------------------------------------------------
# From ColSectionOctIEPin.tcl
#---------------------------------------------------------------------
def BuildOctColPINSection(ColSecTag,  Dcol):
       # OCTAGONAL COLUMN PIN SECTION DEFINITION PROCEDURE
       # -- *UNREINFORCED* INELASTIC PROPERTIES

       ########## Octagonal Column Fiber Cross Section Assignment #########

       # ColSecTag                           -> Column's fiber element tag
       # Dcol                                  -> Width of octagonal column (to flat sides)

       ##### PRINT LINE TO CHECK #####
       # print("Building octagonal column section for column ColSecTag")

       ftfact = 1.0;
       fcfact = 1.0;

       global 'kips', in 'ksi_psi', pi 'fce', ec0 esp 'Ec', Gc 'Es', Esh 'fy', fu 'esh', esu;
       # Read 'in', unit 'and', constant 'variables', as 'well', as 'material', property 'variables'

       IDconc =               (ColSecTag*10+2)
       IDSteel =     (ColSecTag*10+3)
       IDShear =     (ColSecTag*10+4)
       IDTorsion =   (ColSecTag*10+5)

       # Column component dimensions
       tcover = 2.0*inch                                              # 2 inch cover width
       Rcol = Dcol/2.0;                                               # Radius of octagonal column (to flat sides)
       RcolDiag = Dcol*sqrt(4.0+2.0*sqrt(2.0))/(2.0+2.0*sqrt(2.0));   # Radius of octagonal column (to corners)
       Dcore = Dcol-2.0*tcover;                                       # Diameter of circular core
       Rcore = Dcore/2.0;                                             # Radius of circular core

       # Section Area Properties
       Acol = (2.0*Dcol**2)/(1.0+sqrt(2.0));    # Area of octagonal column section
       Jcol = 1.2762*RcolDiag**4;               # Polar moment of inertia for octagonal column section
       I3col = 0.6381*RcolDiag**4;              # Second moment of inertia, 1st transverse direction, for octagonal column section
       I2col = 0.6381*RcolDiag**4;              # Second moment of inertia, 2nd transverse direction, for octagonal column section

       # Material Definition - Steel (Steel02)
       R0 = 18                                  # control the transition from elastic to plastic branches
       cR1 = 0.925;                             # control the transition from elastic to plastic branches
       cR2 = 0.15;                              # control the transition from elastic to plastic branches
       model.uniaxialMaterial('Steel02', IDSteel fy Es 0.02 R0 cR1 cR2)

       # Concrete Properties
       lamda = 0.1;                            # Ratio between unloading slope at eps2 and initial slope Ec
       ftU = 0;                                 # Tensile strength +tension
       Ets = 0;                                 # Tension softening stiffness

       ##### PRINT LINE TO CHECK #####
       # print("ftU=ftU, Ets=Ets")

       # UniaxialMaterial Definition
       model.uniaxialMaterial('Concrete02', IDconc '[expr', -fce] '[expr', -ec0] '[expr', -0.1*fce] '[expr', -esp] lamda ftU Ets);       # Concrete (unconfined)

       # Build Octagonal URC Column Section
       numSubdivCirc = 64; # Determines # of fibers in the circumferential direction, around the entire circumference
       ColMatTag = 10*ColSecTag;    # Column material cross-section tag.  Set equal to 10*Column's fiber model.element(tag.)
       # Define elastic torsional stiffness: Based on SDC-1.6 sec. 5.6.2 and sec. 1.1,
       # fundamental period less than 0.7 sec, reduction is required
       # *****CHECK FUNDAMENTAL PERIOD TO SEE IF REDUCTION STILL REQUIRED
       model.uniaxialMaterial('Elastic', IDTorsion '[expr',  0.2*Gc*Jcol]);
       section.Fiber(ColMatTag torsion=IDTorsion, , [
              patch circ IDconc numSubdivCirc 5 0. 0., (0.000e+00*in), (Rcore/2) 0.000e+00 3.600e+02; # Inner "Core" Patch, 5 radial fibers
              patch circ IDconc numSubdivCirc 10 0. 0., (Rcore/2), (Rcore) 0.000e+00 3.600e+02; # Outer "Core" Patch, 10 radial fibers
              numSlices = 8;       # Determines # of slices in each of the 8 sections of the octagon
              numSubdivIJ = 1;       # Determines # of fibers in the circumferential direction of the cover patch
              numSubdivJK = 4;       # Determines # of fibers in the radial direction of the cover patch
              for {i = 0} {i < 8} {i += 1} { # For 'each', of 'the', 8 sections 'of', the 'octagon'
                  startAngle =       (i*pi/4+pi/8)
                  for {j = 0} {j < numSlices} {j += 1} {
                      phi =              (pi/4/numSlices)
                      sita1 = startAngle + j*phi;       # Slice start angle
                      sita2 = sita1 + phi;                     # Slice end angle
                      yI  =             (Rcore*cos(sita2))
                      zI  =             (Rcore*sin(sita2))
                      yJ  =             (Rcore*cos(sita1))
                      zJ  =             (Rcore*sin(sita1))
                      oR1 =             (Rcol/cos(pi/8 - j*phi))
                      oR2 =             (Rcol/cos(pi/8 - (j+1)*phi))
                      yK  =             (oR1*cos(sita1))
                      zK  =             (oR1*sin(sita1))
                      yL  =             (oR2*cos(sita2))
                      zL  =             (oR2*sin(sita2))
                      patch.quad(IDconc numSubdivIJ numSubdivJK yI zI yJ zJ yK zK yL zL), # Cover 'Patch', connects 'the', circular 'core', to 'the', octagonal 'cover'


       layer.circ(IDSteel 8 3.0 0. 0. 6.0), # Longitudinal 'Bars'

       model.uniaxialMaterial('Elastic', IDShear , (9./10.)*Gc*Acol       )# Define 'elastic', shear 'stiffness'
       section 'Aggregator', ColSecTag IDShear 'Vy', IDShear 'Vz', section=ColMatTag, ;

