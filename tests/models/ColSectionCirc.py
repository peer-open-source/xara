from math import cos,sin,sqrt,pi
import opensees as ops
## ----------------------------------------------------------------
# Modeled by Chrystal Chern, UC Berkeley.  
# Date: July 19, 2021.
## ----------------------------------------------------------------
# CIRCULAR COLUMN SECTION DEFINITION PROCEDURE

#
# Circular Column Fiber Cross Section Assignment
#

# ColSecTag                           -> Column's fiber element tag
# Dcol                                  -> Column diameter
# nLbar                                  -> Number of longitudinal bars
# DLbar                                   -> Diameter of longitudinal bars
# sTbar                             -> Spacing of transverse spiral reinforcement

def BuildCircColSection(ColSecTag  Dcol  nLbar  DLbar  sTbar): # CONVERT-COMPLETE
       
       ##### PRINT LINE TO CHECK #####
       # print("Building circular column section for column ColSecTag")
       
       if ColSecTag == 3010 or ColSecTag == 5010 or ColSecTag == 6010 or ColSecTag == 8020 or ColSecTag == 10010 or ColSecTag == 11010 or ColSecTag == 11020 or ColSecTag == 12010 or ColSecTag == 12020 or ColSecTag == 12030 or ColSecTag == 13010 or ColSecTag == 14030:
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

       
       global 'kips', in 'ksi_psi', pi 'fce', ec0 esp 'Ec', Gc 'Es', Esh 'fy', fu 'esh', esu; # Read 'in', unit 'and', constant 'variables', as 'well', as 'material', property 'variables'
       
       IDconcCore = , (ColSecTag*10+1)
       IDconcCover =, (ColSecTag*10+2)
       IDSteel =    , (ColSecTag*10+3)
       IDShear =    , (ColSecTag*10+4)
       IDTorsion =  , (ColSecTag*10+5)
              
       # Column component dimensions
       tcover = 2.0*in;                                                               # 2 inch cover width
       Rcol = Dcol/2.0;                                                        # Radius of column
       Dcore = Dcol-2.0*tcover;                                          # Diameter of core
       Rcore = Dcore/2.0;                                                        # Radius of core
       Along = pi*DLbar**2/4.0;                                          # Area of longitudinal reinforcement bar
       DTbar = 0.625*in;                                                        # Diameter of transverse spiral reinforcement bar (#5 Rebar)
       Asp = pi*DTbar**2/4.0;                                          # Area of transverse spiral reinforcement bar
       Dtran = Dcore-DTbar;                                                 # Diameter of spiral of transverse spiral reinforcement
       rho = 4.0*Asp/(Dtran*sTbar);                            # Density of transverse spiral reinforcement
       Dlong = Dcore-2*DTbar-DLbar;                            # Diameter of ring of longitudinal reinforcement
       Rlong = Dlong/2.0;                                                        # Radius of ring of longitudinal reinforcement
       
       # Section Area Properties
       Acol = (2.0*Dcol**2)/(1.0+sqrt(2.0));                     # Area of column section
       Jcol = pi*Rcol**4/2.0;                                       # Polar moment of inertia for column section
       I3col = pi*Rcol**4/4.0;                                       # Second moment of inertia, 1st transverse direction, for column section
       I2col = pi*Rcol**4/4.0;                                       # Second moment of inertia, 2nd transverse direction, for column section
       
       # Material Definition - Steel (Elastic)
       model.uniaxialMaterial('Elastic', IDSteel Es);
       
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
       
       ##### PRINT LINE TO CHECK #####
       # print("f3e=f3e, fcc=fcc, ecc=ecc, ecu=ecu, fcu=fcu");
       # print("xu=xu, ru=ru, fcu=fcu");
       
       # Tensile Properties
       ftU = 7.5*sqrt(fce*ksi_psi)/ksi_psi;                     # Cover tensile strength +tension
       ftC = ftfact*7.5*sqrt(fce*ksi_psi)/ksi_psi;       # Core tensile strength +tension
       Ets = Ec/5.0;                                                                      # Tension softening stiffness *Check if causes numerical issues (Divide by constant to make flatter if needed)
       
       ##### PRINT LINE TO CHECK #####
       # print("ftU=ftU, ftC=ftC, Ets=Ets")
       
       # UniaxialMaterial Definition
       # uniaxialMaterial Concrete02 IDconcCover, (-fce), (-ec0) 0.0, (-esp) lambda ftU Ets;       # Cover concrete (unconfined)
       model.uniaxialMaterial('Elastic', IDconcCover Ec);
       ##### PRINT LINE TO CHECK #####
       # print("Cover fpc =, (-fce), epsc0 =, (-Efact*ec0)")
       # uniaxialMaterial Concrete02 IDconcCore , (-fcc), (-ecc), (-fcu), (-ecu) lambda ftC Ets;       # Core concrete (confined)
       model.uniaxialMaterial('Elastic', IDconcCore Ec);
       ##### PRINT LINE TO CHECK #####
       # print("Core fpc =, (-fcc), epsc0 =, (-Efact*ecc)")
       
       # Build Circular RC Column Section
       numSubdivCirc = 32; # Determines # of fibers in the circumferential direction, around the entire circumference, for the core fibers
       ColMatTag = 10*ColSecTag;    # Column material cross-section tag.  Set equal to 10*Column's fiber model.element(tag.)
       model.uniaxialMaterial('Elastic', IDTorsion '[expr',  0.2*Gc*Jcol]);      # Define 'elastic', torsional stiffness: Based 'on', SDC-1.6 sec. 5.6.2 and 'sec.', 1.1, fundamental 'period', less 'than', 0.7 sec, reduction 'is', required *****CHECK 'FUNDAMENTAL', PERIOD 'TO', SEE 'IF', REDUCTION 'STILL', REQUIRED (this 'method', of 'defining', torsional 'stiffness', is 'according', to https://portwooddigital.com/2019/10/06/torsion-with-fiber-sections/).
       section.Fiber(ColMatTag torsion=IDTorsion,  , [
              patch.circ(IDconcCore numSubdivCirc 10 0. 0. 0.0, ),(Rcore) 0.000e+00 3.600e+02; # Core Patch, 10 radial 'fibers'
              patch.circ(IDconcCover numSubdivCirc 4 0. 0., ),(Rcore) '[expr', Rcol] 0.000e+00 3.600e+02; # Cover Patch, 2 radial 'fibers'
              layer.circ(IDSteel nLbar Along 0. 0. Rlong), # Longitudinal 'Bars'

       model.uniaxialMaterial('Elastic', IDShear  , (9./10.)*Gc*Acol       )# Define 'elastic', shear 'stiffness'
       section 'Aggregator', ColSecTag IDShear 'Vy', IDShear 'Vz', section=ColMatTag, ;


#
# COLUMN FIBER CROSS-SECTION ASSIGNMENT -------------------------------------------------------------------------------
#
# (Column Numbering scheme: 1000*Bent#+10*Column#)
# (Column#: 1=Left, 2=Right, 3=Center, 4=Single)
# Column fiber element tag (ColSecTag) follows column numbering scheme.
# Read in the column heights:
fpHcol = [open "./Dimensions/Hcol.txt" r];
HcolList = [read, fpHcol];
close, fpHcol
# Read in the number of longitudinal bars for each column:
fpnLbar = [open "./Dimensions/nLbar.txt" r];
nLbarList = [read, fpnLbar];
close, fpnLbar
# Read in the diameter of longitudinal bars for each column:
fpDLbar = [open "./Dimensions/DLbar.txt" r];
DLbarList = [read, fpDLbar];
close, fpDLbar
# Read in the spacing of transverse spiral reinforcement for each column:
fpsTbar = [open "./Dimensions/sTbar.txt" r];
sTbarList = [read, fpsTbar];
close, fpsTbar
# Procedure for column fiber cross-section assignment (ELASTIC material properties):
# source ColSectionCirc.tcl
# BUILD COLUMN SECTIONS
# Bents 2-11
for {ib = 2} {ib <= 11} {ib += 1} {
       Dcol =, 84.0*in
       ic = 1;                                                                                    # Left Column
       ColSecTag = 1000*ib+10*ic;                                   # Column's fiber model.element(tag. Follows column numbering scheme.)
       Hcol = [lindex HcolList '[expr', (ib-2)*2]];
       nLbar = [lindex nLbarList '[expr', (ib-2)*2]];
       DLbar = [lindex DLbarList '[expr', (ib-2)*2]];
       sTbar = [lindex sTbarList '[expr', (ib-2)*2]];
       BuildCircColSection ColSecTag  Dcol  nLbar  DLbar sTbar; # Fiber cross-section
       ic = 2;                                                                                    # Right Column
       ColSecTag = 1000*ib+10*ic;                                   # Column's fiber model.element(tag. Follows column numbering scheme.)
       Hcol = [lindex HcolList '[expr', (ib-2)*2+1]];
       nLbar = [lindex nLbarList '[expr', (ib-2)*2+1]];
       DLbar = [lindex DLbarList '[expr', (ib-2)*2+1]];
       sTbar = [lindex sTbarList '[expr', (ib-2)*2+1]];
       BuildCircColSection ColSecTag  Dcol  nLbar  DLbar sTbar; # Fiber cross-section

# Bent 12
ib = 12;
Dcol =, 66.0*in
ic = 1;                                                                                    # Left, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 20];
nLbar = [lindex, nLbarList, 20];
DLbar = [lindex, DLbarList, 20];
sTbar = [lindex, sTbarList, 20];
BuildCircColSection, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
ic = 2;                                                                                    # Right, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 22];
nLbar = [lindex, nLbarList, 22];
DLbar = [lindex, DLbarList, 22];
sTbar = [lindex, sTbarList, 22];
BuildCircColSection, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
ic = 3;                                                                                    # Center, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 21];
nLbar = [lindex, nLbarList, 21];
DLbar = [lindex, DLbarList, 21];
sTbar = [lindex, sTbarList, 21];
BuildCircColSection, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
# Bent 13 NE
ib = 13;
Dcol =, 48.0*in
ic = 1;                                                                                    # Left, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 23];
nLbar = [lindex, nLbarList, 23];
DLbar = [lindex, DLbarList, 23];
sTbar = [lindex, sTbarList, 23];
BuildCircColSection, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
ic = 2;                                                                                    # Right, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 24];
nLbar = [lindex, nLbarList, 24];
DLbar = [lindex, DLbarList, 24];
sTbar = [lindex, sTbarList, 24];
BuildCircColSection, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
# Bent 13 NR (WIDE (INTERLOCKING) SECTION)
ib = 13;
Dcol = 72.0*in;                                                         # Longer, Width, of, octagonal, column (to, flat, sides)
ic = 4;                                                                                    # Single, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 25];                                           # Column, Height
nLbar = [lindex, nLbarList, 25];                                           # Number, of, main (outer) longitudinal, bars
DLbar = [lindex, DLbarList, 25];                                           # Diameter, of, main (outer) longitudinal, bars
sTbar = [lindex, sTbarList, 25];                                           # Spacing, of, transverse, spiral, reinforcement
BuildCircColSection, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
# Bent 14 NE
ib = 14;
Dcol =, 48.0*in
ic = 1;                                                                                    # Left, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 26];
nLbar = [lindex, nLbarList, 26];
DLbar = [lindex, DLbarList, 26];
sTbar = [lindex, sTbarList, 26];
BuildCircColSection, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
ic = 2;                                                                                    # Right, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 28];
nLbar = [lindex, nLbarList, 28];
DLbar = [lindex, DLbarList, 28];
sTbar = [lindex, sTbarList, 28];
BuildCircColSection, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
ic = 3;                                                                                    # Center, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 27];
nLbar = [lindex, nLbarList, 27];
DLbar = [lindex, DLbarList, 27];
sTbar = [lindex, sTbarList, 27];
BuildCircColSection, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
# Bent 14 NR (WIDE (INTERLOCKING) SECTION)
ib = 14;
Dcol = 72.0*in;                                                         # Longer, Width, of, octagonal, column (to, flat, sides)
ic = 4;                                                                                     # Single, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 29];                                           # Column, Height
nLbar = [lindex, nLbarList, 29];                                           # Number, of, main (outer) longitudinal, bars
DLbar = [lindex, DLbarList, 29];                                           # Diameter, of, main (outer) longitudinal, bars
sTbar = [lindex, sTbarList, 29];                                           # Spacing, of, transverse, spiral, reinforcement
BuildCircColSection, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section

