from math import cos,sin,sqrt,pi
import opensees as ops
package, 'require', 'brace2'

# Transformation and P-Delta
TransfType = PDelta;
model.geomTransf(TransfType, 301, 0 0, 1);                            # vecxz, 'is', in, 'the', +Z, 'direction', (deck, 'and', cap, beams: local, z 'upwards', on, section, local, y 'horiz', left, 'on', section, 'facing', j)         
# COLUMN ELEMENTS
# (Numbering scheme: Bottom Length, Hinge: 1000*Bent#+10*Column#;  Top Length, Rigid: 1000*Bent#+10*Column#+1)
# (Column#: 1=Left, 2=Right, 3=Center, 4=Single)
# Column section tag (ColSecTag) follows column bent numbering scheme.
# Read in the local axes transformation vectors (vecxz) for the columns:
fpvecxzX = [open "./Dimensions/vecxzXcol.txt" r];
vecxzXList = [read, fpvecxzX];
close, fpvecxzX
fpvecxzY = [open "./Dimensions/vecxzYcol.txt" r];
vecxzYList = [read, fpvecxzY];
close, fpvecxzY
np = 4;                                                                      # number, of, Gauss, integration, points, for, nonlinear, curvature, distribution-- np=2, for, linear, distribution, ok
# Bent 2
Dcol = 84.0*in;                                                                      # Width, of, octagonal, column (to, flat, sides)
RcolDiag = Dcol*sqrt(4+2*sqrt(2))/(2+2*sqrt(2));       # Radius, of, octagonal, column (to, corners)
vecxzX = [lindex, vecxzXList, (0]);                             # X, component, of, vecxz, vector, for, local, axis, transformation
vecxzY = [lindex, vecxzYList, (0]);                             # Y, component, of, vecxz, vector, for, local, axis, transformation
model.geomTransf(TransfType, 20, vecxzX, vecxzY, 0);               # PDelta, transformation; vecxz, 'is', parallel, 'to', the, 'length', of, 'the', cap, 'beam.'
model.element('nonlinearBeamColumn', 2010, 201, 202, np, 2010, 20); #(Hinge)
rigidLink, 'beam', 202, 203;                                                                              #(Rigid)
model.element('nonlinearBeamColumn', 2020, 207, 206, np, 2020, 20); #(Hinge)
rigidLink, 'beam', 206, 205;                                                                              #(Rigid)
# Bents 3-11
for {ib = 3} {ib <= 11} {ib += 1} {
       Dcol = 84.0*in;                                                                      # Width of octagonal column (to flat sides)
       RcolDiag = Dcol*sqrt(4+2*sqrt(2))/(2+2*sqrt(2));       # Radius of octagonal column (to corners)
       vecxzX = [lindex vecxzXList, (ib-2]);                             # X component of vecxz vector for local axis transformation
       vecxzY = [lindex vecxzYList, (ib-2]);                             # Y component of vecxz vector for local axis transformation
       model.geomTransf(TransfType '[expr', 10*ib] vecxzX vecxzY 0);               # PDelta transformation; vecxz 'is', parallel 'to', the 'length', of 'the', cap 'beam.'
       model.element('nonlinearBeamColumn',, (1000*ib)+10) '[expr', 100*ib+1] '[expr', 100*ib+2]       np  , (1000*ib+10) , 10*ib #(Hinge)
       rigidLink 'beam',, (100*ib+2) '[expr', 100*ib+3];                            #(Rigid)
       model.element('nonlinearBeamColumn',, (1000*ib)+20) '[expr', 100*ib+7] '[expr', 100*ib+6]       np  , (1000*ib+20) , 10*ib #(Hinge)
       rigidLink 'beam',, (100*ib+6) '[expr', 100*ib+5];                            #(Rigid)

# Bent 12
Dcol = 66.0*in;                                                                      # Width, of, octagonal, column (to, flat, sides)
RcolDiag = Dcol*sqrt(4+2*sqrt(2))/(2+2*sqrt(2));       # Radius, of, octagonal, column (to, corners)
vecxzX = [lindex, vecxzXList, 10];                                           # X, component, of, vecxz, vector, for, local, axis, transformation
vecxzY = [lindex, vecxzYList, 10];                                           # Y, component, of, vecxz, vector, for, local, axis, transformation
model.geomTransf(TransfType, 120, vecxzX, vecxzY, 0);                                    # PDelta, transformation; vecxz, 'is', parallel, 'to', the, 'length', of, 'the', cap, 'beam.'
model.element('nonlinearBeamColumn', 12010, 1201, 1202, np, 12010, 120); #(Hinge)
rigidLink, 'beam', 1202, 1203;                                                                       #(Rigid)
model.element('nonlinearBeamColumn', 12020, 1209, 1208, np, 12020, 120); #(Hinge)
rigidLink, 'beam', 1208, 1207;                                                                       #(Rigid)
model.element('nonlinearBeamColumn', 12030, 1211, 1210, np, 12030, 120); #(Hinge)
rigidLink, 'beam', 1210, 1205;                                                                       #(Rigid)
# Bent 13 NE
Dcol = 48.0*in;                                                                      # Width, of, octagonal, column (to, flat, sides)
RcolDiag = Dcol*sqrt(4+2*sqrt(2))/(2+2*sqrt(2));       # Radius, of, octagonal, column (to, corners)
vecxzX = [lindex, vecxzXList, 11];                                           # X, component, of, vecxz, vector, for, local, axis, transformation
vecxzY = [lindex, vecxzYList, 11];                                           # Y, component, of, vecxz, vector, for, local, axis, transformation
model.geomTransf(TransfType, 130, vecxzY, vecxzX, 0);                                    # PDelta, transformation; vecxz, 'is', PERPENDICULAR, 'to', the, 'length', of, 'the', cap, 'beam.', (local, y parallel)
model.element('nonlinearBeamColumn', 13010, 1301, 1302, np, 13010, 130); #(Hinge)
rigidLink, 'beam', 1302, 1303;                                                                       #(Rigid)
model.element('nonlinearBeamColumn', 13020, 1307, 1306, np, 13020, 130); #(Hinge)
rigidLink, 'beam', 1306, 1305;                                                                       #(Rigid)
# Bent 13 NR (WIDE (INTERLOCKING) SECTION)
model.element('nonlinearBeamColumn', 13040, 1313, 1314, np, 13040, 130); #(Hinge)
rigidLink, 'beam', 1314, 1315;                                                                       #(Rigid)
# Bent 14 NE
vecxzX = [lindex, vecxzXList, 12];                                           # X, component, of, vecxz, vector, for, local, axis, transformation
vecxzY = [lindex, vecxzYList, 12];                                           # Y, component, of, vecxz, vector, for, local, axis, transformation
model.geomTransf(TransfType, 140, vecxzY, vecxzX, 0);                                    # PDelta, transformation; vecxz, 'is', PERPENDICULAR, 'to', the, 'length', of, 'the', cap, 'beam.', (local, y parallel)
model.element('nonlinearBeamColumn', 14010, 1401, 1402, np, 14010, 140); #(Hinge)
rigidLink, 'beam', 1402, 1403;                                                                       #(Rigid)
model.element('nonlinearBeamColumn', 14020, 1407, 1406, np, 14020, 140); #(Hinge)
rigidLink, 'beam', 1406, 1405;                                                                       #(Rigid)
model.element('nonlinearBeamColumn', 14030, 1409, 1408, np, 14030, 140); #(Hinge)
rigidLink, 'beam', 1408, 1404;                                                                       #(Rigid)
# Bent 14 NR (WIDE (INTERLOCKING) SECTION)
model.element('nonlinearBeamColumn', 14040, 1413, 1414, np, 14040, 140); #(Hinge)
rigidLink, 'beam', 1414, 1415;                                                                       #(Rigid)
# CAP BEAM ELEMENTS
# (Numbering scheme: 10000*Bent#+10,20,30,40)
#source ReadMPR.tcl;                                                                       # Set up ReadMPR procedure for obtaining cap beam section properties
CSDir = "./Dimensions/CapCS/";                                                   # Directory containing cap beam cross section information
CSType = "Cap";                                                                              # Cross section type is cap beam
# Read in the concrete strength for each cap beam:
fpfceCap = [open "./Dimensions/fceCap.txt" r];
fceCapList = [read, fpfceCap];
close, fpfceCap
# Bents 2-10
for {ib = 2} {ib <= 9} {ib += 1} {
       lassign '[ReadMPR', CSDir CSType ib {}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       fceCap = [lindex fceCapList '[expr', ib-2]]
       EcCap = [expr 57000.0*sqrt(fceCap*ksi_psi)/ksi_psi]
       GcCap = [expr EcCap/(2.0*(1.0+Uc))]
       model.element('elasticBeamColumn',, (10000*ib)+10) '[expr', 100*ib+3] '[expr', 100*ib+4] A EcCap GcCap J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       model.element('elasticBeamColumn',, (10000*ib)+20) '[expr', 100*ib+4] '[expr', 100*ib+5] A EcCap GcCap J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);

# Bents 10-11
for {ib = 10} {ib <= 11} {ib += 1} {
       lassign '[ReadMPR', CSDir CSType ib {}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       fceCap = [lindex fceCapList '[expr', ib-2]]
       EcCap = [expr 57000.0*sqrt(fceCap*ksi_psi)/ksi_psi]
       GcCap = [expr EcCap/(2.0*(1.0+Uc))]
       model.element('elasticBeamColumn',, (10000*ib)+10) '[expr', 100*ib+3] '[expr', 100*ib+4] A EcCap GcCap J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       model.element('elasticBeamColumn',, (10000*ib)+20) '[expr', 100*ib+4] '[expr', 100*ib+5] A EcCap GcCap J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);

# Bent 12
lassign, '[ReadMPR', CSDir, CSType, 12 {}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 120010, 1203, 1204, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 120020, 1204, 1205, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 120030, 1205, 1206, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 120040, 1206, 1207, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
# Bent 13 NE
lassign, '[ReadMPR', CSDir, CSType, 13 {}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 130010, 1303, 1304, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 130020, 1304, 1305, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
# Bent 14 NE
lassign, '[ReadMPR', CSDir, CSType, 14 {}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 140010, 1403, 1404, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 140020, 1404, 1405, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
# DECK ELEMENTS
# (Numbering scheme: 10*1stBent#+1stIntermediateNode#)
#source ReadMPR.tcl;                                                                       # Set up ReadMPR procedure for obtaining deck section properties
CSDir = "./Dimensions/DeckCS/";                                                   # Directory containing deck cross section information
CSType = "Deck";                                                                              # Cross section type is deck
# Abut 1 to Bent 2
lassign, '[ReadMPR', CSDir, CSType, 1 {-NodeNum, 0}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 10, 1010, 10001, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 11, 10001, 10002, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 12, 10002, 10003, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 13, 10003, 10004, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 14, 10004, 204, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
# Bent 2 to Bent 5
for {ib = 2} {ib <= 4} {ib += 1} {
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 0}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib) '[expr', 100*ib)+4] '[expr', 10000*ib+1] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 1}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib)+1) '[expr', 10000*ib+1] '[expr', 10000*ib+2] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 2}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib)+2) '[expr', 10000*ib+2] '[expr', 10000*ib+3] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 3}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib)+3) '[expr', 10000*ib+3] '[expr', 10000*ib+4] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 4}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib)+4) '[expr', 10000*ib+4] '[expr', 100*(ib+1)+4] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);

# Bent 5 to Bent 6
lassign, '[ReadMPR', CSDir, CSType, 5 {-NodeNum, 0}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 50, 504, 50001, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 5 {-NodeNum, 1}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 51, 500010, 50002, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 5 {-NodeNum, 2}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 52, 50002, 50003, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 5 {-NodeNum, 3}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 53, 50003, 50004, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 5 {-NodeNum, 4}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 54, 50004, 604, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
# Bent 6 to Bent 8
for {ib = 6} {ib <= 7} {ib += 1} {
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 0}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib) '[expr', 100*ib)+4] '[expr', 10000*ib+1] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 1}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib)+1) '[expr', 10000*ib+1] '[expr', 10000*ib+2] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 2}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib)+2) '[expr', 10000*ib+2] '[expr', 10000*ib+3] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 3}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib)+3) '[expr', 10000*ib+3] '[expr', 10000*ib+4] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 4}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib)+4) '[expr', 10000*ib+4] '[expr', 100*(ib+1)+4] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);

# Bent 8 to Bent 9
lassign, '[ReadMPR', CSDir, CSType, 8 {-NodeNum, 0}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 80, 804, 80001, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 8 {-NodeNum, 1}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 81, 80001, 80002, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 8 {-NodeNum, 2}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 82, 80002, 80003, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 8 {-NodeNum, 3}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 83, 80003, 80004, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 8 {-NodeNum, 4}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 84, 800040, 904, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
# Bent 9 to Bent 11
for {ib = 9} {ib <= 10} {ib += 1} {
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 0}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib) '[expr', 100*ib)+4] '[expr', 10000*ib+1] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 1}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib)+1) '[expr', 10000*ib+1] '[expr', 10000*ib+2] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 2}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib)+2) '[expr', 10000*ib+2] '[expr', 10000*ib+3] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 3}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib)+3) '[expr', 10000*ib+3] '[expr', 10000*ib+4] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);
       lassign '[ReadMPR', CSDir CSType ib {-NodeNum 4}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       model.element('elasticBeamColumn',, (10*ib)+4) '[expr', 10000*ib+4] '[expr', 100*(ib+1)+4] A Ec Gc J Iy Iz 301 -model.mass('[expr', mconc*A] -cMass);

# Bent 11 to Bent 12
lassign, '[ReadMPR', CSDir, CSType, 11 {-NodeNum, 0}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 110, 1104, 110001, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 11 {-NodeNum, 1}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 111, 110001, 110002, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 11 {-NodeNum, 2}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 112, 110002, 110003, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 11 {-NodeNum, 3}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 113, 110003, 110004, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 11 {-NodeNum, 4}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 114, 110004, 1205, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
# Bent 12 to Bent 13 NE
lassign, '[ReadMPR', CSDir, CSType, 12 {-NodeNum, 0}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 120, 1204, 120001, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 12 {-NodeNum, 1}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 121, 1200010, 120002, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 12 {-NodeNum, 2}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 122, 120002, 120003, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 12 {-NodeNum, 3}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 123, 120003, 120004, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 12 {-NodeNum, 4}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 124, 120004, 1304, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
# Bent 12 to Bent 13 NR
lassign, '[ReadMPR', CSDir, CSType, 12 {-NodeNum, 5}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 125, 1206, 120005, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 12 {-NodeNum, 6}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 126, 1200050, 120006, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 12 {-NodeNum, 7}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 127, 120006, 120007, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 12 {-NodeNum, 8}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 128, 120007, 120008, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 12 {-NodeNum, 9}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 129, 120008, 1315, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
# Bent 13 NE to Bent 14 NE
lassign, '[ReadMPR', CSDir, CSType, 13 {-NodeNum, 0}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 130, 1304, 130001, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 13 {-NodeNum, 1}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 131, 130001, 130002, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 13 {-NodeNum, 2}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 132, 130002, 130003, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 13 {-NodeNum, 3}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 133, 130003, 130004, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 13 {-NodeNum, 4}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 134, 130004, 1404, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
# Bent 13 NR to Bent 14 NR
lassign, '[ReadMPR', CSDir, CSType, 13 {-NodeNum, 5}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 135, 1315, 130005, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 13 {-NodeNum, 6}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 136, 130005, 130006, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 13 {-NodeNum, 7}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 137, 130006, 130007, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 13 {-NodeNum, 8}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 138, 130007, 130008, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 139, 130008, 1415, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
# Bent 14 NE to Abut 15 NE
lassign, '[ReadMPR', CSDir, CSType, 14 {-NodeNum, 0}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 140, 1404, 140001, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 14 {-NodeNum, 1}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 141, 140001, 140002, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 14 {-NodeNum, 2}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 142, 140002, 140003, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 14 {-NodeNum, 3}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 143, 140003, 140004, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
lassign, '[ReadMPR', CSDir, CSType, 14 {-NodeNum, 4}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 144, 140004, 15010, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
# Bent 14 NR to Abut 15 NR
lassign, '[ReadMPR', CSDir, CSType, 14 {-NodeNum, 5}], A 'Iy', Iz, J; # Cross, 'Section', 'properties'
model.element('elasticBeamColumn', 145, 1415, 140005, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 146, 140005, 140006, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 147, 140006, 140007, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 148, 140007, 140008, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);
model.element('elasticBeamColumn', 149, 140008, 15040, A Ec, Gc, J Iy, Iz, 301, mass='[expr', mconc*A] -cMass);

