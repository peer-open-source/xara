from math import cos,sin,sqrt,pi
import opensees as ops
#
# GRAVITY (STATIC & MODAL) ANALYSIS -------------------------------------------------------------------------------
#
# COLUMN/ABUTMENT SUPPORT FORCE RECORDERS
gravRec = ""
gravRec.append([ana.recorder('Node', -xml, dataDir)/nodeforceZ.out, time=True, -model.node(1020, 1010, 1030, 201, 207, 301, 307, 401, 407, 501, 507, 601, 607, 701, 707, 801, 807, 901, 907, 1001, 1007, 1101, 1107, 1201, 1209, 1211, 1301, 1307, 1313, 1401, 1407, 1409, 1413, 15020, 15030, 15050, 15060 -dof, 3 reaction]))
gravRec.append([ana.recorder('Node', -xml, dataDir)/nodeforceX.out, time=True, -model.node(1020, 1010, 1030, 201, 207, 301, 307, 401, 407, 501, 507, 601, 607, 701, 707, 801, 807, 901, 907, 1001, 1007, 1101, 1107, 1201, 1209, 1211, 1301, 1307, 1313, 1401, 1407, 1409, 1413, 15020, 15030, 15050, 15060 -dof, 1 reaction]))
gravRec.append([ana.recorder('Node', -xml, dataDir)/nodeforceY.out, time=True, -model.node(1020, 1010, 1030, 201, 207, 301, 307, 401, 407, 501, 507, 601, 607, 701, 707, 801, 807, 901, 907, 1001, 1007, 1101, 1107, 1201, 1209, 1211, 1301, 1307, 1313, 1401, 1407, 1409, 1413, 15020, 15030, 15050, 15060 -dof, 2 reaction]))
gravRec.append([ana.recorder('Node', -xml, dataDir)/nodeMoment.out, time=True, -model.node(1020, 1010, 1030, 201, 207, 301, 307, 401, 407, 501, 507, 601, 607, 701, 707, 801, 807, 901, 907, 1001, 1007, 1101, 1107, 1201, 1209, 1211, 1301, 1307, 1313, 1401, 1407, 1409, 1413, 15020, 15030, 15050, 15060 -dof, 4 5, reaction]))
# COLUMN POTENTIAL PIN X/Y MOMENT RECORDER
gravRec.append([ana.recorder('Node', -xml, dataDir)/nodemomXY.out, time=True, -model.node(201, 207, 302, 306, 402, 406, 502, 506, 702, 706, 802, 806, 902, 906, 1002, 1006, 1102, 1106, 1201, 1211, 1209, 1301, 1307, 1313, 1401, 1407, 1413 -dof, 4 5, reaction]))
# TOP OF COLUMN DISPLACEMENT AND ACCELERATION RECORDERS
gravRec.append([ana.recorder('Node', -file, dataDir)/nodeDisp.out, time=True, -model.node(1408, 1302 -dof, 1 2, 3 disp]))
gravRec.append([ana.recorder('Node', -file, dataDir)/nodeAcc.out, time=True, -model.node(1408, 1302 -dof, 1 2, 3 accel]))
# MODESHAPES AND PERIODS
file mkdir "dataDir/ModeShape"; # Output folder
nPds = 8; # Number, of, periods, to, analyze
# Record Undeformed Original Shape
source, 'PrintNodes.tcl'
# Establish Recorders for Mode Shapes
for {k = 1} {k <= nPds} {k += 1} {
       gravRec.append([recorder Node -file dataDir/ModeShape/modek.txt -time -node 1010 1020 1030 201 202 203 204 205 206 207 301 302 303 304 305 306 307 401 402 403 404 405 406 407 501 502 503 504 505 506 507 601 602 603 604 605 606 607 701 702 703 704 705 706 707 801 802 803 804 805 806 807 901 902 903 904 905 906 907 1001 1002 1003 1004 1005 1006 1007 1101 1102 1103 1104 1105 1106 1107 1201 1202 1203 1204 1205 1206 1207 1208 1209 1210 1211 1301 1302 1303 1304 1305 1306 1307 1313 1314 1315 1401 1402 1403 1404 1405 1406 1407 1408 1409 1413 1414 1415 15010 15020 15030 15040 15050 15060 10001 10002 10003 10004 20001 20002 20003 20004 30001 30002 30003 30004 40001 40002 40003 40004 50001 50002 50003 50004 60001 60002 60003 60004 70001 70002 70003 70004 80001 80002 80003 80004 90001 90002 90003 90004 100001 100002 100003 100004 110001 110002 110003 110004 120001 120002 120003 120004 130001 130002 130003 130004 140001 140002 140003 140004 120005 120006 120007 120008 130005 130006 130007 130008 140005 140006 140007 140008 500010 800040 1200010 1200050 -dof 1 2 3 "eigen k."])

gravRec.append([ana.recorder('Node', -file, dataDir)/allNodesGrav.disp       -precision, 5       -time       -dof, 1 2, 3 4, 5 6, disp])
# Period before gravity
wa = [eigen, nPds];

# source LibIO.tcl
# brace2::io::write_modes "dataDir/ModeShape/modes.yaml" nPds

Periods =        [open, dataDir/ModeShape/PeriodsPreG.txt, w];
print("Fundamental-Period Before Gravity Analysis:")
for {iPd = 1} {iPd <= nPds} {iPd += 1} {
       wwa = [lindex wa iPd-1];
       Ta =, 2*pi/sqrt(wwa)
       print("PeriodiPd= Ta")
       puts Periods "Ta";

close, Periods;
# GRAVITY LOAD
source, 'getVertCosines.tcl'
vecxz301 = "0 0 1";
#puts [getVertCosines 20010 "0 0 1"]
pattern, 'Plain', 999, Linear {
       # Cap Beam Element Loads
       CSDir = "./Dimensions/CapCS/";                                                   # Directory containing cap beam cross section information
       CSType = "Cap";                                                                              # Cross section type is cap beam
       # Bents 2-11
       for {ib = 2} {ib <= 11} {ib += 1} {
              lassign '[ReadMPR', CSDir CSType ib {}] A 'Iy', Iz J; # Cross 'Section', 'properties'
              lassign '[getVertCosines',, (10000*ib+10) vecxz301] 'cosy', cosz 'cosx', curEle;
              eleLoad ele=curEle,  type=True, beamUniform='[expr', -wconc*A*cosy] '[expr', -wconc*A*cosz] '[expr', -wconc*A*cosx];
              lassign '[getVertCosines',, (10000*ib+20) vecxz301] 'cosy', cosz 'cosx', curEle;
              eleLoad ele=curEle,  type=True, beamUniform='[expr', -wconc*A*cosy] '[expr', -wconc*A*cosz] '[expr', -wconc*A*cosx];

       # Bent 12
       lassign '[ReadMPR', CSDir CSType 12 {}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       for {eleCtr = 0; eleNum = "120010 120020 120030 120040";} {eleCtr<=3} {incr eleCtr} {
              lassign '[getVertCosines', [lindex eleNum eleCtr] vecxz301] 'cosy', cosz 'cosx', curEle;
              eleLoad ele=curEle,  type=True, beamUniform='[expr', -wconc*A*cosy] '[expr', -wconc*A*cosz] '[expr', -wconc*A*cosx];

       # Bent 13 NE
       lassign '[ReadMPR', CSDir CSType 13 {}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       lassign '[getVertCosines', 130010 vecxz301] 'cosy', cosz 'cosx', curEle;
       eleLoad ele=curEle,  type=True, beamUniform='[expr', -wconc*A*cosy] '[expr', -wconc*A*cosz] '[expr', -wconc*A*cosx];
       lassign '[getVertCosines', 130020 vecxz301] 'cosy', cosz 'cosx', curEle;
       eleLoad ele=curEle,  type=True, beamUniform='[expr', -wconc*A*cosy] '[expr', -wconc*A*cosz] '[expr', -wconc*A*cosx];
       # Bent 14 NE
       lassign '[ReadMPR', CSDir CSType 14 {}] A 'Iy', Iz J; # Cross 'Section', 'properties'
       lassign '[getVertCosines', 140010 vecxz301] 'cosy', cosz 'cosx', curEle;
       eleLoad ele=curEle,  type=True, beamUniform='[expr', -wconc*A*cosy] '[expr', -wconc*A*cosz] '[expr', -wconc*A*cosx];
       lassign '[getVertCosines', 140020 vecxz301] 'cosy', cosz 'cosx', curEle;
       eleLoad ele=curEle,  type=True, beamUniform='[expr', -wconc*A*cosy] '[expr', -wconc*A*cosz] '[expr', -wconc*A*cosx];
       # Deck Element Loads
       CSDir = "./Dimensions/DeckCS/";                                                   # Directory containing deck cross section information
       CSType = "Deck";                                                                              # Cross section type is deck
       # NE Side
       for {ib = 1} {ib <= 14} {ib += 1} {
              for {ibb = 0} {ibb <= 4} {ibb += 1} {
                     lassign [ReadMPR CSDir CSType ib "-NodeNum ibb"] A Iy Iz J; # Cross Section properties
                     lassign '[getVertCosines',, (10*ib+ibb) vecxz301] 'cosy', cosz 'cosx', curEle;
                     eleLoad ele=curEle,  type=True, beamUniform='[expr', -wconc*A*cosy] '[expr', -wconc*A*cosz] '[expr', -wconc*A*cosx];
                     #eleLoad -ele, (10*ib+ibb) -type -beamUniform 0, -wconc*A # Uniformly distributed element load acting in vertical (local y) direction of element


       for {ib = 12} {ib <= 14} {ib += 1} {
              for {ibb = 5} {ibb <= 9} {ibb += 1} {
                     lassign [ReadMPR CSDir CSType ib "-NodeNum ibb"] A Iy Iz J; # Cross Section properties
                     lassign '[getVertCosines',, (10*ib+ibb) vecxz301] 'cosy', cosz 'cosx', curEle;
                     eleLoad ele=curEle,  type=True, beamUniform='[expr', -wconc*A*cosy] '[expr', -wconc*A*cosz] '[expr', -wconc*A*cosx];
                     #eleLoad -ele, (10*ib+ibb) -type -beamUniform 0, -wconc*A # Uniformly distributed element load acting in vertical (local y) direction of element



   
# Perform static analysis
wipeAnalysis
ana.test('NormDispIncr', 1.0e-8, 10, 0);       
ana.algorithm(Newton);       
ana.integrator('LoadControl', 0.1);
ana.numberer(Plain);
ana.constraints(Transformation);
ana.system(SparseGeneral);
ana.analysis(Static);
print, JSON=True, file=dataDir, /modelDetails.json
ana.analyze(10);
ana.loadConst(time=0.0, );
# Period after gravity
wb = [eigen, fullGenLapack=nPds], ;
Periods =        [open, dataDir/ModeShape/PeriodsPostG.txt, w];
print("Fundamental-Period After Gravity Analysis:")
for {iPd = 1} {iPd <= nPds} {iPd += 1} {
       wwb = [lindex wb iPd-1];
       Tb =, 2*pi/sqrt(wwb)
       print("PeriodiPd= Tb")
       puts Periods "Tb";

close, Periods;
for {ctrRec = 0} {ctrRec<[llength, gravRec]} {incr, ctrRec} {
       ana.remove('recorder', [lindex gravRec ctrRec])


