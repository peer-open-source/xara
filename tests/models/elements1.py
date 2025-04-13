package require brace2
#
# ELEMENTS -------------------------------------------------------------------------------
#
# Transformation and P-Delta
set TransfType PDelta;
geomTransf $TransfType 301 0 0 1;                # vecxz is in the +Z direction (deck and cap beams: local z upwards on section, local y horiz left on section facing j)      
# COLUMN ELEMENTS
# (Numbering scheme: Bottom Length, Hinge: 1000*Bent#+10*Column#;  Top Length, Rigid: 1000*Bent#+10*Column#+1)
# (Column#: 1=Left, 2=Right, 3=Center, 4=Single)
# Column section tag (ColSecTag) follows column bent numbering scheme.
# Read in the local axes transformation vectors (vecxz) for the columns:
set fpvecxzX [open "./Dimensions/vecxzXcol.txt" r];
set vecxzXList [read $fpvecxzX];
close $fpvecxzX
set fpvecxzY   [open "./Dimensions/vecxzYcol.txt" r];
set vecxzYList [read $fpvecxzY];
close $fpvecxzY

set np 4;                                        # number of Gauss integration points for nonlinear curvature distribution-- np=2 for linear distribution ok
# Bent 2
set Dcol [expr 84.0*$in];                                        # Width of octagonal column (to flat sides)
set RcolDiag    [expr $Dcol*sqrt(4+2*sqrt(2))/(2+2*sqrt(2))];    # Radius of octagonal column (to corners)
set vecxzX      [lindex $vecxzXList [expr 0]];                 # X component of vecxz vector for local axis transformation
set vecxzY      [lindex $vecxzYList [expr 0]];                 # Y component of vecxz vector for local axis transformation
geomTransf $TransfType 20 $vecxzX $vecxzY 0;         # PDelta transformation; vecxz is parallel to the length of the cap beam.
element nonlinearBeamColumn 2010   201   202     $np   2010  20; #(Hinge)
rigidLink beam 202 203;                                             #(Rigid)
element nonlinearBeamColumn 2020   207   206      $np   2020  20; #(Hinge)
rigidLink beam 206 205;                                             #(Rigid)
# Bents 3-11
for {set ib 3} {ib <= 11} {incr ib 1} {
    set Dcol [expr 84.0*$in];                                        # Width of octagonal column (to flat sides)
    set RcolDiag    [expr $Dcol*sqrt(4+2*sqrt(2))/(2+2*sqrt(2))];    # Radius of octagonal column (to corners)
    set vecxzX         [lindex $vecxzXList [expr ib-2]];                 # X component of vecxz vector for local axis transformation
    set vecxzY         [lindex $vecxzYList [expr ib-2]];                 # Y component of vecxz vector for local axis transformation
    geomTransf $TransfType [expr 10*ib] $vecxzX $vecxzY 0;         # PDelta transformation; vecxz is parallel to the length of the cap beam.
    element nonlinearBeamColumn [expr 1000*ib+10] [expr 100*ib+1] [expr 100*ib+2]    $np   [expr 1000*ib+10]  [expr 10*ib]; #(Hinge)
    rigidLink beam [expr 100*ib+2] [expr 100*ib+3];                #(Rigid)
    element nonlinearBeamColumn [expr 1000*ib+20] [expr 100*ib+7] [expr 100*ib+6]    $np   [expr 1000*ib+20]  [expr 10*ib]; #(Hinge)
    rigidLink beam [expr 100*ib+6] [expr 100*ib+5];                #(Rigid)
}
# Bent 12
set Dcol [expr 66.0*$in];                                        # Width of octagonal column (to flat sides)
set RcolDiag    [expr $Dcol*sqrt(4+2*sqrt(2))/(2+2*sqrt(2))];    # Radius of octagonal column (to corners)
set vecxzX         [lindex $vecxzXList 10];                         # X component of vecxz vector for local axis transformation
set vecxzY         [lindex $vecxzYList 10];                         # Y component of vecxz vector for local axis transformation
geomTransf $TransfType 120 $vecxzX $vecxzY 0;                     # PDelta transformation; vecxz is parallel to the length of the cap beam.
element nonlinearBeamColumn 12010   1201   1202     $np   12010  120; #(Hinge)
rigidLink beam 1202 1203;                                         #(Rigid)
element nonlinearBeamColumn 12020   1209   1208     $np   12020  120; #(Hinge)
rigidLink beam 1208 1207;                                         #(Rigid)
element nonlinearBeamColumn 12030   1211   1210     $np   12030  120; #(Hinge)
rigidLink beam 1210 1205;                                         #(Rigid)
# Bent 13 NE
set Dcol [expr 48.0*$in];                                        # Width of octagonal column (to flat sides)
set RcolDiag    [expr $Dcol*sqrt(4+2*sqrt(2))/(2+2*sqrt(2))];    # Radius of octagonal column (to corners)
set vecxzX         [lindex $vecxzXList 11];                         # X component of vecxz vector for local axis transformation
set vecxzY         [lindex $vecxzYList 11];                         # Y component of vecxz vector for local axis transformation
geomTransf $TransfType 130 $vecxzY $vecxzX 0;                     # PDelta transformation; vecxz is PERPENDICULAR to the length of the cap beam. (local y parallel)
element nonlinearBeamColumn 13010   1301   1302     $np   13010  130; #(Hinge)
rigidLink beam 1302 1303;                                         #(Rigid)
element nonlinearBeamColumn 13020   1307   1306     $np   13020  130; #(Hinge)
rigidLink beam 1306 1305;                                         #(Rigid)
# Bent 13 NR (WIDE (INTERLOCKING) SECTION)
element nonlinearBeamColumn 13040   1313   1314     $np   13040  130; #(Hinge)
rigidLink beam 1314 1315;                                         #(Rigid)
# Bent 14 NE
set vecxzX         [lindex $vecxzXList 12];                         # X component of vecxz vector for local axis transformation
set vecxzY         [lindex $vecxzYList 12];                         # Y component of vecxz vector for local axis transformation
geomTransf $TransfType 140 $vecxzY $vecxzX 0;                     # PDelta transformation; vecxz is PERPENDICULAR to the length of the cap beam. (local y parallel)
element nonlinearBeamColumn 14010   1401   1402     $np   14010  140; #(Hinge)
rigidLink beam 1402 1403;                                         #(Rigid)
element nonlinearBeamColumn 14020   1407   1406     $np   14020  140; #(Hinge)
rigidLink beam 1406 1405;                                         #(Rigid)
element nonlinearBeamColumn 14030   1409   1408     $np   14030  140; #(Hinge)
rigidLink beam 1408 1404;                                         #(Rigid)
# Bent 14 NR (WIDE (INTERLOCKING) SECTION)
element nonlinearBeamColumn 14040   1413   1414     $np   14040  140; #(Hinge)
rigidLink beam 1414 1415;                                         #(Rigid)
# CAP BEAM ELEMENTS
# (Numbering scheme: 10000*Bent#+10,20,30,40)
#source ReadMPR.tcl;                                         # Set up ReadMPR procedure for obtaining cap beam section properties
set CSDir "./Dimensions/CapCS/";                              # Directory containing cap beam cross section information
set CSType "Cap";                                             # Cross section type is cap beam
# Read in the concrete strength for each cap beam:
set fpfceCap [open "./Dimensions/fceCap.txt" r];
set fceCapList [read $fpfceCap];
close $fpfceCap
# Bents 2-10
for {set ib 2} {ib <= 9} {incr ib 1} {
    lassign [ReadMPR $CSDir $CSType ib {}] A Iy Iz J; # Cross Section properties
    set fceCap [lindex $fceCapList [expr ib-2]]
    set EcCap [expr 57000.0*sqrt($fceCap*$ksi_psi)/$ksi_psi]
    set GcCap [expr $EcCap/(2.0*(1.0+$Uc))]
    element elasticBeamColumn [expr 10000*ib+10] [expr 100*ib+3] [expr 100*ib+4] $A $EcCap $GcCap $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
    element elasticBeamColumn [expr 10000*ib+20] [expr 100*ib+4] [expr 100*ib+5] $A $EcCap $GcCap $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
}
# Bents 10-11
for {set ib 10} {ib <= 11} {incr ib 1} {
    lassign [ReadMPR $CSDir $CSType ib {}] A Iy Iz J; # Cross Section properties
    set fceCap [lindex $fceCapList [expr ib-2]]
    set EcCap [expr 57000.0*sqrt($fceCap*$ksi_psi)/$ksi_psi]
    set GcCap [expr $EcCap/(2.0*(1.0+$Uc))]
    element elasticBeamColumn [expr 10000*ib+10] [expr 100*ib+3] [expr 100*ib+4] $A $EcCap $GcCap $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
    element elasticBeamColumn [expr 10000*ib+20] [expr 100*ib+4] [expr 100*ib+5] $A $EcCap $GcCap $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
}
# Bent 12
lassign [ReadMPR $CSDir $CSType 12 {}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 120010 1203 1204 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 120020 1204 1205 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 120030 1205 1206 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 120040 1206 1207 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# Bent 13 NE
lassign [ReadMPR $CSDir $CSType 13 {}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 130010 1303 1304 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 130020 1304 1305 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# Bent 14 NE
lassign [ReadMPR $CSDir $CSType 14 {}] A Iy Iz J; # Cross Section properties
element elasticBeamColumn 140010 1403 1404 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 140020 1404 1405 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;

# DECK ELEMENTS
# (Numbering scheme: 10*1stBent#+1stIntermediateNode#)
#source ReadMPR.tcl;                                         # Set up ReadMPR procedure for obtaining deck section properties
set CSDir "./Dimensions/DeckCS/";                              # Directory containing deck cross section information
set CSType "Deck";                                             # Cross section type is deck
# Abut 1 to Bent 2
10:   (1010, 10001, (elasticBeamColumn, 1, 0)),
element elasticBeamColumn 11 10001 10002 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 12 10002 10003 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 13 10003 10004 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 14 10004 204 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;

# Bent 2 to Bent 5
for {set ib 2} {ib <= 4} {incr ib 1} {
    10*ib: (100*ib+4, 10000*ib+1, (elasticBeamColumn, ib, 0)),
    10*ib+1: (10000*ib+1, 10000*ib+2, (elasticBeamColumn, ib, 1)),
    10*ib+2: (10000*ib+2, 10000*ib+3, (elasticBeamColumn, ib, 2)),
    10*ib+3: (10000*ib+3, 10000*ib+4, (elasticBeamColumn, ib, 3)),
    10*ib+4: (10000*ib+4, 100*(ib+1)+4, (elasticBeamColumn, ib, 4)),
}
# Bent 5 to Bent 6
50:   (504, 50001, (elasticBeamColumn, 5, 0)),
51: (50001, 50002, (elasticBeamColumn, 5, 1)),
52: (50002, 50003, (elasticBeamColumn, 5, 2)),
53: (50003, 50004, (elasticBeamColumn, 5, 3)),
54: (50004,   604, (elasticBeamColumn, 5, 4)),
# Bent 6 to Bent 8
for {set ib 6} {ib <= 7} {incr ib 1} {
    10*ib:     (100*ib+4,   10000*ib+1, (elasticBeamColumn, ib, 0)),
    10*ib+1: (10000*ib+1,   10000*ib+2, (elasticBeamColumn, ib, 1)),
    10*ib+2: (10000*ib+2,   10000*ib+3, (elasticBeamColumn, ib, 2)),
    10*ib+3: (10000*ib+3,   10000*ib+4, (elasticBeamColumn, ib, 3)),
    10*ib+4: (10000*ib+4, 100*(ib+1)+4, (elasticBeamColumn, ib, 4)),
}
# Bent 8 to Bent 9
80:   (804, 80001, (elasticBeamColumn, 8, 0)),
81: (80001, 80002, (elasticBeamColumn, 8, 1)),
82: (80002, 80003, (elasticBeamColumn, 8, 2)),
83: (80003, 80004, (elasticBeamColumn, 8, 3)),
84: (80004,   904, (elasticBeamColumn, 8, 4)),

# Bent 9 to Bent 11
for {set ib 9} {ib <= 10} {incr ib 1} {
    10*ib:     (100*ib+4,   10000*ib+1, (elasticBeamColumn, ib, 0)),
    10*ib+1: (10000*ib+1,   10000*ib+2, (elasticBeamColumn, ib, 1)),
    10*ib+2: (10000*ib+2,   10000*ib+3, (elasticBeamColumn, ib, 2)),
    10*ib+3: (10000*ib+3,   10000*ib+4, (elasticBeamColumn, ib, 3)),
    10*ib+4: (10000*ib+4, 100*(ib+1)+4, (elasticBeamColumn, ib, 4)),
}
# Bent 11 to Bent 12
110:   (1104, 110001, (elasticBeamColumn, 11, 0)),
111: (110001, 110002, (elasticBeamColumn, 11, 1)),
112: (110002, 110003, (elasticBeamColumn, 11, 2)),
113: (110003, 110004, (elasticBeamColumn, 11, 3)),
114: (110004,   1205, (elasticBeamColumn, 11, 4)),
# Bent 12 to Bent 13 NE
120:   (1204, 120001, (elasticBeamColumn, 12, 0)),
121: (120001, 120002, (elasticBeamColumn, 12, 1)),
122: (120002, 120003, (elasticBeamColumn, 12, 2)),
123: (120003, 120004, (elasticBeamColumn, 12, 3)),
124: (120004,   1304, (elasticBeamColumn, 12, 4)),
# Bent 12 to Bent 13 NR
125:   (1206, 120005, (elasticBeamColumn, 12, 5)),
126: (120005, 120006, (elasticBeamColumn, 12, 6)),
127: (120006, 120007, (elasticBeamColumn, 12, 7)),
128: (120007, 120008, (elasticBeamColumn, 12, 8)),
129: (120008,   1315, (elasticBeamColumn, 12, 9)),
# Bent 13 NE to Bent 14 NE
130:   (1304, 130001, (elasticBeamColumn, 13, 0)),
131: (130001, 130002, (elasticBeamColumn, 13, 1)),
132: (130002, 130003, (elasticBeamColumn, 13, 2)),
133: (130003, 130004, (elasticBeamColumn, 13, 3)),
134: (130004,   1404, (elasticBeamColumn, 13, 4)),
# Bent 13 NR to Bent 14 NR
135:   (1315, 130005, (elasticBeamColumn, 13, 5)),
136: (130005, 130006, (elasticBeamColumn, 13, 6)),
137: (130006, 130007, (elasticBeamColumn, 13, 7)),
138: (130007, 130008, (elasticBeamColumn, 13, 8)),
element elasticBeamColumn 139 130008 1415 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# Bent 14 NE to Abut 15 NE
140:   (1404, 140001, (elasticBeamColumn, 14, 0)),
141: (140001, 140002, (elasticBeamColumn, 14, 1)),
142: (140002, 140003, (elasticBeamColumn, 14, 2)),
143: (140003, 140004, (elasticBeamColumn, 14, 3)),
144: (140004,  15010, (elasticBeamColumn, 14, 4)),
# Bent 14 NR to Abut 15 NR
145:   (1415, 140005, (elasticBeamColumn, 14, 5)),
element elasticBeamColumn 146 140005 140006 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 147 140006 140007 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 148 140007 140008 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
element elasticBeamColumn 149 140008 15040 $A $Ec $Gc $J $Iy $Iz 301 -mass [expr $mconc*$A] -cMass;
# ABUTMENT ELEMENTS
element elasticBeamColumn 1000 1020 1010 [expr 1.0e+9] $Ec $Gc [expr 1.0e+9] [expr 1.0e+9] [expr 1.0e+9] 301;    # abut-1
element elasticBeamColumn 1001 1010 1030 [expr 1.0e+9] $Ec $Gc [expr 1.0e+9] [expr 1.0e+9] [expr 1.0e+9] 301;    # abut-1
element elasticBeamColumn 15000 15020 15010 [expr 1.0e+9] $Ec $Gc [expr 1.0e+9] [expr 1.0e+9] [expr 1.0e+9] 301;    # abut-15NE
element elasticBeamColumn 15001 15010 15030 [expr 1.0e+9] $Ec $Gc [expr 1.0e+9] [expr 1.0e+9] [expr 1.0e+9] 301;    # abut-15NE
element elasticBeamColumn 15002 15050 15040 [expr 1.0e+9] $Ec $Gc [expr 1.0e+9] [expr 1.0e+9] [expr 1.0e+9] 301;    # abut-15NR
element elasticBeamColumn 15003 15040 15060 [expr 1.0e+9] $Ec $Gc [expr 1.0e+9] [expr 1.0e+9] [expr 1.0e+9] 301;    # abut-15NR

