# written: fmk
# units: kip & in
wipe;

set pid [getPID]
set np [getNP]
set counter 0

set plotCounter 1

set pDampSecondary 0.00
set massTheta 1.0e-6;
set massY 1.0e-6;
set massSmallTheta 1.0e-6;
set massSmall 1.0e-6;


# set sinFactor .01
# set sinPeriod 3.797
set dtSin 0.01
set nPtSin 2000
set sinPeriod 3.0
set sinFactor 0.5
set periodStruct 3.797


# load procedures in other files
# source Steel2d.tcl
source Steel2dm.tcl
source ReadRecord.tcl;
source DisplayModel2D.tcl;
source DisplayPlane.tcl;
source PgPy_filtered_MRF_20Story_D_650_CanogaNS.tcl; #file contains gravity Pg/Py ratio for column hinge



# constants
set PI 3.14159
set in 1.0;

# select number of stories <---
set frame 20Story

set testType RelativeNormDispIncr
set Tol 1.0e-4;
set numStep 50;

set n 10.0
set pDamp 0.02; # damping ratio <---

# material properties
set Fyb 50.0;#36.0, 49.2
set Fyc 50.0;#50.0;
set E 29000.
set Ry 1.1;
#set b 0.01
set b 0.03
set g 386.4

set nFactorElem 10.

set matlab [open plotIt.m w]


if {$frame == "20Story" } {

    set numFrameResisting 2.0; # number of lateral load resisting frames
    set percentLoadFrame [expr 1.0/10.0];   # 5 bays in orthog dirn
    set percentLoadPDelta [expr 4.0/10.0];   # 5 bays in orthog dirn


        # set equiFloorDamper NO
    # roof: 90DL, 20LL, 25cladding; floor: 90DL, 50LL, 25cladding

    set firstfloorWeight [expr 0.05396940*2204.6226/$g*6/4]; # nodal mass from Lignos model, assign to each panel zone, not ea end of beam
    set roofWeight       [expr 0.05066770*2204.6226/$g*6/4];
    set floorWeight      [expr 0.05364600*2204.6226/$g*6/4];
#   puts "roof: $roofWeight floor: $floorWeight 1stfloor: $firstfloorWeight"

    set forceColfirstExt [expr 190.7000*0.2248];
    set forceColfirstInt [expr 127.1300*0.2248];
    set forceColtypExt   [expr 188.9800*0.2248];
    set forceColtypInt   [expr 125.9900*0.2248];
    set forceColroofExt  [expr 156.0000*0.2248];
    set forceColroofInt  [expr 103.9700*0.2248];

    set numMode 30

    # frame: dimensions and element properties
    set floorOffsets {180.    156.    156.    156.    156.    156.    156.    156.   156.    156.    156.    156.     156.    156.    156.   156.    156.    156.    156.    156.}
    set colOffsets   {240. 240. 240. }
    #                   1       2       3       4       5       6       7       8       9      10      11      12       13      14      15     16      17      18      19      20
    set colSizes     {W36X529 W36X529 W36X487 W36X487 W36X441 W36X441 W36X395 W36X395 W36X361 W36X361 W36X330 W36X330 W36X262 W36X262 W36X232 W36X232 W27X194 W27X194 W27X129 W27X129};
        # set colExtSizes  $colSizes
    set colExtSizes  {W14X500 W14X500 W14X455 W14X455 W14X455 W14X455 W14X370 W14X370 W14X370 W14X370 W14X311 W14X311 W14X283 W14X283 W14X233 W14X233 W14X193 W14X193 W14X132 W14X132}
    set beamSizes    {W33X130 W33X130 W33X141 W33X141 W33X141 W33X141 W33X152 W33X152 W33X152 W33X152 W33X152 W33X152 W33X152 W33X152 W33X130 W33X130 W33X130 W33X130 W24X68 W24X68};
    set beamExtSizes $beamSizes

    set BeamDepth   {33.1 33.1 33.3 33.3 33.3 33.3 33.5 33.5 33.5 33.5 33.5 33.5 33.5 33.5 33.1 33.1 33.1 33.1 23.7 23.7};
    set ColDepth    {39.8 39.8 39.3 39.3 38.9 38.9 38.4 38.4 38.  38.  37.7 37.7 36.9 36.9 37.1 37.1 28.1 28.1 27.6 27.6};
    set ExtColDepth {19.6 19.6 19.  19.  19.  19.  17.9 17.9 17.9 17.9 17.1 17.1 16.7 16.7 16.0 16.0 15.5 15.5 14.7 14.7};



    #                     1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
    set colSpliceCounter {0 0 1 0 0 0 1 0 0  0  1  0  1  0 1  0  1  0  1  0}
    # set colSpliceCounter {0 0 0 0 0 0 0 0 0  0  0  0  0  0 0  0  0  0  0  0}
    set spliceLength 54.;

    # set colSizesPD     {W14X550 W14X550 W14X550 W14X550 W14X455 W14X455 W14X455 W14X370 W14X370 W14X370 W14X311 W14X311 W14X311 W14X257 W14X257 W14X257 W14X176 W14X176 W14X176 W14X108 W14X108 W14X108};
    # set colExtSizesPD $colSizesPD
    # set beamSizesPD    {W14X22 W14X22  W14X22  W14X22 W14X22 W14X22 W14X22 W14X22 W14X22 W14X22 W14X22 W14X22  W14X22  W14X22  W14X22  W14X22 W14X22  W14X22  W14X22  W12X14};
    # set beamExtSizesPD $beamSizesPD;

    set pDelta YES; #FMK
    set areaPDelta 288000.0; #20ft x 20ft x 5
    set forcepDeltafirst [expr 2692.4940*0.2248];
    set forcepDeltatyp   [expr 2684.4110*0.2248];
    set forcepDeltaroof  [expr 2610.0000*0.2248];
    set factI 0.25
}



### YOU SHOULD NOT HAVE TO EDIT BELOW THIS LINE  #####


#calculated properties
set numFloor [expr [llength $floorOffsets]+1]


set numCline [expr [llength $colOffsets]+1]


for {set i 0; set width 0;} {$i < [expr $numCline-1]} {incr i 1} {set width [expr $width + [lindex $colOffsets $i]]}
# set massAtFloorNode [expr $floorWeight/($g*$numFrameResisting*$numCline*1.0)]
# set massAtRoofNode  [expr $roofWeight/($g*$numFrameResisting*$numCline*1.0)]
set uniformRoofLoad  [expr $roofWeight*$percentLoadFrame/$width]
set uniformFloorLoad [expr $floorWeight*$percentLoadFrame/$width]

# check of list dimensions for errors
#FMK if {[llength $colSizes] != [expr $numFloor-1]} {puts "ERROR: colSizes"; quit}
if {$numCline >= 10} {puts "ERROR: too many column lines, reprogram"; quit}

set tStart [clock clicks -milliseconds]

set g 386.4;
set counter 1

# ======= read external file for axial force in the external HSS columns due to gravity loads ======#
set inFilename ExtEleForcesFloorsRIgravityForceForce.out

if {$frame == "20Story"|$frame == "20StoryPre"} {
  if [catch {open $inFilename r} inFileID] {;                    # Open the input file and check for error
        puts stderr "Cannot open $inFilename for reading";        # output error statement
  } else {

    foreach line [split [read $inFileID] \n] {;             # Look at each line in the file

      if {[llength $line] == 0} {;                            # Blank line --> do nothing
        continue;
      } else {
        set timestep [expr [lindex $line 0]]

        if {$timestep == 1.0} {;                        # read the last step of gravity analysis result

          set counter 0
          for {set colLine 1} {$colLine <= $numCline} {incr colLine [expr $numCline-1]} {
            for {set floor1 1; set floor2 2} {$floor1 < $numFloor} {incr floor1 1; incr floor2 1} {

                set Ngrav [expr [lindex $line [expr 2+$counter*6]]]
                set counter [expr $counter+1]
                lappend HSSgrav $Ngrav      ;     # create list axial forces in external columns

            }
          }
        } else {
        }
      }
    }

    close $inFileID; ;        # Close the input file
  }
}

# ==== set equivalent floor damper properties
set CdampFactor 1.0;
if {$pDamp == 0.05} {
  set CdampValuesK {66.2115   49.1982   47.5585   46.8164   45.4449   46.1448   45.6503   45.0270   44.4242   43.5239   39.8839   38.0886   37.0317  33.0230   32.3421   31.2811   25.6903   23.1134   18.2987   13.2955};
} elseif {$pDamp == 0.02} {
  set CdampValuesK {26.4846   19.6793   19.0234   18.7266   18.1780   18.4579   18.2601   18.0108   17.7697   17.4096   15.9536   15.2354   14.8127  13.2092   12.9368   12.5124   10.2761    9.2454    7.3195    5.3182};
} elseif {$pDamp == 0.03} {
  set CdampValuesK {39.7609   36.9704   34.0338   33.3195   31.6563   31.7220   31.0890   30.5649   29.8126   28.7521   26.3484   25.2405   25.0444   23.0663   22.3096    20.8758   16.9277   14.2868    9.9453    6.2495};
} else {
  set CdampValuesK0 {39.7609   36.9704   34.0338   33.3195   31.6563   31.7220   31.0890   30.5649   29.8126   28.7521   26.3484   25.2405   25.0444   23.0663   22.3096    20.8758   16.9277   14.2868    9.9453    6.2495};
  set CdampValuesK $pDamp*$CdampValuesK0/0.03;
  puts "calc for the specific damping ratio"
}

# max from RHA M002: set FrValues {977.04        894.65        813.11        798.41        765.51        753.98        777.93        772.67        761.32        757.18        755.95        726.71        756.21        739            728.4        698.79        641.37        621.85        478.96        340.69}
# set FrValues {876.54        814.02        775.89        748.06        728.58        728.50        735.09        731.23        715.01        716.83        711.22        695.74        693.15        667.69        632.59        596.33        572.81        510.41        413.25        284.81};
# based on 2% modal
# set FrValues {1844.10        1605.11        1433.93        1326.85        1385.31        1509.90        1338.38        1321.06        1376.57        1392.37        1398.10        1403.03        1272.80        1203.95        1115.97        1184.66        1132.55        959.40        745.08        509.03};
# based on 1% equivFloorDamping under 3636PPhor
# set FrValues {1925.4        1492.0        1371.6        1279.3        1282.8        1298.3        1313.3        1331.9        1344.4        1351.1        1339.9        1271.2        1210.2        1126.5        1036.9        917.5        803.8        654.7        469.4        353.9};
# based on unscaled Simon GM
set FrValues {1388.8        1354.4        1332.3        1311.5        1282.0        1261.2        1245.2        1228.0        1206.0        1183.9        1154.4        1112.7        1060.0        996.2        927.5        852.6        763.1        639.2        477.2        271.1};



set zeta_target $pDamp

set Nfactor    {0.01 0.0788}; # 0.0788=386.4/980.6*0.2
set motionName {CanogaNS};


set dataDir TEST-LA20StoryPZCP-Fye50-$pDamp-LIN-test8
file mkdir $dataDir;
 #2 3 4 5 6 7 8 9 10 11
# equivFloorDampingCapped equivFloorDamping
foreach motionID {1} {

  set nfactor [lindex $Nfactor [expr $motionID -1]]; puts "nfactor = $nfactor"
  set motion [lindex $motionName [expr $motionID -1] ]; set LINfactor 1.0; puts $motion
  # set motion [lindex $motionNameLIN [expr $motionID -1] ]; set LINfactor 0.01
  set EIfactor 1.0; # set this factor only when studying the effect of rigidity of beam

  if {[expr $counter % $np] == $pid} {
    foreach colType {ComponentElement}  beamType {ComponentElement} {
     foreach damp {modal} {
       set totalMass 0


       wipe;
       model BasicBuilder -ndm 2 -ndf 3;  # Define the model builder, ndm = #dimension, ndf = #dofs

       set tStart1 [clock clicks -milliseconds]

       set nodeListMass {};
       set nodeListMassPZ {};
       set nodeListMassValues {};
       set nodeListReaction {};
#===# Build the Nodes===#
       for {set floor 1; set floorLoc 0} {$floor <= $numFloor} {incr floor 1} {
         if {$floor == $numFloor} {
             set massX $roofWeight
             # set massY [expr $massX*$percentLoadFrame]; # gravity cols take vertical
         } elseif {$floor == 2} {
             set massX $firstfloorWeight
             } else {
             set massX $floorWeight
             # set massY [expr $massX*$percentLoadFrame];
         }

         for {set colLine 1; set colLoc 0;} {$colLine <= $numCline} {incr colLine 1} {
           # node $colLine$floor $colLoc $floorLoc -mass $massX $massY $massTheta

           # add nodes for panel zone
           set p01 01; set p02 02; set p03 03; set p04 04; set p05 05; set p06 06; set p07 07; set p08 08; set p09 09; set p10 10;
           set p6 61; set p7 71; set p1 11; set p4 41;

           if {$floor == 1} {
             node $colLine$floor$p7 $colLoc $floorLoc -mass 0.0 0.0 0.0;
             # mass $colLine$floor$p05 0.0 0.0 0.0
             lappend nodeListReaction $colLine$floor$p7
             fix $colLine$floor$p7 1 1 1
           } else {
             # add panel zone nodes
             set aaa [lindex $BeamDepth [expr $floor -2]];
#                    puts "$floor $aaa"
             set pzvert [expr $aaa*0.5];
             set bbb [lindex $ColDepth [expr $floor -2]];
             set pzlat [expr $bbb/2.0];
             set ccc [lindex $ExtColDepth [expr $floor -2]];
             set pzextlat [expr 0.5*$ccc];

             if {$colLine == 1 || $colLine == $numCline} {
               node $colLine$floor$p01 [expr $colLoc - $pzextlat] [expr $floorLoc + $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p02 [expr $colLoc - $pzextlat] [expr $floorLoc + $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p03 [expr $colLoc + $pzextlat] [expr $floorLoc + $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p04 [expr $colLoc + $pzextlat] [expr $floorLoc + $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p05 [expr $colLoc + $pzextlat] [expr $floorLoc]  -mass $massX $massY $massTheta;
               node $colLine$floor$p06 [expr $colLoc + $pzextlat] [expr $floorLoc - $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p07 [expr $colLoc + $pzextlat] [expr $floorLoc - $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p08 [expr $colLoc - $pzextlat] [expr $floorLoc - $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p09 [expr $colLoc - $pzextlat] [expr $floorLoc - $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p10 [expr $colLoc - $pzextlat] [expr $floorLoc] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p6 [expr $colLoc ] [expr $floorLoc - $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p7 [expr $colLoc ] [expr $floorLoc + $pzvert] -mass $massSmall $massSmall $massSmallTheta

             } else {
               node $colLine$floor$p01 [expr $colLoc - $pzlat] [expr $floorLoc + $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p02 [expr $colLoc - $pzlat] [expr $floorLoc + $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p03 [expr $colLoc + $pzlat] [expr $floorLoc + $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p04 [expr $colLoc + $pzlat] [expr $floorLoc + $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p05 [expr $colLoc + $pzlat] [expr $floorLoc]  -mass $massX $massY $massTheta;
               node $colLine$floor$p06 [expr $colLoc + $pzlat] [expr $floorLoc - $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p07 [expr $colLoc + $pzlat] [expr $floorLoc - $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p08 [expr $colLoc - $pzlat] [expr $floorLoc - $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p09 [expr $colLoc - $pzlat] [expr $floorLoc - $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p10 [expr $colLoc - $pzlat] [expr $floorLoc] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p6 [expr $colLoc ] [expr $floorLoc - $pzvert] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p7 [expr $colLoc ] [expr $floorLoc + $pzvert] -mass $massSmall $massSmall $massSmallTheta
             }

               # add cover plate/RBS offset nodes
               # RBS parameter a=0.625bf, b=0.75db, c=0.25bf

             if {$colLine == 1} {
               set theSection [lindex $beamExtSizes [expr $floor -1]]

               set found 0
               foreach {section prop} [array get WSection $theSection] {
                 set propList [split $prop]
                 #AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
                 set db [expr [lindex $propList 2]*$in]
                 set bf [expr [lindex $propList 3]*$in]
                 set found 1
               }
               set phlat [expr 0.625*$bf+0.5*0.75*$db]
               node $colLine$floor$p1 [expr $colLoc + $pzextlat + $phlat] [expr $floorLoc];

             } elseif {$colLine == $numCline} {
               set theSection [lindex $beamExtSizes [expr $floor -1]]
               set found 0
               foreach {section prop} [array get WSection $theSection] {
                   set propList [split $prop]
                   #AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
                   set db [expr [lindex $propList 2]*$in]
                   set bf [expr [lindex $propList 3]*$in]
                   set found 1
               }
               set phlat [expr 0.625*$bf+0.5*0.75*$db]
               node $colLine$floor$p4 [expr $colLoc - $pzextlat - $phlat] [expr $floorLoc] -mass $massSmall $massSmall $massSmallTheta
             } else {
               set theSection [lindex $beamSizes [expr $floor -2]]
               set found 0
               foreach {section prop} [array get WSection $theSection] {
                   set propList [split $prop]
                   #AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
                   set db [expr [lindex $propList 2]*$in]
                   set bf [expr [lindex $propList 3]*$in]
                   set found 1
               }
               set phlat [expr 0.625*$bf+0.5*0.75*$db]
               node $colLine$floor$p1 [expr $colLoc + $pzlat + $phlat] [expr $floorLoc] -mass $massSmall $massSmall $massSmallTheta
               node $colLine$floor$p4 [expr $colLoc - $pzlat - $phlat] [expr $floorLoc] -mass $massSmall $massSmall $massSmallTheta
             }


             set totalMass [expr $totalMass+$massX]
             lappend nodeListMass $colLine$floor$p05
             lappend nodeListMassPZ $colLine$floor$p01 $colLine$floor$p02 $colLine$floor$p03 $colLine$floor$p04 $colLine$floor$p05 $colLine$floor$p06 $colLine$floor$p07 $colLine$floor$p08 $colLine$floor$p09 $colLine$floor$p10 $colLine$floor$p6 $colLine$floor$p7
             lappend nodeListMassValues $massX
             if {$colLine != 1} {
               equalDOF 1$floor$p05 $colLine$floor$p05 1;  # equalDOF applied to floor nodes <======
               # equalDOF 1$floor $colLine$floor 1 2
#                        puts "FMK - removed EQUAL DOF"
             }

           }

           if {$colLine < $numCline} {
               set colLoc [expr $colLoc + [lindex $colOffsets [expr $colLine-1]]]
           }
         }
         if {$floor < $numFloor} {
             set floorLoc [expr $floorLoc + [lindex $floorOffsets [expr $floor-1]]]
         }
       }

#        puts "---> nodal masses:"
#        puts "$nodeListMass";
#        puts "$nodeListMassValues";

       # define material
#        uniaxialMaterial Steel02 1 $Fy $E $b 20 0.925 0.15
       uniaxialMaterial Steel01 1 $Fyc $E $b;  # material for column ForceBeamCol
       uniaxialMaterial Steel01 2 $Fyb $E $b;  # material for beam ForceBeamCol

#===# build the columns
       geomTransf PDelta 1 ; #<====####
       set eleListHinge []
       set eleListComponent []
       set one 1
       set two 2

       set eleColListHinge []

       for {set colLine 1;set colLoc 0;} {$colLine <= $numCline} {incr colLine 1} {
         for {set floor1 1; set floor2 2; set floorLoc 0;} {$floor1 < $numFloor} {incr floor1 1; incr floor2 1} {

           set colSpliceOn [lindex $colSpliceCounter [expr $floor1 -1]]

#------- external columns
           if {$colLine == 1 || $colLine == $numCline} {
             set theSection    [lindex $colExtSizes [expr $floor1 -1]]
             set theSectionlow [lindex $colExtSizes [expr $floor1 -2]]
             if {$colType == "Displacement"} {
               # DispBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $theSection 1 1 -nip 5
             }
             if {$colType == "Force"} {
                 ForceBeamHSS2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection 1 1 -nip 5 ;#-elasticSection $E
             }
             if {$colType == "ForceCap"} {
                 ForceBeamHSS2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection 1 1 -nip 5 -capIt $pDamp $Fyc
             }
             if {$colType == "ForceWithHinges"} {
               # BeamWithHingesWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $theSection 1 1 -nip 5
             }
             if {$colType == "Elastic"} {
             # ElasticBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $theSection $E 1
               ElasticBeamHSSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection $E 1
             }
             if {$colType == "ComponentElement"} {
               # set axial [expr [lindex $HSSgrav [expr $floor1 -1]]]
               # puts $axial
               set location $floor1$floor2$colLine

               if {$colSpliceOn == 1} {
               # add node at column splice
                 node 909$colLine$floor1 [expr $colLoc ] [expr $floorLoc + $pzvert + $spliceLength]
                 # ComponentBeamWSection2d $colLine$floor1$colLine$floor2                  $colLine$floor1$p7 909$colLine$floor1 $theSectionlow $E $Fyc $Ry 1 -matType Bilin -Com_Type other-than-RBS
                 # ComponentBeamWSection2d 909$colLine$floor1$colLine$floor2 909$colLine$floor1 $colLine$floor2$p6 $theSection $E $Fyc $Ry 1 -matType Bilin -Com_Type other-than-RBS

                 set PgPyb [string cat $location b]
                 set PgPyt [string cat $location t]
                 ComponentColWSection2d $colLine$floor1$colLine$floor2                  $colLine$floor1$p7 909$colLine$floor1 $theSectionlow $E $Fyc $Ry [set $PgPyb] 1 -matType Bilin -nFactor $nFactorElem
                 ComponentColWSection2d 909$colLine$floor1$colLine$floor2 909$colLine$floor1 $colLine$floor2$p6 $theSection $E $Fyc $Ry [set $PgPyt] 1 -matType Bilin -nFactor $nFactorElem
                 lappend eleColListHinge $colLine$floor1$colLine$floor2
               } else {
                 set PgPy [string cat $location b]
  #              puts "$colLine$floor1$colLine$floor2"
                 ComponentColWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection $E $Fyc $Ry [set $PgPy] 1 -matType Bilin -nFactor $nFactorElem
                 # ComponentBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection $E $Fyc $Ry 1 -matType Bilin -Com_Type other-than-RBS
                 lappend eleColListHinge $colLine$floor1$colLine$floor2
               }
             }

   #------- internal columns
           } else {
             set theSection [lindex $colSizes [expr $floor1 -1]]
             set theSectionlow [lindex $colSizes [expr $floor1 -2]]
             if {$colType == "Displacement"} {
                 DispBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection 1 1 -nip 5
             }
             if {$colType == "Force"} {
                 ForceBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection 1 1 -nip 5
             }
             if {$colType == "ForceCap"} {
                 ForceBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection 1 1 -nip 5 -capIt $pDamp $Fyc
             }
             if {$colType == "ForceWithHinges"} {
                 BeamWithHingesWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection 1 1 -nip 5
             }
             if {$colType == "Elastic"} {
                 ElasticBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection $E 1
             }
             if {$colType == "Hinge"} {
                 BeamWithSteel01HingesWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection $E $Fyc $b 1 -nFactor $nFactorElem -doRayleigh $doRayleigh
             }

             if {$colType == "ComponentElement"} {
               set location $floor1$floor2$colLine

               if {$colSpliceOn == 1} {
                 # add node at column splice
                 node 909$colLine$floor1 [expr $colLoc ] [expr $floorLoc + $pzvert + $spliceLength]
                 # ComponentBeamWSection2d $colLine$floor1$colLine$floor2                  $colLine$floor1$p7 909$colLine$floor1 $theSectionlow $E $Fyc $Ry 1 -matType Bilin -metric -Com_Type other-than-RBS
                 # ComponentBeamWSection2d 909$colLine$floor1$colLine$floor2 909$colLine$floor1 $colLine$floor2$p6 $theSection $E $Fyc $Ry 1 -matType Bilin -metric -Com_Type other-than-RBS

                 # ComponentColWSection2d {eleTag iNode jNode sectType E Fy Ry PgPy transfTag args}
                 set PgPyb [string cat $location b]
                 set PgPyt [string cat $location t]
                 ComponentColWSection2d $colLine$floor1$colLine$floor2                  $colLine$floor1$p7 909$colLine$floor1 $theSectionlow $E $Fyc $Ry [set $PgPyb] 1 -matType Bilin -nFactor $nFactorElem
                 ComponentColWSection2d 909$colLine$floor1$colLine$floor2 909$colLine$floor1 $colLine$floor2$p6 $theSection $E $Fyc $Ry [set $PgPyt] 1 -matType Bilin -nFactor $nFactorElem
                 lappend eleColListHinge $colLine$floor1$colLine$floor2
#                puts "floor=$floor1 col=$colLine section=$theSectionlow $theSection"
               } else {
                 # ComponentBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection $E $Fyc $Ry 1 -matType Bilin -metric -Com_Type other-than-RBS
                 set PgPy [string cat $location b]
                 ComponentColWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection $E $Fyc $Ry [set $PgPy] 1 -matType Bilin -nFactor $nFactorElem
                 lappend eleColListHinge $colLine$floor1$colLine$floor2
               }
             }
           }

           if {$floor1 < $numFloor} {
             set floorLoc [expr $floorLoc + [lindex $floorOffsets [expr $floor1-1]]]
           }
         }
         if {$colLine < $numCline} {
           set colLoc [expr $colLoc + [lindex $colOffsets [expr $colLine-1]]]
         }
       }

       set nodeListInternal []
#      puts "colHinge $eleColListHinge"

#====# build the beams
       geomTransf Linear 2

       for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
         for {set colLine1  1; set colLine2 2} {$colLine1 < $numCline} {incr colLine1 1; incr colLine2 1} {
           if {$colLine1 == 1 || $colLine2 == $numCline} {
             set theSection [lindex $beamExtSizes [expr $floor -2]]
           } else {
             set theSection [lindex $beamSizes [expr $floor -2]]
           }
 # # ------ # ----- other cases
           if {$beamType == "ComponentElement"} {
             ComponentBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection $E $Fyb $Ry 2 -matType Bilin -Com_Type RBS -nFactor $nFactorElem
             # ComponentBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection $E $Fyb 2 -matType MultiLinear -metric
             lappend eleListHinge $colLine1$floor$colLine2$floor
           }
           if {$beamType == "Displacement"} {
             DispBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection 2 2
             lappend eleListHinge $colLine1$floor$colLine2$floor
           }
           if {$beamType == "Force"} {
             ForceBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection 2 2
             lappend eleListHinge $colLine1$floor$colLine2$floor
           }
           if {$beamType == "ForceCap"} {
             ForceBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection 2 2
             lappend eleListHinge $colLine1$floor$colLine2$floor -capIt $pDamp $Fyb
           }
           if {$beamType == "ForceWithHinge"} {
             BeamWithHingesWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection 2 2
           }
           if {$beamType == "Elastic"} {
             ElasticBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection $E 2
           }
           if {$beamType == "HingeMultiLinear1"  || $beamType == "HingeMultiLinear10"  ||
               $beamType == "HingeMultiLinear50" || $beamType == "HingeMultiLinear100" ||
               $beamType == "HingeML1000" || $beamType == "HingeML1000NoR"} {

             set doRayleigh 1

             if {$beamType == "HingeML1000NoR"} {
               set nFactorElem 1000
               set doRayleigh 0
             }

             BeamWithConcentratedHingesWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection [expr $EIfactor*$E] [expr $EIfactor*$Fyb] 2 -nFactor $nFactorElem -doRayleigh $doRayleigh -matType MultiLinear -metric
             lappend eleListHinge $colLine1$floor$colLine2$floor$one
             lappend eleListHinge $colLine1$floor$colLine2$floor$two
           }

           if {$beamType == "HingeBilin"} {
               BeamWithConcentratedHingesWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection $E $Fyb 2 -nFactor $nFactorElem -doRayleigh $doRayleigh -matType Bilin02 -metric
               lappend eleListHinge $colLine1$floor$colLine2$floor$one
               lappend eleListHinge $colLine1$floor$colLine2$floor$two
           }
           if {$beamType == "HingeMultiLinear"} {
               BeamWithConcentratedHingesWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection $E $Fyb 2 -nFactor $nFactorElem -doRayleigh $doRayleigh -matType MultiLinear -metric
               lappend eleListHinge $colLine1$floor$colLine2$floor$one
               lappend eleListHinge $colLine1$floor$colLine2$floor$two
           }
                 # }
         }
       }


#    puts "build PZ..."
     #====#====# add panel zone elements
     set Apz 1000.0;
     set Ipz 1.0e5;
     #                2   3   4   5    6  7   8   9   10    11     12   13   14   15    16      17    18   19  20     21
     set tdpExtList {0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.125 0.125 0.25  0.25 0.375 0.375 0.3125 0.8125 1.0 1.0  0.0    0.0};
     set tdpIntList {0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0   0.0   0.0   0.0   0.0  0.0   0.125 0.625   1.0 1.0  0.0625 0.0625};
     set aa 1;
     set PZspringList [];
     set PZrigidList [];
     for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
       for {set colLine 1} {$colLine <= $numCline} {incr colLine 1} {
         if {$colLine == 1 || $colLine == $numCline} {
           set theSection [lindex $colExtSizes [expr $floor -1]]
           set tdpExt [lindex $tdpExtList [expr $floor -1]];
           # puts "$theSection"
           set found 0

           foreach {section prop} [array get WSection $theSection] {
             set propList [split $prop]

             #AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
             set A [expr [lindex $propList 1]*$in*$in]
             set dc [expr [lindex $propList 2]*$in]
             set bf_c [expr [lindex $propList 3]*$in]
             set tw [expr [lindex $propList 4]*$in]
             set tf_c [expr [lindex $propList 4]*$in]
             set Ixx [expr [lindex $propList 5]*$in*$in*$in*$in]
             set found 1
             set tp [expr $tw+$tdpExt];
             # puts "$dc $bf_c $tf_c $tp"
           }
         } else {
           set tdpInt [lindex $tdpIntList [expr $floor -1]];
           set theSection [lindex $colSizes [expr $floor -1]];
           # puts "$theSection"
           set found 0

           foreach {section prop} [array get WSection $theSection] {
             set propList [split $prop]

             #AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
             set A [expr [lindex $propList 0]*$in*$in]
             set dc [expr [lindex $propList 1]*$in]
             set bf_c [expr [lindex $propList 2]*$in]
             set tw [expr [lindex $propList 3]*$in]
             set tf_c [expr [lindex $propList 4]*$in]
             set Ixx [expr [lindex $propList 5]*$in*$in*$in*$in]

             set found 1
             set tp [expr $tw + $tdpInt];
             # puts "$dc $bf_c $tf_c $tp"
           }
         }


         #--- rigid element
         # elemPanelZone2D {eleID nodeR E A_PZ I_PZ transfTag}
         elemPanelZone2D 500$colLine$floor$aa $colLine$floor$p01 $E $Apz $Ipz 1;
         lappend PZrigidList [expr 500$colLine$floor$aa+2] [expr 500$colLine$floor$aa+7] [expr 500$colLine$floor$aa+3] [expr 500$colLine$floor$aa+6]

         #--- panel zone spring
         set Ry 1.2;
         set as_PZ 0.1; # J.Hall used 10% strain hardening for panel zone
         set pzspring 00;

         set thebeamSection [lindex $beamSizes [expr $floor -1]];
         # puts "$thebeamSection"
         set found 0

         foreach {section prop} [array get WSection $thebeamSection] {
           set propList [split $prop]

           #AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
           set db [expr [lindex $propList 1]*$in]
           set found 1
         }

         # puts "Joint$colLine$floor: $dc $bf_c $tf_c $tp $db"

         # rotPanelZone2D {eleID nodeR nodeC E Fy dc bf_c tf_c tp db Ry as}
         rotPanelZone2D 5$colLine$floor$pzspring $colLine$floor$p03 $colLine$floor$p04 $E $Fyc $dc $bf_c $tf_c $tp $db $Ry $as_PZ;

         lappend PZspringList 5$colLine$floor$pzspring

       }
     }
#puts "$PZspringList"

# puts "build elastic section btn col and RBS..."
#====#====# add cover plate/RBS offset element
set K44_two 2.0625; #K44 = 6*(1+n)/(2+3*n)
set K11_two 3.9375; #K11 = K33 = (1+2*n)*K44/(1+n)
set K33_two 3.9375;
set K44_one 1.9355; #K44 = 6*n/(1+3*n )
set K11_one 3.9375; #K11 = (1+2*n)*K44/(1+n)
set K33_one 3.8710; #K33 = 2*K44 for one end spring

for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
 for {set colLine1  1; set colLine2 2} {$colLine1 < $numCline} {incr colLine1 1; incr colLine2 1} {

     if {$colLine1 == 1 || $colLine2 == $numCline} {
         set theSection [lindex $beamExtSizes [expr $floor -2]]
     } else {
         set theSection [lindex $beamSizes [expr $floor -2]]
     }

     set found 0

     foreach {section prop} [array get WSection $theSection] {
       set propList [split $prop]

       #AISC_Manual_Label A d bf tw tf Ix Iy Zx Sx rx Zy Sy ry J
       set A [expr [lindex $propList 0]*$in*$in]
       set d [expr [lindex $propList 1]*$in]
       set bf [expr [lindex $propList 2]*$in]
       set Ixx [expr [lindex $propList 5]*$in*$in*$in*$in]

       set found 1
     }
     # build elastic element
     # set coverplatefactor 1.0;
     element elasticBeamColumn 201$colLine1$floor $colLine1$floor$p05 $colLine1$floor$p1 [expr 1.0*$A] $E [expr 1.0*$Ixx] 2;
     element elasticBeamColumn 202$colLine1$floor $colLine2$floor$p4 $colLine2$floor$p10 [expr 1.0*$A] $E [expr 1.0*$Ixx] 2;

     # element ModElasticBeam2d $eleTag $iNode $jNode $A $E $Iz $K11 $K33 $K44 $transfTag <-mass $massDens> <-cMass>
     # element ModElasticBeam2d 201$colLine1$floor $colLine1$floor$p05 $colLine1$floor$p1 $A $E [expr 1.0*$Ixx] $K11_one $K33_one $K44_one 2;
     # element ModElasticBeam2d 202$colLine1$floor $colLine2$floor$p4 $colLine2$floor$p10 $A $E [expr 1.0*$Ixx] $K11_one $K33_one $K44_one 2;

     # set K44_two 2.0625;
     # set K11_two 3.9375;
     # set K33_two 3.9375;
     # set K44_one 1.9355;
     # set K11_one 3.9375;
     # set K33_one 3.8710;

     ## element ModElasticBeam2d $eleTag $iNode $jNode $A $E $Iz $K11 $K33 $K44 $transfTag <-mass $massDens> <-cMass>
     # element ModElasticBeam2d  5002104 21040 2104 $A_Beam12_2 $E [expr (1+(((2672.8250-497.8250)/2672.8250)**2)/10.00)*0.8100*$Ix_Beam12_2] $K11_one $K33_one $K44_one 1;

 }
}

puts "BEFORE"
   constraints Plain
   numberer RCM
   system ProfileSPD; #UmfPack;  # <----
          system SuperLU
   test $testType $Tol $numStep
   algorithm Newton
   integrator Newmark 0.5 0.25
   analysis Transient        ; #Transient

set lambda [eigen $numMode]
puts " Period T1: [expr 2*$PI/sqrt([lindex $lambda 0])] sec [lindex $lambda 0]"
puts " Period T2: [expr 2*$PI/sqrt([lindex $lambda 1])] sec [lindex $lambda 1]"
puts " Period T3: [expr 2*$PI/sqrt([lindex $lambda 2])] sec [lindex $lambda 2]"


# --- # --- Display model for checking
       # DisplayModel2D { {ShapeType nill} {dAmp 5}  {xLoc 10} {yLoc 10} {xPixels 512} {yPixels 384} {nEigen 1} }
#        DisplayModel2D nill 50 10 10 1024 1024 1;
#        after 3000;


#====# === add pDelta columns
       if {$pDelta == "YES"} {

           if {$areaPDelta != 0.} {
               set pFrame 12;
               set ptRoofLoad  $forcepDeltaroof;
               set ptFloorLoad $forcepDeltatyp;
               set ptfirstFloorLoad $forcepDeltafirst;

               # add nodes for PDelta
               set colLine [expr $numCline+1]
               set colLoc [expr $colLoc+50.]
               for {set floor 1; set floorLoc 0} {$floor <= $numFloor} {incr floor 1} {

                   node $colLine$floor $colLoc $floorLoc
                   if {$floor == 1} {
                       fix $colLine$floor 1 1 0
                       lappend nodeListReaction $colLine$floor
                   }  else {
                       equalDOF 1$floor$p05 $colLine$floor 1;# <=====
                       # equalDOF 1$floor $colLine$floor 1 2
                   }
                   if {$floor < $numFloor} {
                       set floorLoc [expr $floorLoc + [lindex $floorOffsets [expr $floor-1]]]
                   }
               }

               # add Pdelta elements
               geomTransf PDelta 3; #<===== ###


               for {set floor1 1; set floor2 2} {$floor1 < $numFloor} {incr floor1 1; incr floor2 1} {
                   element elasticBeamColumn $colLine$floor1$colLine$floor2 $colLine$floor1 $colLine$floor2 $areaPDelta $E 1.0e-1 3
               }
               # add the gravity loads
               pattern Plain 110 Linear {
                   for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
                               if {$floor == $numFloor} {
                                       load $colLine$floor 0. -$ptRoofLoad 0.
                               } elseif {$floor == 2} {
                                       load $colLine$floor 0. -$ptfirstFloorLoad 0.
                               } else {
                                       load $colLine$floor 0. -$ptFloorLoad 0.
                               }
                   }
               }
           } else { ; # NEW PDETA

               set pFrame 12;
               set ptRoofLoad  $forcepDeltaroof;
               set ptFloorLoad $forcepDeltatyp;
               set ptfirstFloorLoad $forcepDeltafirst;

               pattern Plain 110 Linear ; #{
                   # loads added after this will be added to this load pattern
               #}
               set nodeListCLine1Grav []
               for {set floor 1; set floorLoc 0} {$floor <= $numFloor} {incr floor 1} {
                   for {set colLine 1; set colLoc 0;} {$colLine <= $numCline} {incr colLine 1} {
                       set p05 05;

                       if {$floor == 1} {
                           node $colLine$floor$pFrame $colLoc $floorLoc
                           lappend nodeListReaction $colLine$floor$p7
                           fix $colLine$floor$pFrame 1 1 1
                       }  else {
                           node $colLine$floor$pFrame $colLoc $floorLoc
                           equalDOF 1$floor$p05 $colLine$floor$pFrame 1;

                           # add pDelta load
                           if {$floor == $numFloor} {
                               load $colLine$floor$pFrame 0. -$ptRoofLoad 0.
                           } elseif {$floor == 2} {
                               load $colLine$floor$pFrame 0. -$ptfirstFloorLoad 0.
                           } else {
                               load $colLine$floor$pFrame 0. -$ptFloorLoad 0.
                           }
                       }

                       lappend nodeListCLine1Grav 1$floor$pFrame

                       if {$colLine < $numCline} {
                           set colLoc [expr $colLoc + [lindex $colOffsets [expr $colLine-1]]]
                       }
                   }
                   if {$floor < $numFloor} {
                       set floorLoc [expr $floorLoc + [lindex $floorOffsets [expr $floor-1]]]
                   }
               }
               # add gravity frame columns
               for {set colLine 1} {$colLine <= $numCline} {incr colLine 1} {
                   for {set floor1 1; set floor2 2} {$floor1 < $numFloor} {incr floor1 1; incr floor2 1} {

                       #------- external columns
                       if {$colLine == 1 || $colLine == $numCline} {
                           set theSection [lindex $colExtSizesPD [expr $floor1 -1]]
                       } else {
                           set theSection [lindex $colSizesPD [expr $floor1 -1]]
                       }

                       ElasticBeamWSection2d $colLine$floor1$colLine$floor2$pFrame $colLine$floor1$pFrame $colLine$floor2$pFrame $theSection $E 1 -factI $factI
                   }
               }
               # add gravity frame beams
               for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
                   for {set colLine1  1; set colLine2 2} {$colLine1 < $numCline} {incr colLine1 1; incr colLine2 1} {

                       if {$colLine1 == 1 || $colLine2 == $numCline} {
                           set theSection [lindex $beamExtSizesPD [expr $floor -2]]
                       } else {
                           set theSection [lindex $beamSizesPD [expr $floor -2]]
                       }

                       ElasticBeamWSection2d $colLine1$floor$colLine2$floor$pFrame $colLine1$floor$pFrame $colLine2$floor$pFrame $theSection $E 2 -factI $factI
                   }
               }
           }
       }

#===# add beam uniform loads
       # pattern Plain 101 Linear {
           # for {set colLine1  1} {$colLine1 < $numCline} {incr colLine1 1} {
               # set colLine2 [expr $colLine1 + 1]
               # for {set floor 2} {$floor <= $numFloor} {incr floor 1} {
                   # if {$floor == $numFloor} {
                       # eleLoad -ele $colLine1$floor$colLine2$floor -type beamUniform -$uniformRoofLoad
                       # eleLoad -ele 201$colLine1$floor -type beamUniform -$uniformRoofLoad
                       # eleLoad -ele 202$colLine1$floor -type beamUniform -$uniformRoofLoad
                   # } else {
                       # eleLoad -ele $colLine1$floor$colLine2$floor -type beamUniform -$uniformFloorLoad
                       # eleLoad -ele 201$colLine1$floor -type beamUniform -$uniformFloorLoad
                       # eleLoad -ele 202$colLine1$floor -type beamUniform -$uniformFloorLoad
                   # }
               # }
           # }
       # }


#===# add column axial loads
puts "---> gravity loads"

pattern Plain 102 Linear {
       for {set floor2 2} {$floor2 <= $numFloor} {incr floor2 1} {
             for {set colLine  1} {$colLine1 < $numCline} {incr colLine 1} {
                 if {$colLine == 1 || $colLine == $numCline} {
                     if {$floor == $numFloor} {
                             load $colLine$floor2$p6 0. -$forceColroofExt 0.
                     } elseif {$floor == 2} {
                             load $colLine$floor2$p6 0. -$forceColfirstExt 0.
                     } else {
                             load $colLine$floor2$p6 0. -$forceColtypExt 0.
                     }
                 } else {
                     if {$floor == $numFloor} {
                             load $colLine$floor2$p6 0. -$forceColroofInt 0.
                     } elseif {$floor == 2} {
                             load $colLine$floor2$p6 0. -$forceColfirstInt 0.
                     } else {
                             load $colLine$floor2$p6 0. -$forceColtypInt 0.
                     }
                 }
             }
       }

       }
# puts "---> gravity loads"
# puts "roof:$uniformRoofLoad...floor:$uniformFloorLoad...ptRoofLoad:$ptRoofLoad ...ptFloorLoad:$ptFloorLoad"


           # constraints Plain
           # numberer RCM
           # system ProfileSPD; #UmfPack ProfileSPD;  # <----
           # system SuperLU
           # test $testType $Tol $numStep
           # algorithm Newton
           # integrator Newmark 0.5 0.25
           # analysis Transient        ; #Transient
set lambda [eigen $numMode]
puts " Period T1: [expr 2*$PI/sqrt([lindex $lambda 0])] sec [lindex $lambda 0]"
puts " Period T2: [expr 2*$PI/sqrt([lindex $lambda 1])] sec [lindex $lambda 1]"
puts " Period T3: [expr 2*$PI/sqrt([lindex $lambda 2])] sec [lindex $lambda 2]"


# =========== analysis settings ============ #
           # Gravity-analysis: load-controlled static analysis
           constraints Plain;
           numberer RCM;
           system BandGeneral;
           test NormUnbalance 1.0e-6 10;
           algorithm Newton;
           integrator LoadControl 0.1;
           analysis Static;
           analyze 10
               # printA -file stiffnessmatrix.dat; puts "print stiffness matrix to file..."
           # printA;
           # maintain constant gravity loads and reset time to zero
           loadConst -time 0.0

puts "AFTER GRAVITY"
wipeAnalysis
           constraints Plain
           numberer RCM
           system UmfPack; #UmfPack ProfileSPD;  # <----
#            system SuperLU
           test $testType $Tol $numStep
           algorithm Newton
           integrator Newmark 0.5 0.25
           analysis VariableTransient;
set lambda [eigen $numMode]
puts " Period T1: [expr 2*$PI/sqrt([lindex $lambda 0])] sec [lindex $lambda 0]"
puts " Period T2: [expr 2*$PI/sqrt([lindex $lambda 1])] sec [lindex $lambda 1]"
puts " Period T3: [expr 2*$PI/sqrt([lindex $lambda 2])] sec [lindex $lambda 2]"


               # for 10% initial stiffness and a possibility of dt = dt/1000 for some iterations:

               # algorithm NewtonHallM 0.1
               # analysis Transient -numSubLevels 3 -numSubSteps 10




# ----- # ---- add some damping ---- #

# ----- # ---- rayleigh damping using initial stiffness

           if {$damp == "RI"} {

               set numMode 30; #<-----
               set lambda [eigen $numMode]
               # puts "eigenvalues: $lambda"
               set omegaI [expr pow([lindex $lambda 0],0.5)];
               set omegaJ [expr pow([lindex $lambda 2],0.5)];
               set alphaM [expr $pDamp*(2*$omegaI*$omegaJ)/($omegaI+$omegaJ)];
#                set betaK [expr 2.0*$pDamp/$omegaI]
#                puts $betaK

               set betaK [expr 2.*$pDamp/($omegaI+$omegaJ)];

               # puts "a0: $alphaM a1: $betaK"

       # write frequency data
       set fmkOut [open freq.data w]
               for {set fmk 0} {$fmk < $numMode} {incr fmk 1} {
                   puts "T[expr $fmk+1]: [expr 2*$PI/sqrt([lindex $lambda $fmk])]sec w[expr $fmk+1]: [expr sqrt([lindex $lambda $fmk])]rad/sec"
               puts $fmkOut "[expr 2*$PI/sqrt([lindex $lambda $fmk])] sec [expr sqrt([lindex $lambda $fmk])]"
               }
   close $fmkOut

               puts "ALPHA M: $alphaM BETAK $betaK"
               puts " Period T1: [expr 2*$PI/sqrt([lindex $lambda 0])] sec [lindex $lambda 0]"
               puts " Period T2: [expr 2*$PI/sqrt([lindex $lambda 1])] sec [lindex $lambda 1]"
               puts " Period T3: [expr 2*$PI/sqrt([lindex $lambda 2])] sec [lindex $lambda 2]"


       # write first mode shape data ================================================

#                set alphaM 0.
#                set betaK [expr 2.0*$pDamp/$omegaI]
#                puts "ALPHA M: $alphaM BETAK $betaK omegaI: $omegaI omegaJ: $omegaJ"

#                set betaK 0
#                set alphaM [expr 2.0*$pDamp*$omegaI]

          #rayleigh $alphaM $betaK $betaKinit $betaKcomm
               rayleigh $alphaM  0.      $betaK       0.
           }


# ----- # ---- add some damping, assume alphaM=0
           if {$damp == "TI"} {
               set lambda [eigen 3]
               puts "eigenvalues: $lambda"
               set omegaI [expr pow([lindex $lambda 0],0.5)];
               set omegaJ [expr pow([lindex $lambda 2],0.5)];
               set alphaM 0.
               set betaK [expr 2.0*$pDamp/$omegaI]
               puts "ALPHA M: $alphaM BETAK $betaK omegaI: $omegaI omegaJ: $omegaJ"
      #rayleigh $alphaM $betaK $betaKinit $betaKcomm
               rayleigh $alphaM  0.      $betaK       0.
           }



# ----- # ---- rayleigh damping using the current stiffness
           if {$damp == "RC"} {
               set lambda [eigen 3]
               puts "eigenvalues: $lambda"
               set omegaI [expr pow([lindex $lambda 0],0.5)];
               set omegaJ [expr pow([lindex $lambda 2],0.5)];
               set alphaM [expr $pDamp*(2*$omegaI*$omegaJ)/($omegaI+$omegaJ)];
               set betaK [expr 2.*$pDamp/($omegaI+$omegaJ)];

#                set alphaM 0.
#                set betaK [expr 2.0*$pDamp/$omegaI]
#                puts "ALPHA M: $alphaM BETAK $betaK"

#                set betaK 0
#                set alphaM [expr 2.0*$pDamp*$omegaI]

          #rayleigh $alphaM $betaK $betaKinit $betaKcomm
               rayleigh $alphaM $betaK 0. 0.
           }

# ----- # ---- modal damping
           if {$damp == "M" } {
               set numMode 20; #<-----
               set lambda [eigen $numMode]
               puts "LAMBDA $lambda"
               modalDampingQ $pDamp

               puts " Period T1: [expr 2*$PI/sqrt([lindex $lambda 0])] sec [lindex $lambda 0]"
               puts " Period T2: [expr 2*$PI/sqrt([lindex $lambda 1])] sec [lindex $lambda 1]"
               puts " Period T3: [expr 2*$PI/sqrt([lindex $lambda 2])] sec [lindex $lambda 2]"

           }


               # ----- # ---- modal damping
           if {$damp == "MRI" } {

               # set omegaI [expr pow([lindex $lambda 0],0.5)];
               # set omegaJ [expr pow([lindex $lambda 2],0.5)];
               # set alphaM [expr $pDamp*(2*$omegaI*$omegaJ)/($omegaI+$omegaJ)];
               # set betaK [expr 2.*$pDamp/($omegaI+$omegaJ)];
               # puts "a0: $alphaM a1: $betaK"

               set numMode 30; #<-----
               set lambda [eigen $numMode]
               puts "LAMBDA $lambda"
               modalDamping 0.0200 0.0155 0.0200 0.0261 0.0333 0.0415 0.0506 0.0611 0.0721 0.0838 0.0963 0.1100 0.1164 0.1236 0.1243 0.1384 0.1533 0.1555 0.1710 0.1902 0.1913 0.2109 0.2348 0.2352 0.2568 0.2740 0.3062 0.3090 0.3254 0.3484;
           }

               if {$damp == "NoDamping"} {
                       set lambda [eigen 3]
               }

# ----- # ---- inter-story damper
       if {$damp == "equivFloorDamping"} {

           set lambda [eigen 30]

           set colL 3;
           set colR 4;


           set FloorDamperList []

           for {set floor1 1; set floor2 2} {$floor1 < $numFloor} {incr floor1 1; incr floor2 1} {

               set Cdamp [expr [lindex $CdampValuesK [expr $floor1 -1]]];#$CdampValuesK


               set Fr [expr 2*$pDamp*[lindex $FrValues [expr $floor1 -1]]]

               # uniaxialMaterial BilinearOilDamper $matTag $K $Cd <$Fr $p> <$LGap> < $NM $RelTol $AbsTol $MaxHalf>
               # uniaxialMaterial BilinearOilDamper $colL$floor1$colR$floor2 $Kdamp $Cdamp $Fr 0.0;

               #c                                 uniaxialMaterial Elastic $colL$floor1$colR$floor2 $Kdamp

               # element twoNodeLink $eleTag $iNode $jNode -mat $matTags -dir $dirs <-orient <$x1 $x2 $x3> $y1 $y2 $y3>
               # element  twoNodeLink   $colL$floor1$colR$floor2 $colL$floor1 $colR$floor2  -mat  $colL$floor1$colR$floor2 -dir 1


               #======= new dampermaterial ========#
               # Linear F-V relation
               uniaxialMaterial Elastic $colL$colL$floor1$colR$floor2 $Cdamp
               uniaxialMaterial DamperMaterial $colL$floor1$colR$floor2 $colL$colL$floor1$colR$floor2

               if {$floor1 == 1} {
               element zeroLength $colL$floor1$colR$floor2 $colL$floor1$p7 $colR$floor2$p05  -mat  $colL$floor1$colR$floor2 -dir 1
               } else {
               element zeroLength $colL$floor1$colR$floor2 $colL$floor1$p05 $colR$floor2$p05  -mat  $colL$floor1$colR$floor2 -dir 1
               }
               #======= new dampermaterial ========#

               lappend FloorDamperList $colL$floor1$colR$floor2;
           }

           puts "FloorDamperList=$FloorDamperList"

       }

       if {$damp == "equivFloorDampingCapped"} {

           set lambda [eigen 30]

           set colL 3
           set colR 4;

           set FloorDamperList []

           for {set floor1 1; set floor2 2} {$floor1 < $numFloor} {incr floor1 1; incr floor2 1} {

               set Cdamp [expr [lindex $CdampValuesK [expr $floor1 -1]]];#$CdampValuesK

               set Fr [expr 2*$pDamp*[lindex $FrValues [expr $floor1 -1]]]

               #======= new dampermaterial ========#
               # Bilinear F-V relation
               #uniaxialMaterial Steel01 $matTag $Fy $E0 $b <$a1 $a2 $a3 $a4>
               uniaxialMaterial Steel01 $colL$colL$floor1$colR$floor2 $Fr $Cdamp 0.0;

               #uniaxialMaterial ElasticBilin $matTag $EP1 $EP2 $epsP2 <$EN1 $EN2 $epsN2>
#                uniaxialMaterial ElasticBilin $colL$colL$floor1$colR$floor2 $Cdamp 0.0 [expr $Fr/$Cdamp]

               uniaxialMaterial DamperMaterial $colL$floor1$colR$floor2 $colL$colL$floor1$colR$floor2

               if {$floor1 == 1} {
               element zeroLength $colL$floor1$colR$floor2 $colL$floor1$p7 $colR$floor2$p05  -mat  $colL$floor1$colR$floor2 -dir 1
               } else {
               element zeroLength $colL$floor1$colR$floor2 $colL$floor1$p05 $colR$floor2$p05  -mat  $colL$floor1$colR$floor2 -dir 1
               }

               #======= new dampermaterial ========#

               lappend FloorDamperList $colL$floor1$colR$floor2;
           }

#            puts "FloorDamperList=$FloorDamperList"

       }



       # DisplayModel2D 1 50 10 10 1024 1024 2;
       # after 3000;

       # # write first mode shape data ================================================
       # puts "write mode shape"
       # set phiOut [open $dataDir/modeshape.data w]
       # nodeEigenvector $nodeTag $eigenvector <$dof>
                   # set C1 1;
                               # for {set floor 1} {$floor <= $numFloor} {incr floor 1} {

                               # set phihC1($floor) [nodeEigenvector $C1$floor$p05 1 1];
                               # set phivC1($floor) [nodeEigenvector $C1$floor$p05 1 2]
                               # puts "floor$floor $phihC1($floor)"
                               # }
                               # foreach {key value} [array get phihC1] {
                               # puts $phiOut "horizontal $key => $value"
                               # }
                               # foreach {key value} [array get phivC1] {
                               # puts $phiOut "vertical $key => $value"
                               # }
       # close $phiOut





               # puts " Period T1: [expr 2*$PI/sqrt([lindex $lambda 0])] sec [lindex $lambda 0]"
               # puts " Period T2: [expr 2*$PI/sqrt([lindex $lambda 1])] sec [lindex $lambda 1]"
               # puts " Period T3: [expr 2*$PI/sqrt([lindex $lambda 2])] sec [lindex $lambda 2]"
puts "add load"
# ==== create a load pattern for uniform excitation ==== #
# ---- read record
       if {$motion == "elCentro2" || $motion=="elCentro3" || $motion == "elCentro4"  || $motion=="elCentro.25" || $motion == "elCentro.5"} {
           ReadRecord elCentro $motion.g3 dt nPt;
       } elseif { $motion == "se3001" || $motion == "se305"} {
           ReadRecord se30 $motion.g3 dt nPt;
       } elseif {$motion == "se30"} {
           ReadRecord $motion $motion.g3 dt nPt;
               } elseif {$motion == "RSN18501"} {
           ReadRecord RSN185.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN73601"} {
           ReadRecord RSN736.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN83801"} {
           ReadRecord RSN838.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN119301"} {
           ReadRecord RSN1193.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN122301"} {
           ReadRecord RSN1223.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN127701"} {
           ReadRecord RSN1277.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN143001"} {
           ReadRecord RSN1430.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN155101"} {
           ReadRecord RSN1551.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN695301"} {
           ReadRecord RSN6953.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN118301"} {
           ReadRecord RSN1183.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN154101"} {
           ReadRecord RSN1541.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN121301"} {
           ReadRecord RSN1213.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN123701"} {
           ReadRecord RSN1237.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN126101"} {
           ReadRecord RSN1261.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN143601"} {
           ReadRecord RSN1436.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN147101"} {
           ReadRecord RSN1471.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN149001"} {
           ReadRecord RSN1490.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN154501"} {
           ReadRecord RSN1545.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN265501"} {
           ReadRecord RSN2655.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN367001"} {
           ReadRecord RSN3670.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN374701"} {
           ReadRecord RSN3747.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "RSN582501"} {
           ReadRecord RSN5825.AT2 $motion.g3 dt nPt;
               } elseif {$motion == "LA3536PPhor01"} {
           ReadRecord LA3536PPhor.AT2 $motion.g3 dt nPt; puts "dt=$dt nPt=$nPt"
               } elseif {$motion == "CanogaNS01" || $motion == "CanogaNS"} {
                   set dt 0.0100;
                   set nPt 2502;
               } else {
                   # ReadRecord $motion $motion.g3 dt nPt;
                   ReadRecord  $motion.AT2 $motion.g3 dt nPt;
       }
# ---- define time series
       if {$motion == "elCentro"} {
           timeSeries Path 10 -filePath $motion.g3 -dt $dt -factor $g
       } elseif {$motion == "elCentro2"} {
           timeSeries Path 10 -filePath $motion.g3 -dt $dt -factor [expr $g*2.0]
       } elseif {$motion == "elCentro3"} {
           timeSeries Path 10 -filePath $motion.g3 -dt $dt -factor [expr $g*3.0]
       } elseif {$motion == "elCentro4"} {
           timeSeries Path 10 -filePath $motion.g3 -dt $dt -factor [expr $g*4.0]
       } elseif {$motion == "elCentro.25"} {
           timeSeries Path 10 -filePath $motion.g3 -dt $dt -factor [expr $g*0.25]
       } elseif {$motion == "elCentro.5"} {
           timeSeries Path 10 -filePath $motion.g3 -dt $dt -factor [expr $g*0.5]
       } elseif {$motion == "se3001"} {
           timeSeries Path 10 -filePath $motion.g3 -dt $dt -factor [expr ($g*.1)/980.6]
       } elseif {$motion == "se305"} {
           timeSeries Path 10 -filePath $motion.g3 -dt $dt -factor [expr ($g*.5)/980.6]
       } elseif {$motion == "sin"} {
           timeSeries Sine 10 9 12 $sinPeriod -factor [expr $g*$LINfactor*$nfactor]
           set dt $dtSin
           set nPt $nPtSin
       } elseif {$motion == "sin01"} {
           timeSeries Sine 10 9 12 $sinPeriod -factor [expr $g*$LINfactor*$nfactor]
           set dt $dtSin
           set nPt $nPtSin
       } elseif {$motion == "CanogaNS"} {
           timeSeries Path 10 -filePath $motion.txt -dt $dt -factor [expr ($g/9806.0*$scaleFactor)]
       } elseif {$motion == "CanogaNS01"} {
           timeSeries Path 10 -filePath $motion.txt -dt $dt -factor [expr ($g/9806.0*$LINfactor*$nfactor)]
       } elseif {$motion == "LA3536PPhor"} {
           timeSeries Path 10 -filePath $motion.g3 -dt $dt -factor [expr ($LINfactor*$nfactor)]
       } else {
           timeSeries Path 10 -filePath $motion.g3 -dt $dt -factor [expr ($g*$LINfactor*$nfactor)]
#            timeSeries Path 10 -filePath $motion.g3 -dt $dt -factor [expr $g/980.6]; # ss01-40 in units cm/sec2

       }

# -----        # ## set secondary Rayleigh damping
if {$damp == "M" } {
   set numMode 30; #<-----
   set omegaI [expr pow([lindex $lambda 0],0.5)];
   set omegaJ [expr pow([lindex $lambda 2],0.5)];
   set alphaM [expr $pDampSecondary*(2*$omegaI*$omegaJ)/($omegaI+$omegaJ)];
   set betaK [expr 2.0*$pDampSecondary/$omegaI]
   rayleigh $alphaM  0.      $betaK       0.
}
puts "HI -1"
analyze 1 $dt

pattern UniformExcitation 1 1 -accel 10;

puts "HI -2"
analyze 1 $dt
puts "HI -3"

# ----- # ----- define lists for recorder ----- #
           set nodeList []
           set nodeListCLine1 []
               set colFloorList []

           set MbList []
           set one 1
           set two 2
           for {set floor $numFloor} {$floor >=2} {incr floor -1} {
               set beam  [lindex $beamSizes [expr $floor -2]]
               set prop [array get WSection $beam]
               set propList [split $prop]
               set Sx [expr [lindex $propList 8]*$in*$in*$in]
               lappend nodeListCLine1 1$floor$p05
               for {set col 1} {$col <= $numCline} {incr col 1} {
                   lappend nodeList $col$floor$p05
                   lappend MbList [expr $Sx*$Fyb]
                   #                lappend MbList 1.0
               }
           }

           puts "================================"
           set inodeList {1171}; # base node is not in the form of xy05
           set jnodeList {1205};
           for {set floor 2} {$floor < $numFloor} {incr floor 1} {
               lappend inodeList 1$floor$p05
               lappend jnodeList 1[expr $floor+1]$p05
           }
#            puts "$inodeList -- $jnodeList"

           set eleList []
           set floor1 1
           set floor2 2

           for {set col 1} {$col <= $numCline} {incr col 1} {
               lappend eleList $col$floor1$col$floor2
           }

           #            eval "print ele -$eleListHinge"
           #            puts ""
           #            puts $eleListHinge
           #            puts [llength $eleListHinge]

        # ---- output Mb of beam
               set Mbfile [open Mb.out w]
           puts $Mbfile $MbList
           close $Mbfile


       puts "TOTAL MASS: $totalMass"
       set totalMass 0
       foreach mass $nodeListMassValues {
           set totalMass [expr $totalMass+$mass]
       }
       puts "TOTAL MASS: $totalMass"

# =========================== set recorders =============================== #
       set recordAll "YES";  #<-----
       if {$recordAll == "YES"} {
           set cmd "recorder Node -file $dataDir/floorDispCLine1$damp$motion$colType$beamType.out  -time -node $nodeListCLine1 -dof 1 disp"; eval $cmd;
               set cmd "recorder Node -file $dataDir/floorVelCLine1$damp$motion$colType$beamType.out  -time -node $nodeListCLine1 -dof 1 vel"; eval $cmd;
               # set cmd "recorder Node -file $dataDir/floorDispCLine1Grav$damp$motion$colType$beamType.out  -time -node $nodeListCLine1Grav -dof 1 disp"; eval $cmd;

           # puts $cmd

           set cmd "recorder Drift -file $dataDir/floorDrift$damp$motion$colType$beamType.out  -time -iNode $inodeList -jNode $jnodeList  -dof 1 -perpDirn 2 "; eval $cmd;
           # puts $cmd

           set cmd "recorder Node -file $dataDir/NodeMassAccel$damp$motion$colType$beamType.out  -time -node $nodeListMass -dof 1 accel"; eval $cmd;
           set cmd "recorder Node -file $dataDir/NodeMassVel$damp$motion$colType$beamType.out  -time -node $nodeListMass -dof 1 vel"; eval $cmd;
           set cmd "recorder Node -file $dataDir/NodeMassUnbalance$damp$motion$colType$beamType.out  -time -node $nodeListMass -dof 1 unbalance"; eval $cmd;

           set cmd "recorder Node -file $dataDir/NodeMassUnbalanceIncInertia$damp$motion$colType$beamType.out  -time -node $nodeListMass -dof 1 unbalanceIncInertia"; eval $cmd;
           set cmd "recorder Node -file $dataDir/NodeMassRayleighForces$damp$motion$colType$beamType.out  -time -node $nodeListMass -dof 1 nodalRayleighForces"; eval $cmd;
           # puts $cmd

# ---- write matlab file for inertia force, fi ---------------#
           puts $matlab "load $dataDir/NodeMassAccel$damp$motion$colType$beamType.out"
           puts $matlab "a=NodeMassAccel$damp$motion$colType$beamType;"
           set myCounter 2;
           puts $matlab "fi=zeros(size(a,1),1);"
           foreach mass $nodeListMassValues {
               puts $matlab "fi(:,1)=fi(:,1) + $mass*a(:,$myCounter);"
               incr myCounter;
           }
# ------------------------------------------------------------#

           set cmd "recorder Node -file $dataDir/NodeReaction$damp$motion$colType$beamType.out  -time -node $nodeListReaction -dof 1 reaction"; eval $cmd;
           set cmd "recorder Node -file $dataDir/NodeReactionIncInertia$damp$motion$colType$beamType.out  -time -node $nodeListReaction -dof 1 reactionIncInertia"; eval $cmd;
           set cmd "recorder Node -file $dataDir/NodeRayleigh$damp$motion$colType$beamType.out  -time -node $nodeListReaction -dof 1 rayleighForces"; eval $cmd;
           # puts $cmd

           for {set floor1 1; set floor2 2} {$floor1 < $numFloor} {incr floor1 1; incr floor2 1} {
               for {set colLine 1} {$colLine <= [expr $numCline+1]} {incr colLine 1} {
                   lappend colFloorList $colLine$floor1$colLine$floor2
                       lappend colFloorGravList $colLine$floor1$colLine$floor2$pFrame
               }
           }
           # puts $colFloorList
           set cmd "recorder Element -file $dataDir/EleForcesFloors$damp$motion$colType$beamType.out  -time -ele $colFloorList forces"; eval $cmd; puts $cmd
               set cmd "recorder Element -file $dataDir/EleForcesFloorsGrav$damp$motion$colType$beamType.out  -time -ele $colFloorGravList forces"; eval $cmd; puts $cmd


               if {$damp == "equivFloorDampingCapped" || $damp == "equivFloorDamping"} {
           set cmd "recorder Element -file $dataDir/floorDamperForce$damp$motion$colType$beamType.out  -time -ele $FloorDamperList localForce"; eval $cmd;
               set cmd "recorder Element -file $dataDir/floorDamperDeformation$damp$motion$colType$beamType.out  -time -ele $FloorDamperList deformation"; eval $cmd;
               }
       }

       node 1 0 0 0
       fix  1 0 1 1
       pattern MultiSupport 100 {
           groundMotion 1  Series -accel 10
           imposedSupportMotion 1 1 1
       }

       set cmd "recorder Node -file $dataDir/motionAccel$damp$motion$colType$beamType.out  -time -node 1 -dof 1 accel"; eval $cmd;
       # set cmd "recorder Node -file $dataDir/motionAccel$damp$motion$colType$beamType.out -timeSeries 10 -time -node 1 -dof 1 accel; "; eval $cmd;


# -------- write matlab file for stiffness force, fs, and damping force, fd -------------#
       puts $matlab "load $dataDir/NodeReaction$damp$motion$colType$beamType.out"
       # fs = -reactions
       puts $matlab "b=-NodeReaction$damp$motion$colType$beamType;"
       set myCounter 2;
       puts $matlab "fs=zeros(size(a,1),1);"
       foreach nodeTag $nodeListReaction {
           puts $matlab "fs(:,1)=fs(:,1) + b(:,$myCounter);"
           incr myCounter;
       }

       puts $matlab "load $dataDir/motionAccel$damp$motion$colType$beamType.out"
#        puts $matlab "load se30.DATA;"
#        puts $matlab "se30(:,2) =  -((30.75*$g*1.0)/980.6)*se30(:,2);"
       puts $matlab "se30 = motionAccel$damp$motion$colType$beamType"
       puts $matlab "se30(:,2) =  -30.75*se30(:,2);";   #total mass per floor *a
       puts $matlab "fd = se30(:,2)-fs-fi;"
       puts $matlab "subplot(3,1,$plotCounter);";
       puts $matlab "plot(a(:,1),fs,a(:,1),fd(:,1),a(:,1),fi(:,1));"
       puts $matlab "axis(\[0 100 -2000 2000\])"
       puts $matlab "legend('fs','fd','fi')"
       puts $matlab "title('$damp')"
       set plotCounter [expr $plotCounter+1]
# --------------------------------------------------------------------------------------#
           set cmd "recorder EnvelopeNode -file $dataDir/floorDispEnv.out  -node $nodeList -dof 1 disp"; eval $cmd;

       if {$recordAll == "YES"} {
       set cmd "recorder EnvelopeNode -file $dataDir/floorAccEnv.out  -timeSeries 10 -node $nodeList -dof 1 accel"; eval $cmd;
           set cmd "recorder Node -file $dataDir/floorDisp$damp$motion$colType$beamType.out  -time -node $nodeList -dof 1 disp"; eval $cmd;
           set cmd "recorder Node -file $dataDir/floorVel$damp$motion$colType$beamType.out  -time -node $nodeList -dof 1 vel"; eval $cmd;
           set cmd "recorder Node -file $dataDir/floorAccel$damp$motion$colType$beamType.out  -time -node $nodeList -dof 1 accel"; eval $cmd;
           set cmd "recorder Node -file $dataDir/$motion.dat  -time -timeSeries 10 -node 11 -dof 1 accel"; eval $cmd;

           set cmd "recorder Node -file $dataDir/floorRotVel$damp$motion$colType$beamType.out  -time -node $nodeList -dof 3 vel"; eval $cmd;

           set cmd "recorder Node -file $dataDir/floorReaction$damp$motion$colType$beamType.out  -time -node $nodeList -dof 3 reaction"; eval $cmd;

           set cmd "recorder Node -file $dataDir/floorRayleigh$damp$motion$colType$beamType.out  -time -node $nodeList  -dof 3 rayleighForces"; eval $cmd;

           set cmd "recorder Node -file $dataDir/floorRayleigh1$damp$motion$colType$beamType.out  -time -node $nodeList  -dof 1 rayleighForces"; eval $cmd;


           set cmd "recorder Element -file $dataDir/colEleForces$damp$motion$colType$beamType.out  -time -ele $eleList forces"; eval $cmd;

               set cmd "recorder Element -file $dataDir/PZrigidEleForces$damp$motion$colType$beamType.out  -time -ele $PZrigidList forces"; eval $cmd;

               # add if having panel zone element
               set cmd "recorder Element -file $dataDir/PZhingeForce$damp$motion$colType$beamType.out  -time -ele $PZspringList force -dir 6"; eval $cmd;
               set cmd "recorder Element -file $dataDir/PZhingeDeformation$damp$motion$colType$beamType.out  -time -ele $PZspringList deformation"; eval $cmd;



           if {$beamType == "Force" || $beamType == "Displacement"} {
               set cmd "recorder Element -file $dataDir/hingeStressStrain$damp$motion$colType$beamType.out  -time -ele $eleListHinge section forceAndDeformation "; eval $cmd;
               set cmd "recorder Element -file $dataDir/hingeStrain$damp$motion$colType$beamType.out  -time -ele $eleListHinge section deformation "; eval $cmd;
               set cmd "recorder Element -file $dataDir/hingeDefo$damp$motion$colType$beamType.out  -time -ele $eleListHinge basicDeformation "; eval $cmd;
               set cmd "recorder Element -file $dataDir/hingePlasticDefo$damp$motion$colType$beamType.out  -time -ele $eleListHinge plasticDeformation "; eval $cmd;
           } elseif {$beamType == "ComponentElement" } {
               set cmd "recorder Element -file $dataDir/hingeStressStrain$damp$motion$colType$beamType.out  -time -ele $eleListHinge hingeDefoAndForce "; eval $cmd;
               set cmd "recorder Element -file $dataDir/hingeTangent$damp$motion$colType$beamType.out  -time -ele $eleListHinge hingeTangent "; eval $cmd;

               } else {
               set cmd "recorder Element -file $dataDir/hingeTangent$damp$motion$colType$beamType.out  -time -ele $eleListHinge material 1 tangent "; eval $cmd;
               set cmd "recorder Element -file $dataDir/hingeDampingForces$damp$motion$colType$beamType.out  -time -ele $eleListHinge dampingForces "; eval $cmd;
               set cmd "recorder Element -file $dataDir/hingeStressStrain$damp$motion$colType$beamType.out  -time -ele $eleListHinge defoANDforce "; eval $cmd;
               set cmd "recorder Element -file $dataDir/hingePlasticStrain$damp$motion$colType$beamType.out  -time -ele $eleListHinge plasticStrain "; eval $cmd;
           }

           if {$colType == "Force" || $colType == "Displacement"} {
               ;
           } else {
               set cmd "recorder Element -file $dataDir/colHingeStressStrain$damp$motion$colType$beamType.out  -time -ele $eleColListHinge hingeDefoAndForce "; eval $cmd; puts $cmd;
               set cmd "recorder Element -file $dataDir/colPlasticStrain$damp$motion$colType$beamType.out  -time -ele $eleColListHinge plasticStrain "; eval $cmd;
           set cmd "recorder Element -file $dataDir/colHingeTangent$damp$motion$colType$beamType.out  -time -ele $eleColListHinge hingeTangent "; eval $cmd; puts $cmd;
               }
       }

       set tFinal        [expr $dt*$nPt+5.0];        # maximum duration of ground-motion analysis

       # if {$damp == "M" } {
           # set lambda [eigen $numMode]
           # puts "LAMBDA $lambda"
               # puts " Period T1: [expr 2*$PI/sqrt([lindex $lambda 0])] sec [lindex $lambda 0]"
               # puts " Period T2: [expr 2*$PI/sqrt([lindex $lambda 1])] sec [lindex $lambda 1]"
               # puts " Period T3: [expr 2*$PI/sqrt([lindex $lambda 2])] sec [lindex $lambda 2]"
           # exit
           #
               # Damping  $pDamp
       # }


# ============ define solution algorithm =============#
       set ok 0
       set currentTime 0.0
       set countAnalyze 0; set countFail 0

       while {$ok == 0 && $currentTime < $tFinal} {
           incr countAnalyze 1


               # if {$ok != 0} {
               # test $testType [expr $Tol*1e1] 1000 0
               # algorithm NewtonHallM 0.2
               # set ok [analyze 1 [expr $dt/($n)]]
               # test $testType $Tol $numStep
               # #                algorithm NewtonLineSearch
               # algorithm NewtonHallM 0.1
           # }

               # if {$ok != 0} {
               # test $testType [expr $Tol*1e1] 1000 0
               # algorithm NewtonHallM 0.4
               # set ok [analyze 1 [expr $dt/($n)]]
               # test $testType $Tol $numStep
               # #                algorithm NewtonLineSearch
               # algorithm NewtonHallM 0.1
           # }

               # if {$ok != 0} {
               # test $testType [expr $Tol*1e1] 1000 0
               # algorithm NewtonHallM 0.6
               # set ok [analyze 1 [expr $dt/($n)]]
               # test $testType $Tol $numStep
               # #                algorithm NewtonLineSearch
               # algorithm NewtonHallM 0.1
           # }

               # if {$ok != 0} {
               # test $testType [expr $Tol*1e1] 1000 0
               # algorithm NewtonHallM 0.8
               # set ok [analyze 1 [expr $dt/($n)]]
               # test $testType $Tol $numStep
               # #                algorithm NewtonLineSearch
               # algorithm NewtonHallM 0.1
           # }

               set dt_analysis $dt;
               set dt_anal_min [expr $dt/(2.*50.)];
               set dt_anal_max [expr $dt];
               set nr_iter 25;
               set nr_iter_50 50;
               set nr_iter_100 100;
               set nr_analyse 50;
               set tol_0 1.0e-6;
               set tol_1 1.0e-5;
               set tol_2 1.0e-4;
               set tol_3 1.0e-3;
               set tol_4 1.0e-2;
               set main_test RelativeNormDispIncr;
               # set tol_0 1.0e-3;
               # set tol_1 1.0e-2;
               # set tol_2 5.0e-2;
               # set tol_3 1.0e-2;
               # set tol_4 1.0e-1;
               # set main_test RelativeNormUnbalance;
               set show_iter 2;

               set doOld "TRUE"
               if {$doOld == "TRUE"} {
                   while {$ok == 0 && $currentTime < $tFinal} {
                       incr countAnalyze 1

                               set ok [analyze 1 $dt_analysis $dt_anal_min $dt_anal_max $nr_iter];
                               set controlTime [getTime];
                               set krylovflag 1; # do not check again the 1st iteration unconverged Krylov

                               # ------- TOLERANCE_0 ----------------------------
                               if {$ok != 0} {
                                       puts "tol 0";
                                       set currentTolerance $tol_0;
                                       set currentdt [expr $dt_analysis/50.];
                                       source DynamicSolutionAlgorithmSubFile.tcl
                                       set controlTime [getTime];
                               }
                               # ------- TOLERANCE_1 ----------------------------
                               if {$ok != 0} {
                                       puts "tol 1";
                                       set currentTolerance $tol_1;
                                       set currentdt [expr $dt_analysis/50.];
                                       source DynamicSolutionAlgorithmSubFile.tcl
                                       set controlTime [getTime];
                               }
                               # ------- TOLERANCE_2 ----------------------------
                               if {$ok != 0} {
                                       puts "tol 2";
                                       set currentTolerance $tol_2;
                                       set currentdt [expr $dt_analysis/50.];
                                       source DynamicSolutionAlgorithmSubFile.tcl
                                       set controlTime [getTime];
                               }
                               # ------- TOLERANCE_3 ----------------------------
                               if {$ok != 0} {
                                       puts "tol 3";
                                       set currentTolerance $tol_3;
                                       set currentdt [expr $dt_analysis/50.];
                                       source DynamicSolutionAlgorithmSubFile.tcl
                                       set controlTime [getTime];
                               }

                               # ------- TOLERANCE_4 ----------------------------
                               if {$ok != 0} {
                                       puts "tol 4";
                                       set currentTolerance $tol_4;
                                       set currentdt [expr $dt_analysis/50.];
                                       source DynamicSolutionAlgorithmSubFile.tcl
                                       set controlTime [getTime];
                               }
                       set currentTime [getTime]
                       }
               } else {
#                    set nPts [expr int($tFinal/$dt)]
#                    puts "NPTS: $nPts"
#                    set ok [analyze $nPts $dt]
                   set ok 0
                   while {$ok == 0 && $currentTime < $tFinal} {
                       incr countAnalyze 1

                       set ok [analyze 1 $dt $dt_anal_min $dt_anal_max 100]
                       if {$ok != 0} {
                           analysis Transient
                           test $testType [expr $Tol*1e1] 1000 2
                           algorithm Newton
                           set ok [analyze 1 $dt]
                           test $testType $Tol $numStep
                           algorithm Newton
                       }
                       if {$ok != 0} {
                           test $testType [expr $Tol*1e1] 1000 2
                           algorithm ModifiedNewton -initial
                           set ok [analyze 1 [expr $dt/(2*$n)]]
                           test $testType $Tol $numStep
                           #                algorithm NewtonLineSearch
                           algorithm Newton
                       }


                       set currentTime [getTime]
                   }
               }
       }


       set tEnd1 [clock clicks -milliseconds]
       set timeFile [open $dataDir/$frame.time.$pid a]
       if {$ok == 0} {
           set ok SUCCESS
       } else {
           set ok FAILURE
       }

       remove recorders

       puts "\n Results"
       set a [open $dataDir/floorDispEnv.out r]
       set line [gets $a];
       set line [gets $a];
       set line [gets $a]
       set line [lindex [split $line " "] 0]

       puts "$ok DURATION: [expr ($tEnd1-$tStart1)/1000] zeta: $pDamp MAX: $line $motion damp: $damp colType: $colType beamType: $beamType numMode: $numMode  analyzeSteps: $countAnalyze failSteps: $countFail [expr 100*$countFail/($countAnalyze*1.0)]\%"

       puts $timeFile "$ok [getTime] DURATION: [expr ($tEnd1-$tStart1)/1000] zeta: $pDamp MAX: $line $motion damp: $damp colType: $colType beamType: $beamType numMode: $numMode  analyzeSteps: $countAnalyze failSteps: $countFail [expr 100*$countFail/($countAnalyze*1.0)]\% [nodeDisp 12105]"
       close $timeFile
       close $a
       wipe;
      }
    }
  }
  incr counter 1
}

close $matlab
