    # ComponentBeamWSection2d $colLine$floor1$colLine$floor2                  $colLine$floor1$p7 909$colLine$floor1 $theSectionlow $E $Fyc $Ry 1 -matType Bilin -Com_Type other-than-RBS
    # ComponentBeamWSection2d 909$colLine$floor1$colLine$floor2 909$colLine$floor1 $colLine$floor2$p6 $theSection $E $Fyc $Ry 1 -matType Bilin -Com_Type other-than-RBS
    ComponentColWSection2d $colLine$floor1$colLine$floor2                  $colLine$floor1$p7 909$colLine$floor1 $theSectionlow $E $Fyc $Ry [set $PgPyb] 1 -matType Bilin -nFactor $nFactorElem
    ComponentColWSection2d 909$colLine$floor1$colLine$floor2 909$colLine$floor1 $colLine$floor2$p6 $theSection $E $Fyc $Ry [set $PgPyt] 1 -matType Bilin -nFactor $nFactorElem
    ComponentColWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection $E $Fyc $Ry [set $PgPy] 1 -matType Bilin -nFactor $nFactorElem
    # ComponentBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection $E $Fyc $Ry 1 -matType Bilin -Com_Type other-than-RBS
    # ComponentBeamWSection2d $colLine$floor1$colLine$floor2                  $colLine$floor1$p7 909$colLine$floor1 $theSectionlow $E $Fyc $Ry 1 -matType Bilin -metric -Com_Type other-than-RBS
    # ComponentBeamWSection2d 909$colLine$floor1$colLine$floor2 909$colLine$floor1 $colLine$floor2$p6 $theSection $E $Fyc $Ry 1 -matType Bilin -metric -Com_Type other-than-RBS
    ComponentColWSection2d $colLine$floor1$colLine$floor2                  $colLine$floor1$p7 909$colLine$floor1 $theSectionlow $E $Fyc $Ry [set $PgPyb] 1 -matType Bilin -nFactor $nFactorElem
    ComponentColWSection2d 909$colLine$floor1$colLine$floor2 909$colLine$floor1 $colLine$floor2$p6 $theSection $E $Fyc $Ry [set $PgPyt] 1 -matType Bilin -nFactor $nFactorElem
    # ComponentBeamWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection $E $Fyc $Ry 1 -matType Bilin -metric -Com_Type other-than-RBS
    ComponentColWSection2d $colLine$floor1$colLine$floor2 $colLine$floor1$p7 $colLine$floor2$p6 $theSection $E $Fyc $Ry [set $PgPy] 1 -matType Bilin -nFactor $nFactorElem
  ComponentBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection $E $Fyb $Ry 2 -matType Bilin -Com_Type RBS -nFactor $nFactorElem
  # ComponentBeamWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection $E $Fyb 2 -matType MultiLinear -metric
  BeamWithConcentratedHingesWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection [expr $EIfactor*$E] [expr $EIfactor*$Fyb] 2 -nFactor $nFactorElem -doRayleigh $doRayleigh -matType MultiLinear -metric
  BeamWithConcentratedHingesWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection $E $Fyb 2 -nFactor $nFactorElem -doRayleigh $doRayleigh -matType Bilin02 -metric
  BeamWithConcentratedHingesWSection2d $colLine1$floor$colLine2$floor $colLine1$floor$p1 $colLine2$floor$p4 $theSection $E $Fyb 2 -nFactor $nFactorElem -doRayleigh $doRayleigh -matType MultiLinear -metric
