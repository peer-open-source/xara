# Plane 3bar Truss Example
#
# planar 3 bar system, all bars same A & E, unit load P
#

set A 10.0
set E 3000.
set L 200.0
set alpha 30.0
set P 200.0

set sigmaYP 60.0

set pi [expr 2.0*asin(1.0)]
set alphaRad [expr $alpha*$pi/180.]
set cosA [expr cos($alphaRad)]
set sinA [expr sin($alphaRad)]

set dX [expr $L*tan($alphaRad)]


#EXACT
# Exact per Popov
set PA [expr ($sigmaYP*$A) * (1.0+2*$cosA*$cosA*$cosA)]
set PB [expr ($sigmaYP*$A) * (1.0+2*$cosA)]

# create the finite element model for nonlinear case

wipe

model Basic -ndm 2 -ndf 2

node 1    0.0          0.0
node 2    $dX          0.0
node 3 [expr 2.0*$dX]  0.0
node 4    $dX         -$L     

fix 1 1 1
fix 2 1 1
fix 3 1 1

uniaxialMaterial ElasticPP 1 $E [expr $sigmaYP/$E]
element Truss 1 1 4 $A 1
element Truss 2 2 4 $A 1
element Truss 3 3 4 $A 1

eval "timeSeries Path 1 -dt 1.0 -values {0.0 $PA $PB $PB}"
pattern Plain 1 1 {
    load 4 0. -1.0
}

print -json
