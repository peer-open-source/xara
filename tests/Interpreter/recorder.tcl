# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# Portal Frame Example 3.2
# ------------------------
#  Nonlinear pushover analysis using Portal Frame Example 1
#  as starting point
# 
# Units: kips, in, sec
#
# Written: GLF/MHS/fmk
# Date: January 2001

# ----------------------------------------------------
# Start of Model Generation & Initial Gravity Analysis
# ----------------------------------------------------

set width    360
set height   144

model basic -ndm 2 -ndf 3
# Create nodes
#    tag        X       Y 
node  1       0.0     0.0 
node  2    $width     0.0 
node  3       0.0 $height
node  4    $width $height

# Fix supports at base of columns
#    tag   DX   DY   RZ
fix   1     1    1    1
fix   2     1    1    1


# Define materials for nonlinear columns
# ------------------------------------------
# CONCRETE                  tag   f'c        ec0   f'cu        ecu
# Core concrete (confined)
uniaxialMaterial Concrete01  1  -6.0  -0.004   -5.0     -0.014

# Cover concrete (unconfined)
uniaxialMaterial Concrete01  2  -5.0   -0.002   0.0     -0.006

# STEEL
# Reinforcing steel 
pset fy 60.0;      # Yield stress
pset E 30000.0;    # Young's modulus

#                        tag  fy E0  b
uniaxialMaterial Steel01  3  $fy $E 0.01

# Define cross-section for nonlinear columns
# ------------------------------------------

# set some parameters
set colWidth 15
set colDepth 24 

set cover  1.5
set As    0.60;     # area of no. 7 bars

# some variables derived from the parameters
set y1 [expr $colDepth/2.0]
set z1 [expr $colWidth/2.0]

section Fiber 1 {
    # Create the concrete core fibers
    patch rect 1 10 1 [expr $cover-$y1] [expr $cover-$z1] [expr $y1-$cover] [expr $z1-$cover]

    # Create the concrete cover fibers (top, bottom, left, right)
    patch rect 2 10 1  [expr -$y1] [expr $z1-$cover] $y1 $z1
    patch rect 2 10 1  [expr -$y1] [expr -$z1] $y1 [expr $cover-$z1]
    patch rect 2  2 1  [expr -$y1] [expr $cover-$z1] [expr $cover-$y1] [expr $z1-$cover]
    patch rect 2  2 1  [expr $y1-$cover] [expr $cover-$z1] $y1 [expr $z1-$cover]

    # Create the reinforcing fibers (left, middle, right)
    layer straight 3 3 $As [expr $y1-$cover] [expr $z1-$cover] [expr $y1-$cover] [expr $cover-$z1]
    layer straight 3 2 $As 0.0 [expr $z1-$cover] 0.0 [expr $cover-$z1]
    layer straight 3 3 $As [expr $cover-$y1] [expr $z1-$cover] [expr $cover-$y1] [expr $cover-$z1]
}    


# Define column elements
# ----------------------

# Geometry of column elements
#                tag 
geomTransf Corotational 1  

# Number of integration points along length of element
set np 5

# Create the coulumns using beam-column elements
#               e       tag ndI ndJ nsecs secID transfTag
element ForceBeamColumn  1   1   3   $np    1       1 
element ForceBeamColumn  2   2   4   $np    1       1 

# Define beam element
# -----------------------------

# Geometry of column elements
#                tag 
geomTransf Linear 2  

# Create the beam element
#                          tag ndI ndJ     A       E    Iz   transfTag
element ElasticBeamColumn   3   3   4    360    4030  8640    2

# Define gravity loads
# --------------------

# Set a parameter for the axial load
set P 180;                # 10% of axial capacity of columns

# Create a Plain load pattern with a Linear TimeSeries
# driving point loads at nodes 3 and 4
pattern Plain 1 Linear {
        #    nd    FX          FY  MZ 
        load  3   0.0  [expr -$P] 0.0
        load  4   0.0  [expr -$P] 0.0
}
# initialize in case we need to do an initial stiffness iteration
initialize

#
# Start of analysis generation
#
# Configure the system of equations, a sparse solver with partial pivoting
system ProfileSPD

# Configure the constraint handler, the transformation method
constraints Transformation

# Configure the convergence test, the norm of the residual with a tolerance of 
# 1e-12 and a max number of iterations of 10
test NormDispIncr 1.0e-12  10 3

# Configure the integration scheme, the LoadControl scheme using steps of 0.1 
integrator LoadControl 0.1

# Configure the analysis
analysis Static

#
# Create Recorders
#
recorder Node    -xml  out/nodeGravity.xml -time -node 3 4 -dof 1 2 3 disp
recorder Node    -txt  out/nodeGravity.txt -time -node 3 4 -dof 1 2 3 disp
recorder Element -file out/elemGravity.out -ele 1 section  forces

#
# Finally perform the analysis
#
analyze 10
puts "Gravity load analysis completed";

# Set the gravity loads to be constant & reset the time in the domain
loadConst -time 0.0

# ----------------------------------------------------
# End of Model Generation & Initial Gravity Analysis
# ----------------------------------------------------


# ----------------------------------------------------
# Impose lateral loads
# ----------------------------------------------------

# Set some parameters
set H 10.0;                # Reference lateral load

# Set lateral load pattern with a Linear TimeSeries
pattern Plain 2 "Linear" {

        # Create nodal loads at nodes 3 & 4
        #    nd    FX  FY  MZ 
        load 3 $H 0.0 0.0 
        load 4 $H 0.0 0.0 
}

# End of additional modelling for lateral loads
# ----------------------------------------------------



# ----------------------------------------------------
# Configure push over analysis
# ----------------------------------------------------

# Change the integration scheme to be displacement control

set dU 0.1;                 # Displacement increment
#                             node dof init Jd min max
integrator DisplacementControl  3   1   $dU  1 $dU $dU


#
# Define recorders
#
# Stop the old recorders by destroying them
remove recorders

recorder Node  -file out/node32.out -time -node 3 4 -dof 1 2 3 disp
recorder Node  -txt  out/node32.txt -time -node 3 4 -dof 1 2 3 disp
recorder Node  -xml  out/node32.xml -time -node 3 4 -dof 1 2 3 disp
recorder Node  -csv  out/node32.csv -time -node 3 4 -dof 1 2 3 disp

recorder EnvelopeElement -file out/ele32.out -time -ele 1 2 localForce
recorder EnvelopeElement -txt  out/ele32.txt -time -ele 1 2 localForce
recorder EnvelopeElement -csv  out/ele32.csv -time -ele 1 2 localForce
recorder EnvelopeElement -xml  out/ele32.xml -time -ele 1 2 localForce

#
# Finally perform the analysis
#

# Set some parameters
set maxU 15.0;                # Max displacement
set numSteps [expr int($maxU/$dU)]

# Perform the analysis
set ok [analyze $numSteps]
if {$ok != 0} {
    set currentDisp [nodeDisp 3 1]
    set ok 0
    while {$ok == 0 && $currentDisp < $maxU} {

        set ok [analyze 1]

        # if the analysis fails try initial tangent iteration
        if {$ok != 0} {
            puts "... Newton failed, trying an initial stiffness"
            test NormUnbalance 1.0  1000 5
            algorithm ModifiedNewton -initial
            set ok [analyze 1]
            if {$ok == 0} {puts "... that worked, back to regular newton"}
            test NormDispIncr 1.0e-12  10 
            algorithm Newton
        }
        set currentDisp [nodeDisp 3 1]
    }
}

puts "";
if {$ok == 0} {
    puts "Pushover analysis completed SUCCESSFULLY";
} else {
    puts "Pushover analysis FAILED";    
}

print -file out/print_node_3.out node 3

print -file out/print.out

