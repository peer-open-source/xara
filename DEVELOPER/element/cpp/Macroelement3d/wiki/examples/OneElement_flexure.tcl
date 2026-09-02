# EXAMPLE 1: One element subjected to cyclic in-plane shear
#
# @file   OneElement_flexure.tcl
#
# @author: Francesco Vanin <francesco.vanin@epfl.ch>
#
# @date creation: 15 Feb 2019
#
# @date modification: 15 Feb 2019
#
# @brief: this input file produces the response of a single cantilever wall subjected to in-plane cyclic 
# loading. Several macroelement definitions are available (lines 72-94): comment/uncomment the desired 
# definition (-pier/-tremuri/-fiberSection) to switch from one to another.
#
# @section LICENSE
#
# Copyright (©) 2013-2018 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (EESD - Earthquake Engineering and Structural Dynamics Laboratory)
# More information can be found on our website <https://eesd.epfl.ch/>.
#
# All files are distributed under the GNU Lesser General Public License.
#
# Detailed information can be found in <http://www.gnu.org/licenses/>.
#
# --------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------- 
#MODEL --------------------------------------------------------------------------------------------- 
# -------------------------------------------------------------------------------------------------- 

# Create a model with three-dimensions (-ndm) and 6 DOF/node (-ndf)
wipe;
model basic -ndm 3 -ndf 6;     

#VARIABLES ----------------------------------------------------------------------------------------- 
# Dimensions: m, s, N, Pa
# wall geometry
set  h     0.90;                # total height (m)
set  L     0.90;                # length (m)
set  t     0.20;                # thickness (m)

# Material properties
set E          4200.0e+06;      # Young's modulus (Pa)
set G          1260.0e+06;      # Shear modulus (Pa)
set fc         10.0e+06;        # Compressive strength of masonry (Pa)
set c          0.050e+06;       # Equivalent cohesion (Pa)
set Gc         6.0;             # Parameter Gc, defines softening of the response in the pre-peak branch (-)
set mu0        0.47;            # Friction coefficient at peak force (-)
set muR        0.16;            # Residual friction coefficient (-) - only for "standard" implementation
set dropDrift  0.0087;          # drift ratio at 20% force drop, for shear failure  - only for "standard" implementation
set beta       0.4;             # Parameter beta, defines the post peak response (-) - only for "tremuri" implementation 

# drift model, if applied (constant drift capacity for all axial load ratios) 
set driftF     0.04;            # drift at loss of lateral force capacity, flexural failure (-)
set driftS     0.02;            # drift at loss of lateral force capacity, shear failure (-)

#NODES ---------------------------------------------------------------------------------------------
#       tag     X         Y         Z  
node      1     0.000     0.000     0.000 
node      2     0.000     0.000     [expr $h/2.]
node      3     0.000     0.000     $h


#MACROELEMENTS --------------------------------------------------------- 
# --- standard implementation ---
#                         eTag    nodesI,J,E    axis vector     out-of-plane vector    
element Macroelement3d    1       1  3  2       0.0  0.0  1.0   -1.0  0.0  0.0        -pier  $h $L $t $E $G $fc $mu0 $c $Gc $dropDrift $muR 
# ---

# --- alternative1: standard implementation with drift model ---
#                         eTag    nodesI,J,E    axis vector     out-of-plane vector    
#element Macroelement3d    1       1  3  2       0.0  0.0  1.0   -1.0  0.0  0.0        -pier  $h $L $t $E $G $fc $mu0 $c $Gc $dropDrift $muR  -driftFlexure $driftF 1.0 -driftShear $driftS 0.0
# ---

# --- alternative2: tremuri-equivalent element ---
#                         eTag    nodesI,J,E    axis vector     out-of-plane vector
#element Macroelement3d    1       1  3  2       0.0  0.0  1.0   -1.0  0.0  0.0        -tremuri $h $L $t $E $G $fc $mu0 $c $Gc $beta  
# ---

# --- alternative3: fiber section ---
#set numfibersIP       20;
#set numfibersOOP       5;
#set fiberSect        100;
#set tagMaterial       10;

#uniaxialMaterial   CompressionDamage1d    $tagMaterial   $E   $fc
#section  Fiber     $fiberSect   { patch rect $tagMaterial $numfibersIP $numfibersOOP  [expr -$L/2.0] [expr -$t/2.0]  [expr $L/2.0] [expr $t/2.0]  } 
#  #                       eTag    nodesI,J,E    axis vector     out-of-plane vector   
#element Macroelement3d    1       1  3  2       0.0  0.0  1.0   -1.0  0.0  0.0       -fiberSection  $fiberSect $fiberSect $fiberSect  $h $L $t $E $G $fc $mu0 $c $Gc $beta
# ---
puts "Elements defined." 


#CONSTRAINTS --------------------------------------------------------------------------------------- 
fix       1    1   1   1   1   1   1 


#RECORDERS ----------------------------------------------------------------------------------------- 
recorder Node    -file AllDispl.out      -time -node 3        -dof 2 3         disp
recorder Node    -file AllReactions.out  -time -nodeRange 1 3 -dof 1 2 3 4 5 6 reaction 
recorder Element -file AllDrifts.out     -time -ele 1  drift 
recorder Element -file Alpha.out         -time -ele 1  shear state 


# -------------------------------------------------------------------------------------------------- 
#ANALYSIS ------------------------------------------------------------------------------------------ 
# -------------------------------------------------------------------------------------------------- 

#AXIAL LOAD----------------------------------------------------------------------------------------- 
pattern Plain 1 "Linear" {
	# axial load ratio 15%
	load      3   0.0  0.0  -270.0e3  0.0   0.0   0.0 
}

# Define analysis parameters
system BandGeneral; 
numberer Plain;
constraints Transformation;	
integrator LoadControl 1;
test NormUnbalance 1.0e1  50 0;
algorithm Newton;
analysis Static;
analyze 1;
puts "Vertical load applied."
	
loadConst -time 0.0;
record;


#CYCLIC PUSHOVER, y direction----------------------------------------------------------------------- 
set controlled_node   3
set controlled_dof    2

pattern Plain 10 "Linear" {
    load      3    0.0  1.0e3  0.0   0.0   0.0   0.0
}

wipeAnalysis;
test NormUnbalance 1.0e1  50  0
algorithm Newton;
constraints Transformation;
numberer Plain;
integrator LoadControl 1;
system BandGeneral;
analysis Static;

# Loading history (mm) --------------------------------------------------------------------------- 	
set dispHistory [list  2 -2 4 -4 6 -6 8 -8 12 -12 16 -16 20 -20 30 -30 0]
	
	
# displacement increment (m)
set incr 0.10e-3;  

set node2 0.0
set dispOld 0.0

foreach dMax $dispHistory {
	set dMax [expr $dMax/1000.0]
	
	if {$dMax > $dispOld} {
		# positive loading
		while {$node2<$dMax} {	
			algorithm Newton 
			integrator    DisplacementControl     $controlled_node      $controlled_dof     [expr 1.0*$incr] 
			set ok [analyze 1]
		
			if {$ok != 0} {
				set stringToPlot ""
				append stringToPlot "Newton-Raphson failed at displacement " $node2;
				puts $stringToPlot;

				algorithm Newton -initial
				integrator    DisplacementControl     $controlled_node      $controlled_dof     [expr 1.0*$incr] 
				set ok [analyze 1]
				if {$ok != 0} {
					set stringToPlot ""
					append stringToPlot "Modified Newton-Raphson failed at displacement " $node2;
					puts $stringToPlot;
					break;
				}
			}
			set node2 [nodeDisp $controlled_node  $controlled_dof ]
		}
	
	} else {
	
		#negative loading
		set dMin [expr $dMax]
		while {$node2 > $dMin} {
			algorithm Newton 
			integrator    DisplacementControl     $controlled_node      $controlled_dof     [expr -1.0*$incr]  
			set ok [analyze 1]
		
			if {$ok != 0} {
				set stringToPlot ""
				append stringToPlot "Newton-Raphson failed at displacement " $node2;
				puts $stringToPlot;
						
				algorithm Newton -initial
				integrator    DisplacementControl     $controlled_node      $controlled_dof     [expr -1.0*$incr*0.1] 
				set ok [analyze 10]
				if {$ok != 0} {
					set stringToPlot ""
					append stringToPlot "Modified Newton-Raphson failed at displacement " $node2;
					puts $stringToPlot;
					break;
				}
			}
			set node2 [nodeDisp $controlled_node      $controlled_dof ]
		}
	}
	set dispOld $node2
}


puts "Analysis completed."
wipe
exit