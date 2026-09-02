# EXAMPLE 2d: One gable element subjected to out-of-plane loading
#
# @file   OneElement_OOP_case4.tcl
#
# @author: Francesco Vanin <francesco.vanin@epfl.ch>
#
# @date creation: 15 Feb 2019
#
# @date modification: 15 Feb 2019
#
# @brief: this input file produces the response of a single cantilever gable subjected to out-of-plane 
# monotonic loading, including P-Delta effects. The gable is loaded with triagular load pattern both in  
# vertical and horizontal directions. Results can be compared to a the lateral load multiplier 
# calculated through limit analysis. 
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
set  hTot      3.00;            # total height (m)
set  L         6.00;            # length (m)
set  t         0.30;            # thickness (m)

# Material properties
set E          2000.0e+06;      # Young's modulus (Pa)
set G          700.0e+06;       # Shear modulus (Pa)
set fc         3.0e+06;         # Compressive strength of masonry (Pa)
set c          0.150e+06;       # Equivalent cohesion (Pa)
set Gc         3.0;             # Parameter Gc, defines softening of the response in the pre-peak branch (-)
set mu0        0.4;             # Friction coeffcient at peak force (-)
set muR        0.4;             # Residual friction coeffcient (-) 
set dropDrift  0.008;           # drift ratio at 20% force drop, for shear failure 
set rho        2200.;           # material density (kg/m3) 

set g          9.81;            # gravity acceleration (m/s2)


#NODES ---------------------------------------------------------------------------------------------
#       tag     X         Y         Z  
node      1     0.000     0.000     0.000 
node      2     0.000     0.000     [expr $hTot/2.]
node      3     0.000     0.000     $hTot

#CONSTRAINTS ---------------------------------------------------------------------------------------
fix       1    1   1   1   1   1   1 


#MACROELEMENTS ------------------------------------------------------------------------------------- 
# --- standard implementation ---
#                         eTag    nodesI,J,E    axis vector     out-of-plane vector    
element Macroelement3d    1       1  3  2       0.0  0.0  1.0   1.0  0.0  0.0        -gable  $hTot $L $t $E $G $fc $mu0 $c $Gc $dropDrift $muR  -pDelta -density $rho 
# ---


#RECORDERS ----------------------------------------------------------------------------------------- 
recorder Node -file dTop.out   -time -node  3  -dof 1 3  disp


# -------------------------------------------------------------------------------------------------- 
#ANALYSIS ------------------------------------------------------------------------------------------ 
# -------------------------------------------------------------------------------------------------- 

#AXIAL LOAD----------------------------------------------------------------------------------------- 
set mTot     [expr $rho*$t*$L*$hTot]
pattern Plain 10 "Linear" {     
	 eleLoad -ele 1  -type -selfWeight 0.0 0.0 -$g
}
	
# Define analysis parameters
system SparseGEN;
numberer Plain;
constraints Transformation; 
integrator LoadControl 1
test NormUnbalance 1.0e-1  50 1
algorithm Newton
analysis Static
analyze 1
loadConst -time 0.0
record

#IMPOSED LATERAL DISPLACEMENT, x direction---------------------------------------------------------- 	
pattern Plain 1 "Linear" {
	eleLoad -ele 1  -type -selfWeight 1.0 0.0 0.0  
}

wipeAnalysis;
test NormUnbalance 1.0e-1 50  0;
algorithm Newton;
constraints Transformation;
numberer Plain;
integrator LoadControl 1;
system SparseGEN;
analysis Static;
	
set initialDispl [expr 0.95*$t*1.5]
set incr 0.1e-3;   
set nSteps  [expr int(abs($initialDispl/$incr))] 

set controlled_node   3
set controlled_dof    1

test NormUnbalance 1.0e-1  50  2
algorithm Newton 
integrator    DisplacementControl     $controlled_node      $controlled_dof     [expr $incr] 

set ok [analyze $nSteps]


puts "Analysis completed."
wipe
exit