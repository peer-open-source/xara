# EXAMPLE 0: Sectional response for combined in-plane/out-of-plane loading
#
# @file   MomentCurvature_rotated.tcl
#
# @author: Francesco Vanin <francesco.vanin@epfl.ch>
#
# @date creation: 15 Feb 2019
#
# @date modification: 15 Feb 2019
#
# @brief: The model consists in a zero-length section equipped with the analytically integrated 
# NoTensionSection3d model, on which a moment-curvature analysis is performed. Crushing is included
# in the formulation. The orientation of the moment acting on the section is controlled by the 
# variable momentAngle (line 50). The moment orientation (not the curvature, in the nonlinear 
# range) is maintained constant in the analysis. 
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
wipe
model BasicBuilder -ndm 3 -ndf 6

#VARIABLES ----------------------------------------------------------------------------------------- 
# Dimensions: m, s, N, Pa
# section geometry
set L  		   	 1.0;     # length (m)
set t            0.2;     # thickness (m)
set oopSlices     20;     # number of slices considered in the out-of-plane direction

# material properties
set E          6.0e9;     # Young's modulus (MPa)
set G          2.0e9;     # Shear modulus (MPa) - controls torsion
set fc         5.0e6;     # Compressive strength (MPa)

# loading options
set ALR         0.25;     # axial load ratio
set  pi         3.14;
set momentAngle  30.;     # moment orientation (°deg)

set secTag        10;

#NODES ---------------------------------------------------------------------------------------------
#     tag     X     Y     Z  
node 	1 	0.0   0.0   0.0
node 	2 	0.0   0.0   0.0

# Fix all degrees of freedom except axial and bending
fix 	1   1 1 1 1 1 1
fix 	2   0 1 1 1 0 0


#ZERO-LENGTH SECTION--------------------------------------------------------------------------------
section NoTensionSection3d $secTag $E $G $L $t -1.0  $fc  $oopSlices

# --- alternative: fiber section with the same material model ---
#set  matTag  1;
#set numfibersIP       50;
#set numfibersOOP       5;
#uniaxialMaterial CompressionDamage1d   $matTag   $E   $fc
#section   Fiber  $secTag  { patch rect $matTag  $numfibersIP $numfibersOOP  [expr -$L/2.0] [expr -$t/2.0]  [expr $L/2.0] [expr $t/2.0]  } 
# ---

# Define element
set theta    [expr $momentAngle/180.*$pi]
#                         tag ndI ndJ   section tag
element zeroLengthSection   1   1   2  $secTag       -orient  1.0 0.0 0.0    0.0 [expr cos($theta)] [expr sin($theta)]


#RECORDERS ----------------------------------------------------------------------------------------- 
recorder Node -file sectionResponse.out -time -node 2 -dof 6 1 disp


# -------------------------------------------------------------------------------------------------- 
#ANALYSIS ------------------------------------------------------------------------------------------ 
# -------------------------------------------------------------------------------------------------- 

#AXIAL LOAD----------------------------------------------------------------------------------------- 
pattern Plain 1 "Constant" {
	load    2      [expr -1.0*$fc*$L*$t*$ALR]  0.0  0.0   0.0  0.0  0.0
}

# Define analysis parameters
integrator LoadControl 0.0
system SparseGeneral -piv;	
test NormUnbalance 1.0e-1 100 5 
numberer Plain
constraints Plain
algorithm Newton 
analysis Static

analyze 1

#CYCLIC CURVATURE ---------------------------------------------------------------------------------- 
pattern Plain 2 "Linear" {
	load 2   0.0  0.0  0.0  0.0  0.0  1.0
}

set numIncr   300;	             # Number of analysis increments to reach max curvature
set maxK     0.03;
set dK [expr $maxK/$numIncr];    # maxCurvature



set numIncr [expr int($numIncr/3)];

integrator DisplacementControl 2 6 [expr +1.0*$dK ]
analyze [expr $numIncr]
		
integrator DisplacementControl 2 6 [expr -1.0*$dK ]
analyze [expr $numIncr*2]
		
integrator DisplacementControl 2 6 [expr 1.0*$dK ]
analyze [expr $numIncr*3]

integrator DisplacementControl 2 6 [expr -1.0*$dK ]
analyze [expr $numIncr*4]

integrator DisplacementControl 2 6 [expr 1.0*$dK]
analyze [expr $numIncr*5]

integrator DisplacementControl 2 6 [expr -1.0*$dK]
analyze [expr $numIncr*6]

integrator DisplacementControl 2 6 [expr 1.0*$dK]
analyze [expr $numIncr*3]
	
puts "Analysis completed."
wipe
exit