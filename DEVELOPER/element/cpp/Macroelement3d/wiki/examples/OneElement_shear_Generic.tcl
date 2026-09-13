# EXAMPLE 1b: One element subjected to cyclic in-plane shear, generic uniaxial shear model
#
# @file   OneElement_shear_genericUniaxial.tcl
#
# @author: Francesco Vanin <francesco.vanin@epfl.ch>
#
# @date creation: 15 Feb 2019
#
# @date modification: 15 Feb 2019
#
# @brief: this input file produces the response of a single wall subjected to in-plane cyclic 
# loading. Base DOF's and top rotations are restrained. With the imposed axial load ratio the 
# element fails in shear. The shear response is defined by a generic uniaxial material model,
# uncoupled though from the axial response.
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
model basic -ndm 3 -ndf 6

#VARIABLES ----------------------------------------------------------------------------------------- 
set  h     0.90;                # total height (m)
set  L     0.90;                # length (m)
set  t     0.20;                # thickness (m)

# Material properties
set E          3000.0e+06;      # Young's modulus (Pa)
set G          1000.0e+06;      # Shear modulus (Pa)
set fc         5.0e+06;         # Compressive strength of masonry (Pa)
set c          0.30e+06;        # Equivalent cohesion (Pa)
set Gc         3.0;             # Parameter Gc, defines softening of the response in the pre-peak branch (-)
set mu         0.40;            # Friction coeffcient at peak force (-)
set muR        0.25;            # Residual friction coeffcient (-)
set dropDrift  0.004;           # drift ratio at 20% force drop, for shear failure 

set g         9.81; 
set pi   3.1415927; 



#NODES --------------------------------------------------------------------------------------------- 
#       tag     X         Y         Z   
node      1     0.000     0.000     0.000 
node      2     0.000     0.000     [expr $h/2.]
node      3     0.000     0.000     $h



#MACROELEMENTS ------------------------------------------------------------------------------------
set sec 100;
uniaxialMaterial CompressionDamage1d  10  $E  $fc
section Fiber $sec -GJ [expr $G*$L*$t*$t*$t/4.]  {
     patch rect 10  51 21   [expr -$L/2.] [expr -$t/2.] [expr $L/2.] [expr $t/2.]  }
	 
set shearIP  101;
set shearOOP 102;
	 
# translation of IKM model parameter for an equivalent shear response 	 
set N    -150.0e3;	 
set Vmax [expr $c*$L*$t - $mu*$N];     #shear force capacity, Mohr-Coulomb criterion
set Vy   [expr -$mu*$N];	           #shear force at cracking limit state       

set K0	    [expr 5./6.*$G*$L*$t/$h];  #elastic stiffness
set delta_max     [expr (-$muR*$N+(1.+$Gc)*($Vmax+$muR*$N))/$K0/$h]; # drift at peak force

set theta_p_Plus  [expr $h*$delta_max-$Vy/$K0];    #pre-capping rotation for positive loading direction (often noted as plastic rotation capacity)
set theta_p_Neg	  $theta_p_Plus;                   #pre-capping rotation for negative loading direction (often noted as plastic rotation capacity) (must be defined as a positive value)

puts $theta_p_Plus

set theta_pc_Plus [expr 5.*$h*($dropDrift-$delta_max)];    #post-capping rotation for positive loading direction
set theta_pc_Neg  $theta_pc_Plus;                          #post-capping rotation for negative loading direction (must be defined as a positive value)

puts $theta_pc_Plus

set theta_u_Plus  [expr 4.*$theta_pc_Plus];  #ultimate rotation capacity for positive loading direction
set theta_u_Neg	  $theta_u_Plus ;            #ultimate rotation capacity for negative loading direction (must be defined as a positive value)

set as_Plus	[expr ($Vmax-$Vy)/($theta_p_Plus*$K0)];    #strain hardening ratio for positive loading direction
set as_Neg	$as_Plus;                                  #strain hardening ratio for negative loading direction

set Lamda_S	0.0;     #Cyclic deterioration parameter for strength deterioration [E_t=Lamda_S*M_y, see Lignos and Krawinkler (2011); set Lamda_S = 0 to disable this mode of deterioration]
set Lamda_C	1.0;     #Cyclic deterioration parameter for post-capping strength deterioration [E_t=Lamda_C*M_y, see Lignos and Krawinkler (2011); set Lamda_C = 0 to disable this mode of deterioration]
set Lamda_A	0.0;     #Cyclic deterioration parameter for accelerated reloading stiffness deterioration [E_t=Lamda_A*M_y, see Lignos and Krawinkler (2011); set Lamda_A = 0 to disable this mode of deterioration]
set Lamda_K	0.0;     #Cyclic deterioration parameter for unloading stiffness deterioration [E_t=Lamda_K*M_y, see Lignos and Krawinkler (2011); set Lamda_K = 0 to disable this mode of deterioration]

set c_S	    1.0;     #rate of strength deterioration. The default value is 1.0.
set c_C	    1000.0;  #rate of post-capping strength deterioration. The default value is 1.0.
set c_A	    1.0;     #rate of accelerated reloading deterioration. The default value is 1.0.
set c_K	    1.0;     #rate of unloading stiffness deterioration. The default value is 1.0.

set Res_Pos	      [expr $muR/$mu];  #residual strength ratio for positive loading direction
set Res_Neg	      $Res_Pos;         #residual strength ratio for negative loading direction (must be defined as a positive value)

set D_Plus	      1.0;   #rate of cyclic deterioration in the positive loading direction (this parameter is used to create assymetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.
set D_Neg	      1.0;   #rate of cyclic deterioration in the negative loading direction (this parameter is used to create assymetric hysteretic behavior for the case of a composite beam). For symmetric hysteretic response use 1.0.

uniaxialMaterial ModIMKPeakOriented $shearIP $K0 $as_Plus $as_Neg $Vy -$Vy $Lamda_S $Lamda_C $Lamda_A $Lamda_K $c_S $c_C $c_A $c_K $theta_p_Plus $theta_p_Neg $theta_pc_Plus $theta_pc_Neg $Res_Pos $Res_Neg $theta_u_Plus $theta_u_Neg $D_Plus $D_Neg
uniaxialMaterial Elastic   $shearOOP  [expr 5./6.*$G*$L*$t/$h]

element Macroelement3d  1    1  3   2   0.000  0.000  1.000  -1.000  0.000  0.000 -fiberSectionShearModel1d $sec $sec $sec $shearIP $shearOOP $h $driftF $driftS -intWeights 0.1667  0.6666 0.1667


#CONSTRAINTS --------------------------------------------------------------------------------------- 
fix       1    1   1   1   1   1   1 
fix       3    0   0   0   1   1   1 


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
	load      3   0.0  0.0  -207.0e3  0.0   0.0   0.0 
}

# Define analysis parameters
system BandGeneral; 
numberer Plain;
constraints Transformation;	
integrator LoadControl 1;
test NormUnbalance 1.0e-1  50 0;
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
test NormUnbalance 1.0e-1  50  0
algorithm Newton;
constraints Transformation;
numberer Plain;
integrator LoadControl 1;
system BandGeneral;
analysis Static;

# Loading history (mm) --------------------------------------------------------------------------- 	
set dispHistory [list  1 -1 2 -2 3 -3 4 -4 5 -5 0]

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