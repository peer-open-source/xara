# 	Example 4:
# 	MODEL FOR CHECKING TripleFrictionPendulum ELEMENT
# 	TYPE: BI-DIRECTIONAL DYNAMIC ANALYSIS, CONSTANT FRICTION COEFFICIENTS
# 	CREATED BY: NHAN DAO, UNR
################################################
# 	INPUT INFORMATION
# 	INPUT INFORMATION
# 	INPUT INFORMATION

set GMFilex		TCU065-E.ATH;	# Input ground motion
set GMFiley		TCU065-N.ATH;	# Input ground motion
set nSteps 		14000; 			# Number of steps to analyze
set GMFac		9.81;			# GM factor = g = 9.81 m/s^2
set dt 			0.005; 			# Time step of ground motion (s)
set OutDir 		Output;			# Output folder
set OutFile1 	Disp.txt; 		# Output file
set OutFile2	Reaction.txt;
# bearing information
set L1 	0.36; # effective length (m)
set L2 	1.25;
set L3 1.25;
set mu1 0.012; # friction coefficient
set mu2 0.052;
set mu3 0.14;
set d1 0.1; # pendulum displacement limit (m)
set d2 0.2;
set d3 0.2;
set uy 0.001; # displacement where sliding starts (m)
set kvc 1000000.; # vertical compression stiffness (N/m)
set kvt 1.; # vertical tension stiffness (N/m)
set minFv 0.1; # minimum compression force in the bearing (N)

set W 1000.; # static weight supported by the bearing (N)
set tol 1.e-5; # relative tolerance for checking convergence
set g 	9.81; # gravity acceleration (m/s^2)

# 	END OF INPUT INFORMATION
# 	END OF INPUT INFORMATION
# 	END OF INPUT INFORMATION
################################################
# CREATE MODEL AND ANAYLIZE
wipe;
model basic -ndm 3 -ndf 6;
#----------------------------------------------
# Creating nodes
node 1 0 0 0; # End i
node 2 0 0 0 -mass [expr $W/$g] [expr $W/$g] [expr $W/$g] 1. 1. 1.; # End j
#---------------------------------------------
# Applying kinematic boundary condition
fix 1 	1 1 1 1 1 1;
#-----------------------------------------------
# Creating friction models
# frictionModel Coulomb tag mu
frictionModel Coulomb 1 $mu1
frictionModel Coulomb 2 $mu2
frictionModel Coulomb 3 $mu3
#-----------------------------------------------
# Creating material for compression and rotation behaviors
uniaxialMaterial 	Elastic 1 $kvc;
uniaxialMaterial 	Elastic 2 	10.;
#-----------------------------------------------
# Creating TripleFrictionPendulum element
# element TripleFrictionPendulum $eleTag $iNode $jNode $frnTag1 $frnTag2 $frnTag3 $vertMatTag $rotZMatTag $rotXMatTag $rotYMatTag $L1 $L2 $L3 $d1 $d2 $d3 $W $uy $kvt $minFv $tol
element TripleFrictionPendulum 1 1 2  1 2 3 1 2 2 2  $L1 $L2 $L3 $d1 $d2 $d3 $W $uy $kvt $minFv $tol;
#----------------------------------------------
# Applying static vertical load
pattern Plain 1 Linear {
    load 2 0. 0. [expr -$W] 0. 0. 0.;
}
#-----------------------------------------------
# Creating analysis object
numberer RCM;
system BandGeneral;
constraints Transformation;
test NormDispIncr 1e-5 100;
algorithm Newton;
integrator LoadControl 1.;
analysis Static;
analyze 1;
loadConst -time 0.0;
#-------------------------------------------
# apply damping
set Damp 0.01;
set	EigVal	[eigen	-genBandArpack	1];
set	w1	[expr	pow($EigVal,0.5)];
set a1 [expr 2*$Damp/$w1];
rayleigh	0.0	0.0	0.0	$a1; # converged tangent stiffness proportional damping
#-------------------------------------------------
# Creating analysis object
integrator Newmark 0.5 0.25;
analysis Transient
#-----------------------------------------------------
# create recorder object
file mkdir $OutDir;
recorder Node -file $OutDir/$OutFile1 -time -nodes 2 -dof 1 2 3 disp
recorder Node -file $OutDir/$OutFile2 -time -nodes 1 -dof 1 2 3 reaction
#----------------------------------------------------
set comp1 "Series -dt $dt -filePath $GMFilex -factor $GMFac"
set comp2 "Series -dt $dt -filePath $GMFiley -factor $GMFac"

pattern UniformExcitation 2 1 -accel $comp1
pattern UniformExcitation 3 2 -accel $comp2

#Create the analysis parameters
set DtAnalysis		[expr	$dt]
set	testTypeDynamic			NormDispIncr
set	TolDynamic				1.e-5;
set	maxNumIterDynamic		50;

puts "Dynamics analysis starts..."
analyze $nSteps $DtAnalysis;
wipe
exit
