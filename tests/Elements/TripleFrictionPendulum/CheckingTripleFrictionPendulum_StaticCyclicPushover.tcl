# 	MODEL FOR CHECKING TripleFrictionPendulum ELEMENT
# 	TYPE: STATIC CYCLIC PUSHOVER - CONSTANT FRICTION COEFFICIENTS
# 	CREATED BY: NHAN DAO, UNR
################################################
# 	INPUT INFORMATION
# 	INPUT INFORMATION
# 	INPUT INFORMATION

set Drtn 		0.; 		# Angle of direction to check, degree, -90 <= Drtn <=90
set OutDir 		Output;		# Output folder
set OutFile 	Disp.txt; 	# Output file
# bearing information
set L1 	0.36; # effective length
set L2 	1.25;
set L3 1.25;
set mu1 0.012; # friction coefficient
set mu2 0.052;
set mu3 0.14;
set d1 0.1; # pendulum displacement limit
set d2 0.2;
set d3 0.2;
set uy 0.001; # displacement where sliding starts
set kvc 1000.; # vertical compression stiffness
set kvt 1.; # vertical tension stiffness
set minFv 0.1; # minimum compression force in the bearing

set W 1000.; # static weight supported by the bearing
set tol 1.e-5; # relative tolerance for checking convergence

# 	END OF INPUT INFORMATION
# 	END OF INPUT INFORMATION
# 	END OF INPUT INFORMATION
################################################
# CREATE MODEL AND ANAYLIZE
wipe;
model basic -ndm 3 -ndf 6;
#----------------------------------------------
# Creating nodes
node 1 0. 0. 0.;
node 2 0. 0. 0.;
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
# Creating materials for compression and rotation behaviors
uniaxialMaterial 	Elastic 1 $kvc;
uniaxialMaterial 	Elastic 2 	10.;
#-----------------------------------------------
# Creating TripleFrictionPendulum element
# element TripleFrictionPendulum $eleTag $iNode $jNode $frnTag1 $frnTag2 $frnTag3 $vertMatTag $rotZMatTag $rotXMatTag $rotYMatTag $L1 $L2 $L3 $d1 $d2 $d3 $W $uy $kvt $minFv $tol
element TripleFrictionPendulum 1 1 2  1 2 3 1 2 2 2  $L1 $L2 $L3 $d1 $d2 $d3 $W $uy $kvt $minFv $tol;
#---------------------------------------------------
# Applying static vertical load
pattern Plain 1 Linear {
    load 2 0. 0. [expr -$W] 0. 0. 0.
}
#----------------------------------------------
# Creating analysis object
numberer RCM
system BandGeneral
constraints Transformation
test NormDispIncr 1e-8 100;
algorithm Newton;
integrator LoadControl 1.;
analysis Static
analyze 1
loadConst -time 0.0;
#--------------------------------------------
# Applying unit load for pushover
pattern Plain 2 Linear {
    load 2 [tcl::mathfunc::cos [expr $Drtn*3.141592653598793/180]] [tcl::mathfunc::sin [expr $Drtn*3.141592653598793/180]] [expr 0.] 0. 0. 0.
}
#----------------------------------------------
# Creating analysis object
integrator DisplacementControl 2 1 .1;
analysis Static
#-----------------------------------------------------
# create recorder object
file mkdir $OutDir;
recorder Node -file $OutDir/$OutFile -time -nodes 2 -dof 1 2 disp
#----------------------------------------------------
# Cyclic pushover
foreach d {.2 -.4 .4	-.4	.4	-.4	.4	-.4	.4	-.4	.4	-.4	.6	-.8	.8	-.8	 1.	-1.2	1.2}   {
#foreach d {.49}  
	incr cnt 1;
	puts "step: $cnt , d= $d";
	if {$d >= 0} {
		set du 	.001;		
	} else {	
		set du 	-.001;		
	}
	set nstep [expr int($d/$du)];
	if {$Drtn >= -45 && $Drtn <= 45} {
		set du [expr $du*[tcl::mathfunc::cos [expr $Drtn*3.141592653598793/180]]];
		integrator DisplacementControl 2 1 $du;
	} else {
		set du [expr $du*[tcl::mathfunc::sin [expr $Drtn*3.141592653598793/180]]];
		integrator DisplacementControl 2 2 $du;
	}
	set ok 0;
	set temp 1;
	while {$temp <= $nstep && $ok == 0} {
		incr temp 1;
		set ok [analyze 1]
		if {$ok != 0} {
			puts "Trying NewtonWithLineSearch .."
			algorithm NewtonLineSearch 
			set ok [analyze 1]
			algorithm Newton
		}
		if {$ok != 0} {
			puts "Trying KrylovNewton .."
			test EnergyIncr 1.0e-10 100 0
			algorithm KrylovNewton
			set ok [analyze 1]
			test EnergyIncr 1.0e-10 100 0
			algorithm Newton
		}
		if {$ok != 0} {
			puts "Trying Broyden .."
			algorithm Broyden 100
			set ok [analyze 1]
			algorithm Newton
		}
	}; # end while loop	
}
wipe
exit
