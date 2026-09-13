# Units: m, s, N
wipe

# Paramters for the analysis (change)
set   maxK            [expr 0.002];       	# maximum sectional deformation
set   numIncr           250;             	# number of increments
set   controlled_dof    1;    				# controlled degree of freedom
set   alpha 0.5;
set   n_cycles          1.0;

# Set parameters for the loading conditions
set NormalLoad    -0.5e6;
set ShearLoad     0.0e3;

# Create ModelBuilder (with two-dimensions (-ndm) and 2 DOF/node (-ndf))
model basic -ndm 3 -ndf 6

# -----------------------------------------------------------------------------------
# Definition of the geometry
# Create nodes
#     tag       X       Y       Z 
node  1001      0.0     0.0     0.0
node  1002      0.0     0.0     0.0

# -----------------------------------------------------------------------------------
# Definition of the materials
# MASONRY      
set    E    [expr 3000e6]
set    G    [expr 1000e6]
set c		0.03e6;
set GfII		0.2e1;
set mu    	0.4;
set h       0.01;

#Declare material
#nDMaterial ElasticIsotropic 10 [expr $E] 0.25
nDMaterial FrictionCohesion3d  10   $G  $mu  $c  $GfII -imposedSigma -0.1e6
puts "Material defined."

# -----------------------------------------------------------------------------------
# Definition of the interface
#            

#                      	$eleTag	$iNode	$jNode 	$matTag <$uniTag> 	<-orient $x1 $x2 $x3 $yp1 $yp2 $yp3>
#element ZeroLengthND2DOF 	2001 	1001 	1002 	10	 1.0 0.01  -orient 0.0  1.0  0.0   -1.0  0.0  0.0
element zeroLengthND     1    1001  1002    10   -orient 0. 0. 1.   1.  0.  0.  			                  
puts "Interface defined."

# Constraints
#     tag       DX   DY   DZ   RX   RY   RZ
fix   1001      1    1    1    1    1    1
fix   1002      0    0    1    1    1    1  

# -----------------------------------------------------------------------------------
# Recorders
recorder    Node      -file node_displ.out         -time     -node 1002     -dof  1 2 3  disp;
recorder    Node      -file node_force.out         -time     -node 1001     -dof  1 2 3  reaction;
#recorder 	Element	  -file interface_stress.out        -precision  4  -ele  2001   material     stress;
#recorder 	Element	  -file interface_strain.out        -precision  4  -ele  2001   material     strain;
#recorder 	Element	  -file interface_stiffness.out     -precision  4  -ele  2001   material     stiffness;
#recorder 	Element	  -file interface_mode.out          -precision  4  -ele  2001   material     mode;

# -----------------------------------------------------------------------------------
# CONSTANT LOADS
# Define constant loads
#pattern  Plain   3001   "Constant" {
#	load  1002   0.0   $NormalLoad
#}
# Define analysis parameters
integrator LoadControl 0.0 1 0 0
system SparseGeneral -piv;
test  NormUnbalance  1.0e-4 100   0
numberer Plain
constraints Plain
algorithm Newton
analysis Static
#analyze 1
#loadConst -time 0.0
#puts "Vertical load applied"


# -----------------------------------------------------------------------------------
# INCREMENTAL ANALYSIS 
# Define reference force
pattern Plain 3002 "Linear" {
		load  1002      0.  1.  0.    0.  0.  0.
}

#history in microstrain
#set dispHistory [list -0.2 0.2 -0.4 0.4 -0.6 0.6]
#set dispHistory [list 0.50 -0.50]
set dispHistory [list 0.10 -0.10 0.20 -0.20 0.30 -0.3 0.0]

set d_n 0.0
set node2 0.0
set dispOld 0.0
# increment in microstrain
set incr 0.01
set incr [expr $incr/1000.0]  
set loadingDir    2
foreach dMax $dispHistory {
	set dMax [expr $dMax/1000.0]
	
	if {$dMax > $dispOld} {
	# positive loading
	while {$node2<$dMax} {
		set d_n [expr $d_n+$incr]
		#                                    $node $dof       $incr  
		integrator    DisplacementControl     1002     $loadingDir    [expr 1.0*$incr] 
        set ok [analyze 1]
		if {$ok != 0} {
			puts "Failed at displacement " $node2;
			break
		}
		set node2 [nodeDisp 1002 $loadingDir]
	}
	
	} else {
	
	#negative loading
	set dMin [expr $dMax]
	# first negative loading
	while {$node2 > $dMin} {
		set d_n [expr $d_n-$incr]
		#                                    $node $dof       $incr  
		integrator    DisplacementControl     1002      $loadingDir     [expr -1.0*$incr] 
        set ok [analyze 1]
	if {$ok != 0} {
		puts "Failed at displacement " $node2;
		break
	}
	set node2 [nodeDisp 1002 $loadingDir]
	}
	
	}
	set dispOld $node2
}
puts "Analysis completed."
wipe
#exit