wipe 
model basic -ndm 3 -ndf 6

#nDMaterial ElasticIsotropic $matTag $E $v

nDMaterial ElasticIsotropic 1 432000000 0

nDMaterial PlateFiber 601 1
section PlateFiber 1 601 0.25

node 1 0 0 5.84889

node 2 8.5505 0 4.3412
node 3 16.0697 0 0

node 4 0 12.5 5.84889

node 5 8.5505 12.5 4.3412

node 6 16.0697 12.5 0

node 7 0 25 5.84889

node 8 8.5505 25 4.3412

node 9 16.0697 25 0


element ShellDKGT 1 1 2 4 1
element ShellDKGT 2 4 2 5 1
element ShellDKGT 3 5 2 3 1
element ShellDKGT 4 5 3 6 1
element ShellDKGT 5 7 4 5 1
element ShellDKGT 6 7 5 8 1
element ShellDKGT 7 8 5 6 1
element ShellDKGT 8 8 6 9 1


fixY  0  1 0 1 0 1 0
fixX  0  1 0 0 0 1 1
fixY 25  0 1 0 1 0 1

recorder Node -file disp.txt -time -node 9 -dof 3  disp

pattern Plain 1 Linear {

load 1 0 0 -2454.369261 0 0 0
load 2 0 0 -4908.73852  0 0 0
load 3 0 0 -2454.369261 0 0 0
load 4 0 0 -4908.73852  0 0 0
load 5 0 0 -9817.477044 0 0 0
load 6 0 0 -4908.73852  0 0 0
load 7 0 0 -2454.369261 0 0 0
load 8 0 0 -4908.73852  0 0 0
load 9 0 0 -2454.369261 0 0 0




}


constraints Plain
test NormDispIncr 1.0e-6 2000  2
algorithm KrylovNewton
numberer RCM
system SparseGEN
integrator LoadControl 0.001   
analysis Static    
analyze 1000; 



