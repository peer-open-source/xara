wipe 
model basic -ndm 3 -ndf 6

#nDMaterial ElasticIsotropic $matTag $E $v

nDMaterial ElasticIsotropic 1 432000000 0

nDMaterial PlateFiber 601 1
section PlateFiber 1 601 0.25

node 1 0 0 5.84889
node 2 4.3412 0 5.46908
node 3 8.5505 0 4.3412
node 4 12.5 0 2.49952
node 5 16.0697 0 0
node 6 0 6.25 5.84889
node 7 4.3412 6.25 5.46908
node 8 8.5505 6.25 4.3412
node 9 12.5 6.25 2.49952
node 10 16.0697 6.25 0
node 11 0 12.5 5.84889
node 12 4.3412 12.5 5.46908
node 13 8.5505 12.5 4.3412
node 14 12.5 12.5 2.49952
node 15 16.0697 12.5 0
node 16 0 18.75 5.84889
node 17 4.3412 18.75 5.46908
node 18 8.5505 18.75 4.3412
node 19 12.5 18.75 2.49952
node 20 16.0697 18.75 0
node 21 0 25 5.84889
node 22 4.3412 25 5.46908
node 23 8.5505 25 4.3412
node 24 12.5 25 2.49952
node 25 16.0697 25 0


element ShellDKGT 1 1 2 6 1
element ShellDKGT 2 2 7 6 1
element ShellDKGT 3 2 3 7 1
element ShellDKGT 4 3 8 7 1
element ShellDKGT 5 3 4 8 1
element ShellDKGT 6 4 9 8 1
element ShellDKGT 7 4 5 9 1
element ShellDKGT 8 5 10 9 1
element ShellDKGT 9 6 7 11 1
element ShellDKGT 10 7 12 11 1
element ShellDKGT 11 7 8 12 1
element ShellDKGT 12 8 13 12 1
element ShellDKGT 13 8 9 13 1
element ShellDKGT 14 9 14 13 1
element ShellDKGT 15 9 10 14 1
element ShellDKGT 16 10 15 14 1
element ShellDKGT 17 11 12 16 1
element ShellDKGT 18 12 17 16 1
element ShellDKGT 19 12 13 17 1
element ShellDKGT 20 13 18 17 1
element ShellDKGT 21 13 14 18 1
element ShellDKGT 22 14 19 18 1
element ShellDKGT 23 14 15 19 1
element ShellDKGT 24 15 20 19 1
element ShellDKGT 25 16 17 21 1
element ShellDKGT 26 17 22 21 1
element ShellDKGT 27 17 18 22 1
element ShellDKGT 28 18 23 22 1
element ShellDKGT 29 18 19 23 1
element ShellDKGT 30 19 24 23 1
element ShellDKGT 31 19 20 24 1
element ShellDKGT 32 20 25 24 1

fixY  0  1 0 1 0 1 0
fixX  0  1 0 0 0 1 1
fixY 25  0 1 0 1 0 1

recorder Node -file disp.txt -time -node 25 -dof 3  disp

pattern Plain 1 Linear {

load 1 0 0 -613.5923152 0 0 0
load 2 0 0 -1227.18463 0 0 0
load 3 0 0 -1227.18463 0 0 0
load 4 0 0 -1227.18463 0 0 0
load 5 0 0 -613.5923152 0 0 0
load 6 0 0 -1227.18463 0 0 0
load 7 0 0 -2454.369261 0 0 0
load 8 0 0 -2454.369261 0 0 0
load 9 0 0 -2454.369261 0 0 0
load 10 0 0 -1227.18463 0 0 0
load 11 0 0 -1227.18463 0 0 0
load 12 0 0 -2454.369261 0 0 0
load 13 0 0 -2454.369261 0 0 0
load 14 0 0 -2454.369261 0 0 0
load 15 0 0 -1227.18463 0 0 0
load 16 0 0 -1227.18463 0 0 0
load 17 0 0 -2454.369261 0 0 0
load 18 0 0 -2454.369261 0 0 0
load 19 0 0 -2454.369261 0 0 0
load 20 0 0 -1227.18463 0 0 0
load 21 0 0 -613.5923152 0 0 0
load 22 0 0 -1227.18463 0 0 0
load 23 0 0 -1227.18463 0 0 0
load 24 0 0 -1227.18463 0 0 0
load 25 0 0 -613.5923152 0 0 0


}


constraints Plain
test NormDispIncr 1.0e-6 2000  2
algorithm KrylovNewton
numberer RCM
system SparseGEN
integrator LoadControl 0.001   
analysis Static    
analyze 1000; 



