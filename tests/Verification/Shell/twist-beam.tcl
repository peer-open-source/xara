wipe 
model basic -ndm 3 -ndf 6

#nDMaterial ElasticIsotropic $matTag $E $v

nDMaterial ElasticIsotropic 1 29000000 0.22

nDMaterial PlateFiber 601 1
section PlateFiber 1 601 0.32

node 1 0 0.55 1
node 2 0 0 1
node 3 0 -0.55 1
node 4 -0.0717894 0.545295 2
node 5 0 0 2
node 6 0.0717894 -0.545295 2
node 7 -0.14235 0.531259 3
node 8 0 0 3
node 9 0.14235 -0.531259 3
node 10 -0.210476 0.508134 4
node 11 0 0 4
node 12 0.210476 -0.508134 4
node 13 -0.275 0.476314 5
node 14 0 0 5
node 15 0.275 -0.476314 5
node 16 -0.334819 0.436344 6
node 17 0 0 6
node 18 0.334819 -0.436344 6
node 19 -0.388909 0.388909 7
node 20 0 0 7
node 21 0.388909 -0.388909 7
node 22 -0.436344 0.334819 8
node 23 0 0 8
node 24 0.436344 -0.334819 8
node 25 -0.476314 0.275 9
node 26 0 0 9
node 27 0.476314 -0.275 9
node 28 -0.508134 0.210476 10
node 29 0 0 10
node 30 0.508134 -0.210476 10
node 31 -0.531259 0.14235 11
node 32 0 0 11
node 33 0.531259 -0.14235 11
node 34 -0.545295 0.0717894 12
node 35 0 0 12
node 36 0.545295 -0.0717894 12
node 37 -0.55 0 13
node 38 0 0 13
node 39 0.55 0 13
element ShellDKGT 1 1 4 2 1
element ShellDKGT 2 2 4 5 1
element ShellDKGT 3 2 5 3 1
element ShellDKGT 4 3 5 6 1
element ShellDKGT 5 4 7 5 1
element ShellDKGT 6 5 7 8 1
element ShellDKGT 7 5 8 6 1
element ShellDKGT 8 6 8 9 1
element ShellDKGT 9 7 10 8 1
element ShellDKGT 10 8 10 11 1
element ShellDKGT 11 8 11 9 1
element ShellDKGT 12 9 11 12 1
element ShellDKGT 13 10 13 11 1
element ShellDKGT 14 11 13 14 1
element ShellDKGT 15 11 14 12 1
element ShellDKGT 16 12 14 15 1
element ShellDKGT 17 13 16 14 1
element ShellDKGT 18 14 16 17 1
element ShellDKGT 19 14 17 15 1
element ShellDKGT 20 15 17 18 1
element ShellDKGT 21 16 19 17 1
element ShellDKGT 22 17 19 20 1
element ShellDKGT 23 17 20 18 1
element ShellDKGT 24 18 20 21 1
element ShellDKGT 25 19 22 20 1
element ShellDKGT 26 20 22 23 1
element ShellDKGT 27 20 23 21 1
element ShellDKGT 28 21 23 24 1
element ShellDKGT 29 22 25 23 1
element ShellDKGT 30 23 25 26 1
element ShellDKGT 31 23 26 24 1
element ShellDKGT 32 24 26 27 1
element ShellDKGT 33 25 28 26 1
element ShellDKGT 34 26 28 29 1
element ShellDKGT 35 26 29 27 1
element ShellDKGT 36 27 29 30 1
element ShellDKGT 37 28 31 29 1
element ShellDKGT 38 29 31 32 1
element ShellDKGT 39 29 32 30 1
element ShellDKGT 40 30 32 33 1
element ShellDKGT 41 31 34 32 1
element ShellDKGT 42 32 34 35 1
element ShellDKGT 43 32 35 33 1
element ShellDKGT 44 33 35 36 1
element ShellDKGT 45 34 37 35 1
element ShellDKGT 46 35 37 38 1
element ShellDKGT 47 35 38 36 1
element ShellDKGT 48 36 38 39 1

fix 1 1 1 1 1 1 1 
fix 2 1 1 1 1 1 1
fix 3 1 1 1 1 1 1

recorder Node -file disp.txt -time -nodeRange 37 39 -dof 1  disp
recorder Node -file disp2.txt -time -nodeRange 37 39 -dof 2  disp;

pattern Plain 1 Linear {
load 38 1 0 0 0 0 0 
load 38 0 0 0 0 0 0 

}

pragma analysis off
export twisted-beam.vtk

constraints Plain
test NormDispIncr 1.0e-6 200 
algorithm Newton -count 100
numberer RCM
system SparseGEN
integrator LoadControl 0.001			
analysis Static				
analyze 1000;	


