source truss.tcl

# build the components for the analysis object
system BandSPD
constraints Plain
integrator LoadControl 1.0
algorithm Linear
numberer RCM

# create the analysis object 
analysis Static 


# perform the analysis
analyze 1

# print the results at node and at all elements
print node 4
print ele


if {[getNumElements] != 3} {
  puts "FAILED - getNumElements"
} else {
  puts "PASSED - getNumElements"
}

getNodeTags

print

remove sp 2
remove sp 3 2 ; # sp with tag 5

print
