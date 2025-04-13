
proc reset {} {
  wipeAnalysis 
  numberer RCM
  algorithm Linear
  constraints Plain
  system "FullGeneral"
  # analysis "Transient"
}

 
# Mass
reset
integrator "GimmeMCK" 1.0 0.0 0.0
analysis "Transient"
analyze 1 1.0  
printA "-file" "M.out"
 
# Stiffness
reset
integrator "GimmeMCK" 0.0 0.0 1.0
analysis "Transient"
analyze 1 1.0
printA "-file" "K.out"
 
# Damping
reset
integrator "GimmeMCK" 0.0 1.0 0.0
analysis "Transient"
analyze 1 1.0
printA "-file" "C.out"
