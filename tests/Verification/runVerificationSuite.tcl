# script to run all the verification scripts
# PASS/FAILURE results in file README.md when run

# open results file such that it is cleared out of any data
set results [open README.md w]
puts $results "| Status | Notes |\n|--------|------------------------------|"
close $results

cd Basic
source sdofTransient.tcl
source SmallEigen.tcl
source NewmarkIntegrator.tcl
source mdofModal.tcl
cd ..

source Truss/PlanarTruss.tcl
source Truss/PlanarTruss.Extra.tcl

source Frame/PortalFrame2d.tcl
source Frame/EigenFrame.tcl
source Frame/EigenFrame.Extra.tcl
source Frame/AISC25.tcl

source Plane/PlaneStrain.tcl
source Plane/QuadBending.tcl

# Shells
source Shell/PinchedCylinder.tcl
source Shell/PlanarShearWall.tcl

