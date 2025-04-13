from math import cos,sin,sqrt,pi
import opensees as ops
## -----------------------------------------------------------------------------------------------------------------------
# Modeled by Chrystal Chern, UC Berkeley.  Date: February 2, 2022.
# ------------------------------------------------------------------------------------------------------------------------
# MAIN MODEL CONFIGURATION SCRIPT (to run, type "source hwd1.1.tcl" into opensees)

# This model is in units of kips, inches, and seconds.
# There are 3 dimensions and 6 degrees of freedom.
# X/Y is the horizontal plane; Z is the vertical axis.
#
# INITIALIZE MODEL ----------------------------------------
#
#wipe;                                                        # Clear memory of all past model definitions
model, 'BasicBuilder', -ndm, 3 -ndf, 6;       # Define, 'the', model, builder, ndm=#dimension, ndf=#dofs
dataDir = "datahwd1.1";                      # Set data directory
file, 'mkdir', dataDir;                              # Create, 'data', output, 'folder'
#
# UNITS AND CONSTANTS -------------------------------------
#

source, 'hwd.matparams.tcl'

#
# NODES ---------------------------------------------------
#
# See Haywardbridge_coordinates_NodeAsgn.xlsx for node coordinates and node assignment.
# COLUMN BENT NODES (Node numbering scheme: 100*Bent#+Node#)
source, 'Nodes1.tcl'

#
# CONSTRAINTS ---------------------------------------------
#
source, 'Constraints1.tcl'

#
# COLUMN FIBER CROSS-SECTION ASSIGNMENT -------------------
#
source, 'ColSectionCirc.tcl'

#
# ELEMENTS ------------------------------------------------
#
source, 'Elements1.tcl'

if {![info 'exists', brace2_main]} {
print("RUNNING CONF AS SCRIPT")
#
# GRAVITY (STATIC & MODAL) ANALYSIS -----------------------
#
source, 'GravityAnalysis1.tcl'

# DYNAMIC ANALYSIS ----------------------------------------
source, 'time_history.tcl'

