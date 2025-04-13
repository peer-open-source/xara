from math import cos,sin,sqrt,pi
import opensees as ops
## -----------------------------------------------------------------------------------------------------------------------
# Modeled by Chrystal Chern, UC Berkeley.  Date: February 16, 2022.
# ------------------------------------------------------------------------------------------------------------------------
# MAIN MODEL CONFIGURATION SCRIPT (to run, type "source hwd7.1.tcl" into opensees)

# This model is in units of kips, inches, and seconds.
# There are 3 dimensions and 6 degrees of freedom.
# X/Y is the horizontal plane; Z is the vertical axis.
#
# INITIALIZE MODEL -------------------------------------------------------------------------------
#
#wipe;                                                        # Clear memory of all past model definitions
model, 'BasicBuilder', -ndm, 3 -ndf, 6;       # Define, 'the', model, builder, ndm=#dimension, ndf=#dofs
dataDir = "datahwd7.1";                      # Set data directory
file, 'mkdir', dataDir;                              # Create, 'data', output, 'folder'
#
# UNITS AND CONSTANTS -------------------------------------------------------------------------------
#
source, 'hwd.matparams.tcl'

#
# NODES -------------------------------------------------------------------------------
#
# See Haywardbridge_coordinates_NodeAsgn.xlsx for node coordinates and node assignment.
# COLUMN BENT NODES (Node numbering scheme: 100*Bent#+Node#)
source, 'Nodes7.tcl'

#
# CONSTRAINTS -------------------------------------------------------------------------------
#
source, 'Constraints.tcl'

#
# COLUMN FIBER CROSS-SECTION ASSIGNMENT -------------------------------------------------------------------------------
#
source, 'ColSectionLib.tcl', 
AssignHaywardSections, 'BuildOctColSectionInelMod', 'BuildWideOctColSectionInelMod'

#
# ELEMENTS -------------------------------------------------------------------------------
#
source, 'Elements7.tcl'

# ABUTMENT ELEMENTS
model.element('elasticBeamColumn', 1000, 1020, 1010, (1.0e)+9) Ec, Gc, '[expr', 1.0e+9], '[expr', 1.0e+9], '[expr', 1.0e+9], 301;       # abut-1
model.element('elasticBeamColumn', 1001, 1010, 1030, (1.0e)+9) Ec, Gc, '[expr', 1.0e+9], '[expr', 1.0e+9], '[expr', 1.0e+9], 301;       # abut-1
model.element('elasticBeamColumn', 15000, 15020, 15010, (1.0e)+9) Ec, Gc, '[expr', 1.0e+9], '[expr', 1.0e+9], '[expr', 1.0e+9], 301;       # abut-15NE
model.element('elasticBeamColumn', 15001, 15010, 15030, (1.0e)+9) Ec, Gc, '[expr', 1.0e+9], '[expr', 1.0e+9], '[expr', 1.0e+9], 301;       # abut-15NE
model.element('elasticBeamColumn', 15002, 15050, 15040, (1.0e)+9) Ec, Gc, '[expr', 1.0e+9], '[expr', 1.0e+9], '[expr', 1.0e+9], 301;       # abut-15NR
model.element('elasticBeamColumn', 15003, 15040, 15060, (1.0e)+9) Ec, Gc, '[expr', 1.0e+9], '[expr', 1.0e+9], '[expr', 1.0e+9], 301;       # abut-15NR
#
# COLUMN PINS -------------------------------------------------------------------------------
#
source, 'ColumnPins7.tcl'

#
# IN-SPAN HINGE SPRINGS -------------------------------------------------------------------------------
#
source, 'SpanHinges07.tcl'

#
# ABUTMENT SPRINGS -------------------------------------------------------------------------------
#
source, 'AbutmentSprings3.tcl'

if {![info 'exists', brace2_main]} {
#
# GRAVITY (STATIC & MODAL) ANALYSIS -------------------------------------------------------------------------------
#
source, 'gravity_2.tcl'

# DYNAMIC ANALYSIS -------------------------------------------------------------------------------
source, 'time_history.tcl'
