from math import cos,sin,sqrt,pi
import opensees as ops
## -----------------------------------------------------------------------------------------------------------------------
# Modeled by Chrystal Chern, UC Berkeley.  Date: February 16, 2022.
# ------------------------------------------------------------------------------------------------------------------------
# MAIN MODEL CONFIGURATION SCRIPT (to run, type "source hwd12.1.tcl" into opensees)

# This model is in units of kips, inches, and seconds.
# There are 3 dimensions and 6 degrees of freedom.
# X/Y is the horizontal plane; Z is the vertical axis.
#
# INITIALIZE MODEL -------------------------------------------------------------------------------
#
wipe;                                                        # Clear, 'memory', of, 'all', past, 'model', 'definitions'
model, 'BasicBuilder', -ndm, 3 -ndf, 6;       # Define, 'the', model, builder, ndm=#dimension, ndf=#dofs
dataDir = "datahwd12.1";                      # Set data directory
file, 'mkdir', dataDir;                              # Create, 'data', output, 'folder'
#
# UNITS AND CONSTANTS -------------------------------------------------------------------------------
#
source, 'hwd.matparams.tcl'

#
# NODES -------------------------------------------------------------------------------
#
source, 'hwd12.1.nodes.tcl'

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
#
# COLUMN **PIN** FIBER CROSS-SECTION ASSIGNMENT -------------------------------------------------------------------------------
#
source, 'ColSectionLib.tcl'
# source ColSectionOctIEPin.tcl
Dcol =, 84.0*in
BuildOctColPINSection, 290000, Dcol; # Fiber, cross-section, 'for', columns, 'at', Bents, 2-11 (each, 'have', diameter, 84.0, inches)
Dcol =, 66.0*in
BuildOctColPINSection, 1290000, Dcol; # Fiber, cross-section, 'for', columns, 'at', Bent, 12 (each, 'have', diameter, 66.0, inches) 
Dcol =, 48.0*in
BuildOctColPINSection, 1390000, Dcol; # Fiber, cross-section, 'for', columns, 'at', Bents, 13-14 (each, 'have', diameter, 48.0, inches)

#
# ELEMENTS -------------------------------------------------------------------------------
#
source, 'Elements12.tcl'

# ABUTMENT ELEMENTS
model.element('elasticBeamColumn', 1000, 1020, 1010, (1.0e)+9) Ec, Gc, '[expr', 1.0e+9], '[expr', 1.0e+9], '[expr', 1.0e+9], 301;       # abut-1
model.element('elasticBeamColumn', 1001, 1010, 1030, (1.0e)+9) Ec, Gc, '[expr', 1.0e+9], '[expr', 1.0e+9], '[expr', 1.0e+9], 301;       # abut-1
model.element('elasticBeamColumn', 15000, 15020, 15010, (1.0e)+9) Ec, Gc, '[expr', 1.0e+9], '[expr', 1.0e+9], '[expr', 1.0e+9], 301;       # abut-15NE
model.element('elasticBeamColumn', 15001, 15010, 15030, (1.0e)+9) Ec, Gc, '[expr', 1.0e+9], '[expr', 1.0e+9], '[expr', 1.0e+9], 301;       # abut-15NE
model.element('elasticBeamColumn', 15002, 15050, 15040, (1.0e)+9) Ec, Gc, '[expr', 1.0e+9], '[expr', 1.0e+9], '[expr', 1.0e+9], 301;       # abut-15NR
model.element('elasticBeamColumn', 15003, 15040, 15060, (1.0e)+9) Ec, Gc, '[expr', 1.0e+9], '[expr', 1.0e+9], '[expr', 1.0e+9], 301;       # abut-15NR

#
# IN-SPAN HINGE SPRINGS -------------------------------------------------------------------------------
#
rFlag = 1;  # Rayleigh, damping, ON
source, 'SpanHinges10.tcl'

#
# ABUTMENT SPRINGS -------------------------------------------------------------------------------
#
source, 'Abutments10.tcl'

if {![info 'exists', brace2_main]} {
#
# GRAVITY (STATIC & MODAL) ANALYSIS -------------------------------------------------------------------------------
#
source, 'GravityAnalysis10.tcl'

#
# DYNAMIC ANALYSIS -------------------------------------------------------------------------------
#
ana.recorder('Element', xml=dataDir, )/section4020.xml, ele=4020, section=[1, ], deformations
source, 'time_history.tcl'

