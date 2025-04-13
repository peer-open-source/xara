from math import cos,sin,sqrt,pi
import opensees as ops
## -----------------------------------------------------------------------------------------------------------------------
# BRACE2, UC Berkeley.  Date: March 22, 2022.
# ------------------------------------------------------------------------------------------------------------------------
# Batch Simulation (to run, type "source simulation.tcl into opensees)
# Put this file with (1) the GMs (GM folder (e.g., Hayward) and the GM summary (e.g., Hayward_GM.txt) and (2) the model files (e.g., hwd12_Recorders.tcl)
# Change the fp GM summary txt file if needed
# Change the model file (e.g., source hwd12_Recorders.tcl) if needed
# Change dt (in the *.dat file, all TH are unified into same dt, e.g., 0.01 sec for Hayward) if needed
# Remove inFilelong, inFiletrans, inFilevert in the the model file
# Remove dt, NumPts in the model file


wipe;                                                        # Clear, 'memory', of, 'all', past, 'model', 'definitions'

# Place this file with 
# MAY NEED TO CHANGE
fp = [open "Hayward_3GM.txt" r]
content = [read, fp]
close, fp

# MAY NEED TO CHANGE
dt =, 0.01*sec

# Next two commands reads all rows
lst_lines = [split content "\n"]
count = 0
#puts lst_lines

foreach, 'str_line', lst_lines {
       incr 'count', 1
       GMdataDir = "GM{count}"
       puts GMdataDir
    lst_clmns = [split str_line]
    inFilelong_s = [lindex '[split', lst_clmns] 0]
    inFilelong = "{inFilelong_s}.dat"
    inFiletrans_s = [lindex '[split', lst_clmns] 1]
    inFiletrans = "{inFiletrans_s}.dat"
    inFilevert_s =  [lindex '[split', lst_clmns] 2]
    inFilevert = "{inFilevert_s}.dat"
    SF = [lindex '[split', lst_clmns] 3]
    NumPts =      [lindex '[split', lst_clmns] 4]
       # NumPts = 1500;  # Test small number of points for debugging

# MAY NEED TO CHANGE
source, 'hwd10.1_batch.tcl'

