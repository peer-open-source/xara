from math import cos,sin,sqrt,pi
import opensees as ops
#
# COLUMN FIBER CROSS-SECTION ASSIGNMENT -------------------------------------------------------------------------------
#
# (Column Numbering scheme: 1000*Bent#+10*Column#)
# (Column#: 1=Left, 2=Right, 3=Center, 4=Single)
# Column fiber element tag (ColSecTag) follows column numbering scheme.
# Read in the column heights:
fpHcol = [open "./Dimensions/Hcol.txt" r];
HcolList = [read, fpHcol];
close, fpHcol
# Read in the number of longitudinal bars for each column:
fpnLbar = [open "./Dimensions/nLbar.txt" r];
nLbarList = [read, fpnLbar];
close, fpnLbar
# Read in the diameter of longitudinal bars for each column:
fpDLbar = [open "./Dimensions/DLbar.txt" r];
DLbarList = [read, fpDLbar];
close, fpDLbar
# Read in the spacing of transverse spiral reinforcement for each column:
fpsTbar = [open "./Dimensions/sTbar.txt" r];
sTbarList = [read, fpsTbar];
close, fpsTbar
# Procedure for column fiber cross-section assignment (ELASTIC material properties):

source, 'ColSectionLib.tcl'
colcmd = "BuildOctColSectionElastic"
colwidecmd = "BuildWideOctColSectionElastic"
# BUILD COLUMN SECTIONS
# Bents 2-11
for {ib = 2} {ib <= 11} {ib += 1} {
       Dcol =, 84.0*in
       ic = 1;                                                                                    # Left Column
       ColSecTag = 1000*ib+10*ic;                                   # Column's fiber model.element(tag. Follows column numbering scheme.)
       Hcol = [lindex HcolList '[expr', (ib-2)*2]];
       nLbar = [lindex nLbarList '[expr', (ib-2)*2]];
       DLbar = [lindex DLbarList '[expr', (ib-2)*2]];
       sTbar = [lindex sTbarList '[expr', (ib-2)*2]];
       colcmd ColSecTag  Dcol  nLbar  DLbar sTbar; # Fiber cross-section
       ic = 2;                                                                                    # Right Column
       ColSecTag = 1000*ib+10*ic;                                   # Column's fiber model.element(tag. Follows column numbering scheme.)
       Hcol = [lindex HcolList '[expr', (ib-2)*2+1]];
       nLbar = [lindex nLbarList '[expr', (ib-2)*2+1]];
       DLbar = [lindex DLbarList '[expr', (ib-2)*2+1]];
       sTbar = [lindex sTbarList '[expr', (ib-2)*2+1]];
       colcmd ColSecTag  Dcol  nLbar  DLbar sTbar; # Fiber cross-section

# Bent 12
ib = 12;
Dcol =, 66.0*in
ic = 1;                                                                                    # Left, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 20];
nLbar = [lindex, nLbarList, 20];
DLbar = [lindex, DLbarList, 20];
sTbar = [lindex, sTbarList, 20];
colcmd, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
ic = 2;                                                                                    # Right, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 22];
nLbar = [lindex, nLbarList, 22];
DLbar = [lindex, DLbarList, 22];
sTbar = [lindex, sTbarList, 22];
colcmd, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
ic = 3;                                                                                    # Center, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 21];
nLbar = [lindex, nLbarList, 21];
DLbar = [lindex, DLbarList, 21];
sTbar = [lindex, sTbarList, 21];
colcmd, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
# Bent 13 NE
ib = 13;
Dcol =, 48.0*in
ic = 1;                                                                                    # Left, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 23];
nLbar = [lindex, nLbarList, 23];
DLbar = [lindex, DLbarList, 23];
sTbar = [lindex, sTbarList, 23];
colcmd, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
ic = 2;                                                                                    # Right, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 24];
nLbar = [lindex, nLbarList, 24];
DLbar = [lindex, DLbarList, 24];
sTbar = [lindex, sTbarList, 24];
colcmd, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
# Bent 13 NR (WIDE (INTERLOCKING) SECTION)
ib = 13;
Dcol = 48.0*in;                                                         # Shorter, Width, of, octagonal, column (to, flat, sides)
ic = 4;                                                                                    # Single, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 25];                                           # Column, Height
nLbar = [lindex, nLbarList, 25];                                           # Number, of, main (outer) longitudinal, bars
DLbar = [lindex, DLbarList, 25];                                           # Diameter, of, main (outer) longitudinal, bars
sTbar = [lindex, sTbarList, 25];                                           # Spacing, of, transverse, spiral, reinforcement
Wcol = 72.0*in;                                                         # Longer, Width, of, octagonal, column (to, flat, sides)
nLbar2 = 8;                                                                              # Number, of, secondary (inner) longitudinal, bars
DLbar2 = 0.625*in;                                                  # Diameter, of, secondary (inner) longitudinal, bars (#5, rebar)
colwidecmd, ColSecTag, Dcol, Wcol, nLbar, nLbar2, DLbar, DLbar2, sTbar; # Fiber, cross-section
# Bent 14 NE
ib = 14;
Dcol =, 48.0*in
ic = 1;                                                                                    # Left, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 26];
nLbar = [lindex, nLbarList, 26];
DLbar = [lindex, DLbarList, 26];
sTbar = [lindex, sTbarList, 26];
colcmd, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
ic = 2;                                                                                    # Right, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 28];
nLbar = [lindex, nLbarList, 28];
DLbar = [lindex, DLbarList, 28];
sTbar = [lindex, sTbarList, 28];
colcmd, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
ic = 3;                                                                                    # Center, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 27];
nLbar = [lindex, nLbarList, 27];
DLbar = [lindex, DLbarList, 27];
sTbar = [lindex, sTbarList, 27];
colcmd, ColSecTag, Dcol, nLbar, DLbar, sTbar; # Fiber, cross-section
# Bent 14 NR (WIDE (INTERLOCKING) SECTION)
ib = 14;
Dcol = 48.0*in;                                                         # Shorter, Width, of, octagonal, column (to, flat, sides)
ic = 4;                                                                                     # Single, Column
ColSecTag = 1000*ib+10*ic;                                   # Column's, fiber, model.element(tag., Follows, column, numbering, scheme.)
Hcol = [lindex, HcolList, 29];                                           # Column, Height
nLbar = [lindex, nLbarList, 29];                                           # Number, of, main (outer) longitudinal, bars
DLbar = [lindex, DLbarList, 29];                                           # Diameter, of, main (outer) longitudinal, bars
sTbar = [lindex, sTbarList, 29];                                           # Spacing, of, transverse, spiral, reinforcement
Wcol = 72.0*in;                                                         # Longer, Width, of, octagonal, column (to, flat, sides)
nLbar2 = 8;                                                                              # Number, of, secondary (inner) longitudinal, bars
DLbar2 = 0.625*in;                                                  # Diameter, of, secondary (inner) longitudinal, bars (#5, rebar)
colwidecmd, ColSecTag, Dcol, Wcol, nLbar, nLbar2, DLbar, DLbar2, sTbar; # Fiber, cross-section

