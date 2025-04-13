
from math import cos,sin,sqrt,pi
import numpy as np
import opensees as ops
from opensees.units.english import inch

class HaywardSections:
    def __init__(self, base_type, wide_type):
        self._base_type = base_type
        self._wide_type = wide_type
        self.sections = {}

    def add_basic(self, tag, *args): # diameter, nlong, dlong, stran):
        if args in self.sections:
            return self.sections[args]

    def assign(self):
       #
       # COLUMN CROSS-SECTION ASSIGNMENT ------------------------------
       #
       # (Column Numbering scheme: 1000*Bent#+10*Column#)
       # (Column#: 1=Left, 2=Right, 3=Center, 4=Single)
       # Column fiber element tag (ColSecTag) follows column numbering scheme.
       # global 'kips', in 'ksi_psi', pi ;

       # Read in the number of longitudinal bars for each column:
       # Read in the diameter of longitudinal bars for each column:
       # Read in the spacing of transverse spiral reinforcement for each column:
       bar_list = np.loadtxt("columns.txt", delimiter="\t")

       # BUILD COLUMN SECTIONS
       # Bents 2-11
       Dcol = 84.0*inch
       for ib in range(2, 11+1):
           DLbar, nLbar, sTbar = bar_list[(ib-2)*2]
           for ic in 1, 2: # Left and right columns
               base_type(1000*ib+10*ic, Dcol, nLbar, DLbar, sTbar); # Fiber cross-section

       # Bent 12
       ib = 12;
       Dcol = 66.0*inch
       ic = 1                                                       # Left Column
       DLbar, nLbar, sTbar = bar_list[20] # 20, 20, 20
       base_type(1000*ib+10*ic, Dcol, nLbar, DLbar, sTbar); # Fiber cross-section
       ic = 2                                                       # Right Column
       DLbar, nLbar, sTbar = bar_list[22] # 22, 22, 22
       base_type(1000*ib+10*ic, Dcol, nLbar, DLbar, sTbar); # Fiber cross-section
       ic = 3                                                       # Center Column
       DLbar, nLbar, sTbar = bar_list[21] # 21, 21, 21
       base_type(1000*ib+10*ic, Dcol, nLbar, DLbar, sTbar); # Fiber cross-section

       # Bent 13 NE
       ib = 13;
       Dcol = 48.0*inch
       ic = 1                                         # Left Column
       DLbar, nLbar, sTbar = bar_list[23] # 23, 23, 23
       base_type(1000*ib+10*ic, Dcol, nLbar, DLbar, sTbar); # Fiber cross-section

       ic = 2                                         # Right Column
       DLbar, nLbar, sTbar = bar_list[24] # 24, 24, 24
       base_type(1000*ib+10*ic, Dcol, nLbar, DLbar, sTbar); # Fiber cross-section

       # Bent 13 NR (WIDE (INTERLOCKING) SECTION)
       ib = 13;
       Dcol = 48.0*inch;                # Shorter Width of octagonal column (to flat sides)
       ic = 4                           # Single Column
       DLbar, nLbar, sTbar = bar_list[25] # 25, 25, 25
       Wcol = 72.0*inch;                # Longer Width of octagonal column (to flat sides)
       nLbar2 = 8                       # Number of secondary (inner) longitudinal bars
       DLbar2 = 0.625*inch;             # Diameter of secondary (inner) longitudinal bars (#5 rebar)
       wide_type(1000*ib+10*ic, Dcol, Wcol, nLbar, nLbar2, DLbar, DLbar2, sTbar); #, Fiber, cross-section

       # Bent 14 NE
       ib = 14;
       Dcol = 48.0*inch
       ic = 1                                                       # Left Column
       DLbar, nLbar, sTbar = bar_list[26] # 26, 26, 26
       base_type(1000*ib+10*ic, Dcol, nLbar, DLbar, sTbar); # Fiber cross-section
       ic = 2                                                       # Right Column
       DLbar, nLbar, sTbar = bar_list[28] # 28, 28, 28
       base_type(1000*ib+10*ic, Dcol, nLbar, DLbar, sTbar); # Fiber cross-section
       ic = 3                                                       # Center Column
       DLbar, nLbar, sTbar = bar_list[27] # 27, 27, 27
       base_type(1000*ib+10*ic, Dcol, nLbar, DLbar, sTbar); # Fiber cross-section
       # Bent 14 NR (WIDE (INTERLOCKING) SECTION)
       ib = 14;
       Dcol = 48.0*inch;                       # Shorter Width of octagonal column (to flat sides)
       ic = 4                                  # Single Column
       DLbar, nLbar, sTbar = bar_list[29] # 29, 29, 29
       Wcol = 72.0*inch;                       # Longer Width of octagonal column (to flat sides)
       nLbar2 = 8                              # Number of secondary (inner) longitudinal bars
       DLbar2 = 0.625*inch;                    # Diameter of secondary (inner) longitudinal bars (#5 rebar)
       wide_type(1000*ib+10*ic, Dcol, Wcol, nLbar, nLbar2, DLbar, DLbar2, sTbar); #, Fiber, cross-section

AssignHaywardSections(print, print)

