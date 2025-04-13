from math import cos,sin,sqrt,pi
import opensees as ops
#
# COLUMN PINS -------------------------------------------------------------------------------
#
rFlag = 1;  # Rayleigh, damping, ON
model.uniaxialMaterial('Elastic', 21             , 1.0e)+15*kips/in                      # RIGID, 'TRANSLATIONAL', AND, 'TORSIONAL', 'STIFFNESS'
model.uniaxialMaterial('Elastic', 22       , 0.0*kips/in                             )# FREE, X 'and/or', Y, 'ROTATIONAL', 'STIFFNESS'
# Procedure for UNREINFORCED column fiber cross-section assignment (INELASTIC material properties):
# source ColSectionOctIEPin.tcl
source, 'ColSectionLib.tcl'
Dcol =, 84.0*in
BuildOctColPINSection, 290000, Dcol; # Fiber, cross-section, 'for', columns, 'at', Bents, 2-11 (each, 'have', diameter, 84.0, inches)
yp1 =        [lindex, vecxzXList, 0];
yp2 =        [lindex, vecxzYList, 0];
model.element('zeroLengthSection', 290000, 201, 2010, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 2, Left, Column, Base, DIR = {}., Modeled, 'as', a, 'zerolength', octagonal, 'section', of, 'unreinforced', concrete., Orientation, 'vector', s = 'the', local, x 'axis', (perpendicular, 'to', the, section) in, 'the', global, Z (vertical) direction, and, 'the', local, y 'axis', parallel, 'to', the, 'length', of, 'the', cap, 'beam.'
model.element('zeroLengthSection', 290001, 207, 2070, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 2, Right, Column, Base, DIR = {Y}
yp1 =        [lindex, vecxzXList, 1];
yp2 =        [lindex, vecxzYList, 1];
model.element('zeroLengthSection', 390000, 302, 3020, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 3, Left, Column, Top, DIR = {X}
model.element('zeroLengthSection', 390001, 306, 3060, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 3, Right, Column, Top, DIR = {X,Y}
yp1 =        [lindex, vecxzXList, 2];
yp2 =        [lindex, vecxzYList, 2];
model.element('zeroLengthSection', 490000, 402, 4020, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 4, Left, Column, Top, DIR = {X}
model.element('zeroLengthSection', 490001, 406, 4060, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 4, Right, Column, Top, DIR = {X,Y}
yp1 =        [lindex, vecxzXList, 3];
yp2 =        [lindex, vecxzYList, 3];
model.element('zeroLengthSection', 590000, 502, 5020, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 5, Left, Column, Top, DIR = {Y}
model.element('zeroLengthSection', 590001, 506, 5060, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 5, Right, Column, Top, DIR = {X}
yp1 =        [lindex, vecxzXList, 4];
yp2 =        [lindex, vecxzYList, 4];
model.element('zeroLengthSection', 690000, 602, 6020, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 6, Left, Column, Top, DIR = {X}
model.element('zeroLengthSection', 690001, 606, 6060, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 6, Right, Column, Top, DIR = {X}
yp1 =        [lindex, vecxzXList, 5];
yp2 =        [lindex, vecxzYList, 5];
model.element('zeroLengthSection', 790000, 702, 7020, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 7, Left, Column, Top, DIR = {X}
model.element('zeroLengthSection', 790001, 706, 7060, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 7, Right, Column, Top, DIR = {X,Y}
yp1 =        [lindex, vecxzXList, 6];
yp2 =        [lindex, vecxzYList, 6];
model.element('zeroLengthSection', 890000, 802, 8020, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 8, Left, Column, Top, DIR = {X}
model.element('zeroLengthSection', 890001, 806, 8060, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 8, Right, Column, Top, DIR = {X,Y}
yp1 =        [lindex, vecxzXList, 7];
yp2 =        [lindex, vecxzYList, 7];
model.element('zeroLengthSection', 990000, 902, 9020, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 9, Left, Column, Top, DIR = {X}
model.element('zeroLengthSection', 990001, 906, 9060, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 9, Right, Column, Top, DIR = {X}
yp1 =        [lindex, vecxzXList, 8];
yp2 =        [lindex, vecxzYList, 8];
model.element('zeroLengthSection', 1090000, 1002, 10020, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 10, Left, Column, Top, DIR = {X,Y}
model.element('zeroLengthSection', 1090001, 1006, 10060, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 10, Right, Column, Top, DIR = {X,Y}
yp1 =        [lindex, vecxzXList, 9];
yp2 =        [lindex, vecxzYList, 9];
model.element('zeroLengthSection', 1190000, 1102, 11020, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 11, Left, Column, Top, DIR = {X}
model.element('zeroLengthSection', 1190001, 1106, 11060, 290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 11, Right, Column, Top, DIR = {X}
Dcol =, 66.0*in
BuildOctColPINSection, 1290000, Dcol; # Fiber, cross-section, 'for', columns, 'at', Bent, 12 (each, 'have', diameter, 66.0, inches) 
yp1 =        [lindex, vecxzXList, 10];
yp2 =        [lindex, vecxzYList, 10];
model.element('zeroLengthSection', 1290000, 1201, 12010, 1290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 12, Left, Column, Base, DIR = {Y}
model.element('zeroLengthSection', 1290001, 1209, 12090, 1290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 12, Right, Column, Base, DIR = {Y}
model.element('zeroLengthSection', 1290002, 1211, 12110, 1290000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 12, Center, Column, Base, DIR = {Y}
Dcol =, 48.0*in
BuildOctColPINSection, 1390000, Dcol; # Fiber, cross-section, 'for', columns, 'at', Bents, 13-14 (each, 'have', diameter, 48.0, inches)
yp1 =        [lindex, vecxzXList, 11];
yp2 =        [lindex, vecxzYList, 11];
model.element('zeroLengthSection', 1390000, 1301, 13010, 1390000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 13, Left, Column, Base, DIR = {X}
model.element('zeroLengthSection', 1390001, 1307, 13070, 1390000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 13, Right, Column, Base, DIR = {X,Y}
yp1 =        [lindex, vecxzXList, 12];
yp2 =        [lindex, vecxzYList, 12];
model.element('zeroLengthSection', 1490000, 1401, 14010, 1390000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 14, Left, Column, Base, DIR = {}
model.element('zeroLengthSection', 1490001, 1407, 14070, 1390000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 14, Right, Column, Base, DIR = {Y}
model.element('zeroLengthSection', 1490002, 1409, 14090, 1390000 -orient, 0 0, 1 yp1, yp2, 0 -doRayleigh, rFlag); # PIN, 'Bent', 14, Center, Column, Base, DIR = {Y}

