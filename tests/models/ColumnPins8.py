from math import cos,sin,sqrt,pi
import opensees as ops
#
# COLUMN PINS -------------------------------------------------------------------------------
#
rFlag = 1;  # Rayleigh, damping, ON
model.uniaxialMaterial('Elastic', 21             , 1.0e)+15*kips/in                      # RIGID, 'TRANSLATIONAL', AND, 'TORSIONAL', 'STIFFNESS'
model.uniaxialMaterial('Elastic', 22       , 0.0*kips/in                             )# FREE, X 'and/or', Y, 'ROTATIONAL', 'STIFFNESS'
model.element('zeroLength',        290000, 201, 2010 -mat, 21, 21, 21, 21, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 2, Left, Column, Base, DIR = {}
model.element('zeroLength',        290001, 207, 2070 -mat, 21, 21, 21, 21, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 2, Right, Column, Base, DIR = {Y}
model.element('zeroLength',        390000, 302, 3020 -mat, 21, 21, 21, 22, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 3, Left, Column, Top, DIR = {X}
model.element('zeroLength',        390001, 306, 3060 -mat, 21, 21, 21, 22, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 3, Right, Column, Top, DIR = {X,Y}
model.element('zeroLength',        490000, 402, 4020 -mat, 21, 21, 21, 22, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 4, Left, Column, Top, DIR = {X}
model.element('zeroLength',        490001, 406, 4060 -mat, 21, 21, 21, 22, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 4, Right, Column, Top, DIR = {X,Y}
model.element('zeroLength',        590000, 502, 5020 -mat, 21, 21, 21, 21, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 5, Left, Column, Top, DIR = {Y}
model.element('zeroLength',        590001, 506, 5060 -mat, 21, 21, 21, 22, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 5, Right, Column, Top, DIR = {X}
model.element('zeroLength',        690000, 602, 6020 -mat, 21, 21, 21, 22, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 6, Left, Column, Top, DIR = {X}
model.element('zeroLength',        690001, 606, 6060 -mat, 21, 21, 21, 22, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 6, Right, Column, Top, DIR = {X}
model.element('zeroLength',        790000, 702, 7020 -mat, 21, 21, 21, 22, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 7, Left, Column, Top, DIR = {X}
model.element('zeroLength',        790001, 706, 7060 -mat, 21, 21, 21, 22, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 7, Right, Column, Top, DIR = {X,Y}
model.element('zeroLength',        890000, 802, 8020 -mat, 21, 21, 21, 22, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 8, Left, Column, Top, DIR = {X}
model.element('zeroLength',        890001, 806, 8060 -mat, 21, 21, 21, 22, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 8, Right, Column, Top, DIR = {X,Y}
model.element('zeroLength',        990000, 902, 9020 -mat, 21, 21, 21, 22, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 9, Left, Column, Top, DIR = {X}
model.element('zeroLength',        990001, 906, 9060 -mat, 21, 21, 21, 22, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 9, Right, Column, Top, DIR = {X}
model.element('zeroLength',        1090000, 1002, 10020 -mat, 21, 21, 21, 22, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 10, Left, Column, Top, DIR = {X,Y}
model.element('zeroLength',        1090001, 1006, 10060 -mat, 21, 21, 21, 22, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 10, Right, Column, Top, DIR = {X,Y}
model.element('zeroLength',        1190000, 1102, 11020 -mat, 21, 21, 21, 22, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 11, Left, Column, Top, DIR = {X}
model.element('zeroLength',        1190001, 1106, 11060 -mat, 21, 21, 21, 22, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 11, Right, Column, Top, DIR = {X}
model.element('zeroLength',        1290000, 1201, 12010 -mat, 21, 21, 21, 21, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 12, Left, Column, Base, DIR = {Y}
model.element('zeroLength',        1290001, 1209, 12090 -mat, 21, 21, 21, 21, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 12, Right, Column, Base, DIR = {Y}
model.element('zeroLength',        1290002, 1211, 12110 -mat, 21, 21, 21, 21, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 12, Center, Column, Base, DIR = {Y}
model.element('zeroLength',        1390000, 1301, 13010 -mat, 21, 21, 21, 22, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 13, Left, Column, Base, DIR = {X}
model.element('zeroLength',        1390001, 1307, 13070 -mat, 21, 21, 21, 22, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 13, Right, Column, Base, DIR = {X,Y}
model.element('zeroLength',        1490000, 1401, 14010 -mat, 21, 21, 21, 21, 21, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 14, Left, Column, Base, DIR = {}
model.element('zeroLength',        1490001, 1407, 14070 -mat, 21, 21, 21, 21, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 14, Right, Column, Base, DIR = {Y}
model.element('zeroLength',        1490002, 1409, 14090 -mat, 21, 21, 21, 21, 22, 21 -dir, 1 2, 3 4, 5 6 -doRayleigh, rFlag); # PIN, 'Bent', 14, Center, Column, Base, DIR = {Y}

