from math import sin, cos, pi, isclose
from opensees import section, patch
"""
This test validates section properties computed from assemblies
of patches layers and fibers.

  __ __________________ __
 4   |    |  b   |    |
  __ |    |______|    |
     | a  |      | c  | 16
     |    |      |    |
     |    |      |    |
     |____|      |____| __

     | 4  |  8   |  4 |


"""

def test_rectangles():
    sect = section.FiberSection(areas=[
        patch.rect(corners=[[-8,  0],[-4, 16]]),
        patch.rect(corners=[[-4, 12],[ 4, 16]]),
        patch.rect(corners=[[ 4,  0],[ 8, 16]]),
    ])

    assert sect.areas[0].ixc == 4*16**3/12
    assert sect.areas[1].ixc ==  8*4**3/12
    assert sect.areas[2].ixc == 4*16**3/12

    assert sect.centroid[1] == 9.2

    # centroid about x-x
    assert (14*8*4+8*2*16*4)/(16*8+4*8) == sect.centroid[1]

    #                           Ia        d      area         Ib
    assert sect.ixc == 2 * (4*16**3/12 + 1.2**2*(16*4)) + (8*4**3/12 + 4*8*4.8**2)

def test_polygon():
    import math
    Re, Ri = 20, 15
    for n in 4, 6, 7, 8:
        A = section.ConfinedPolygon(n,Re,s=10).area - section.ConfinedPolygon(n,Ri,s=10).area
        assert math.isclose(section.PolygonRing(n, Re, Ri).area, A)

def test_sector():
    r = 20.
    a = pi/7
    A = a*r**2
    ixc = ixo = r**4/4*(a - sin(a)*cos(a))
    iyo = r**4/4*(a + sin(a)*cos(a))
    yc = 2*r*sin(a)/(3*a)
    iyc = ixo + A*yc**2

    sect = patch.circ(extRad=r, startAng=pi/2-a, endAng=pi/2+a)

    assert isclose(sect.centroid[1], yc)
    assert isclose(sect.centroid[0], 0.)
    assert isclose(sect.area, A), f"{sect.area, A}"
    assert isclose(sect.ixc, ixc), f"{sect.ixc, ixc, iyc}"
    assert isclose(sect.iyc, iyc)




if __name__ == "__main__":
    test_rectangles()
    test_sector()
    test_polygon()

