from math import sin, cos, pi
from opensees.section import FiberSection
from opensees import patch, layer

def _oct_outline(Rcol):
    n = 8
    phi =  2*pi/n
    R = Rcol/cos(phi/2)
    region = [
        layer.line(vertices=[
            [R*cos(i*phi-phi/2),  R*sin(i*phi-phi/2)],
            [R*cos(i*phi+phi/2),  R*sin(i*phi+phi/2)]
        ]) for i in range(n)
    ]
    sect = FiberSection(areas=region)
    sect.extRad = Rcol
    sect.intRad = 0.0
    return sect

def ConfiningPolygon(n, extRad=None, intRad=None, divs=None, diameter=None, s=4, material=None):
    psi = 2*pi/n
    phi = psi/s
    collection = []
    if divs is None:
        divs = 8,8 # divisions in each slice of the cover

    iR1, iR2 = [intRad]*2

    for i in range(n):
        startAngle =  (i - 1/2)*psi
        for j in range(s):
            sita1  =  startAngle + j*phi   # Slice start angle
            sita2  =  sita1 + phi          # Slice end angle
            oR1    =  extRad/cos(pi/n -  j*phi)
            oR2    =  extRad/cos(pi/n - (j+1)*phi)
            # Cover Patch connects the circular core to the polygonal cover
            collection.append(
              patch.quad(None, divs,
                vertices = [
                  [   iR1*cos(sita1),    iR1*sin(sita1)],
                  [   oR1*cos(sita1),    oR1*sin(sita1)],
                  [   oR2*cos(sita2),    oR2*sin(sita2)],
                  [   iR2*cos(sita2),    iR2*sin(sita2)],
                ]
            ))
    sect = FiberSection(areas=collection, material=material)
    sect.extRad = extRad
    sect.intRad = intRad
    return sect

def ConfinedPolygon(
    n:      int,
    extRad: float,
    intRad: float = None,
    sdivs : tuple = None, # sector divs
    cover : float = None,
    s:      int   = 4,    # number of patches per edge
    DLbar         = 4,
    core_conc     = None,
    cover_conc    = None,
    ColMatTag     = None,
    units         = None,
    diameter      = None
):
    """
    Dcol     :     Width of octagonal column (to flat sides)
    nLbar    :     Number of longitudinal bars
    DLbar    :     Diameter of longitudinal bars
    sTbar    :     Spacing of transverse spiral reinforcement
    """
    #
    # Column component dimensions
    #
    if intRad == extRad:
        return _oct_outline(extRad)

    if intRad is None:
        intRad = extRad - 2.0

    # assert intRad == 0.0 and extRad > 0.0

    inch    =  1.0
    Dcol    =  2*extRad
    # tcover  =  2.0*inch           # 2 inch cover width
    Rcol    =  Dcol/2.0           # Radius of octagonal column (to flat sides)
    Dcore   =  2*intRad
    Rcore   =  Dcore/2.0          # Radius of circular core
    # Along   =  pi*DLbar**2/4.0    # Area of longitudinal reinforcement bar
    DTbar   =  0.625*inch         # Diameter of transverse spiral bar (#5 Rebar)
    Asp     =  pi*DTbar**2/4.0    # Area of transverse spiral reinforcement bar
    Dtran   =  Dcore - DTbar      # Diameter of spiral of transverse spiral reinforcement

    # Density of transverse spiral reinforcement
    # rho       =  4.0* Asp/(Dtran*sTbar)
    # Diameter of ring of longitudinal reinforcement
    Dlong     =   Dcore - 2*DTbar - DLbar
    # Rlong     =   Dlong/2.        # Radius of ring of longitudinal reinforcement

    # Build Octagonal RC Column Section
    cdivs = 64        # # fibers around the entire circumference

    numSlices =  8    # # slices in each of the 8 sections of the polygon
    #if sdivs is None: sdivs = [1, 2]
    sect = FiberSection(
      material = ColMatTag,
      GJ       = 1e12,
      areas = [
        # Inner Core Patch, 5 radial fibers
        patch.circ(core_conc, [cdivs,  5], [0., 0.],     0.0, Rcore/2, 0.0, 2*pi),
        # Outer Core Patch, 10 radial fibers
        patch.circ(core_conc, [cdivs, 10], [0., 0.], Rcore/2, Rcore,   0.0, 2*pi)
    ])

    sect.add_patches(ConfiningPolygon(n, extRad, intRad, divs=sdivs, s=s).patches)

    #layer circ  long_steel  nLbar  Along 0. 0.  Rlong; # Longitudinal Bars
    sect.extRad = extRad
    sect.intRad = intRad
    return sect

