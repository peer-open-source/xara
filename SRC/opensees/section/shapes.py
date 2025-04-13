from . import patch
layer = patch.layer
from opensees.section import FiberSection
import numpy as np
from numpy import pi, sin, cos

def RegularPolygon():
    pass

def WideFlange():
    pass


def PolygonRing(n, extRad, intRad):
    """
    Create a polygon annulus.
    """
    psi = 2*pi/n
    phi = psi
    collection = []
    cover_divs = 1,2      # divisions in each slice of the cover
    iR1, iR2 = [intRad/cos(pi/n)]*2
    oR1, oR2 = [extRad/cos(pi/n)]*2
    j = 0
    for i in range(n):
        startAngle =  (i - 1/2)*psi
        sita1  =  startAngle + j*phi   # Slice start angle
        sita2  =  sita1 + phi          # Slice end angle
        # Cover Patch connects the circular core to the polygonal cover
        collection.append(
          patch.quad(None, cover_divs,
            vertices = [
              [   iR1*cos(sita1),    iR1*sin(sita1)],
              [   oR1*cos(sita1),    oR1*sin(sita1)],
              [   oR2*cos(sita2),    oR2*sin(sita2)],
              [   iR2*cos(sita2),    iR2*sin(sita2)],
            ]
        ))
    sect = FiberSection(shapes=collection)
    sect.extRad = extRad
    sect.intRad = intRad
    return sect


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
    sect = FiberSection(shapes=region)
    sect.extRad = Rcol
    sect.intRad = 0.0
    return sect

def ConfiningPolygon(n, extRad=None, intRad=None, divs=None, diameter=None, s=4, material=None):
    psi = 2*pi/n
    phi = psi/s
    collection = []
    if divs is None: divs = 8,8 # divisions in each slice of the cover
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
    sect = FiberSection(shapes=collection, material=material)
    sect.extRad = extRad
    sect.intRad = intRad
    return sect

def ConfinedPolygon(
    n:      int,
    extRad: float,
    intRad: float = None,
    sdivs : tuple = None, # sector divs
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
    Dlong   =   Dcore - 2*DTbar - DLbar
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




def sect2shapely(section):
    """
    Generate `shapely` geometry objects
    from `opensees` patches or a FiberSection.
    """
    import numpy as np
    import shapely.geometry
    from shapely.ops import unary_union
    shapes = []
    if hasattr(section, "patches"):
        patches = section.patches
    else:
        patches = [section]

    for patch in patches:
        name = patch.__class__.__name__.lower()
        if name in ["quad", "poly", "rect", "_polygon"]:
            points = np.array(patch.vertices)
            width,_ = points[1] - points[0]
            _,height = points[2] - points[0]
            shapes.append(shapely.geometry.Polygon(points))
        else:
            n = 64
            x_off, y_off = 0.0, 0.0
            # calculate location of the point
            external = [[
                0.5 * patch.extRad * np.cos(i*2*np.pi*1./n - np.pi/8) + x_off,
                0.5 * patch.extRad * np.sin(i*2*np.pi*1./n - np.pi/8) + y_off
                ] for i in range(n)
            ]
            if patch.intRad > 0.0:
                internal = [[
                    0.5 * patch.intRad * np.cos(i*2*np.pi*1./n - np.pi/8) + x_off,
                    0.5 * patch.intRad * np.sin(i*2*np.pi*1./n - np.pi/8) + y_off
                    ] for i in range(n)
                ]
                shapes.append(shapely.geometry.Polygon(external, [internal]))
            else:
                shapes.append(shapely.geometry.Polygon(external))

    if len(shapes) > 1:
        return unary_union(shapes)
    else:
        return shapes[0]

def sect2gmsh(sect, size, **kwds):
    import pygmsh
    import numpy as np
    if isinstance(size, int):
        size = [size]*2

    shape = sect2shapely(sect)

    with pygmsh.geo.Geometry() as geom:
        geom.characteristic_length_min = size[0]
        geom.characteristic_length_max = size[1]
        coords = np.array(shape.exterior.coords)
        holes = [
            geom.add_polygon(np.array(h.coords)[:-1], size[0], make_surface=False).curve_loop
            for h in shape.interiors
        ]
        if len(holes) == 0:
            holes = None

        poly = geom.add_polygon(coords[:-1], size[1], holes=holes)
        # geom.set_recombined_surfaces([poly.surface])
        mesh = geom.generate_mesh(**kwds)

    mesh.points = mesh.points[:,:2]
    for blk in mesh.cells:
        blk.data = blk.data.astype(int)
    # for cell in mesh.cells:
    #     cell.data = np.roll(np.flip(cell.data, axis=1),3,1)
    return mesh


class TorsionalConstantAnalysis:
    pass

class MomentCurvatureAnalysis:
    @staticmethod
    def solve_eps(sect, kap, axial: float, eps0, tol=1e-6, maxiter=25):
        # Newton-Raphson iteration
        eps = eps0
        s = sect.getStressResultant([eps, kap], False)
        for i in range(maxiter):
            if abs(s[0] - axial) < tol:
                return eps
            s = sect.getStressResultant([eps, kap], False)
            eps -= (s[0] - axial)/sect.getSectionTangent()[0,0]
        return eps

    def __init__(self, axial):
        pass

def MomentSearch(section, ):
    pass

class MomentAxialLocus:
    def __init__(self, section, axial):
        self.axial = axial
        self.section = section

    def plot(self):
        pass

    def analyze(self, nstep = 30, incr=5e-6):
        import matplotlib.pyplot as plt
        import numpy as np
        fig, ax = plt.subplots(1,2, constrained_layout=True)
        sect = self.section
        axial = self.axial

        if sect.name is None:
            sect.name = 1

        solve_eps = MomentCurvatureAnalysis.solve_eps

        # Curvature increment
        dkap = incr
        for P in axial:
            with sect as s:
                k0 = 0.0
                e0 = solve_eps(s,  k0,  P,  solve_eps(s,  k0,  P,  0.0))
                PM = [
                    s.getStressResultant([e0, k0], True),
                    s.getStressResultant([solve_eps(s, k0+dkap, P, e0), k0+dkap], True),
                ]
                e = e0
                kap = 2*dkap
                for _ in range(nstep):
                    if abs(PM[-1][1]) < 0.995*abs(PM[-2][1]):
                        break
                    e = solve_eps(s, kap, P, e)
                    PM.append(s.getStressResultant([e, kap], True))
                    kap += dkap

            p, m = zip(*PM)

            ax[0].plot(np.linspace(0.0, kap, len(m)), m)

            ax[1].scatter(m, p, s=0.2, color="k")

        ax[1].set_ylabel("Axial force, $P$")
        ax[1].set_xlabel("Moment, $M$")

        plt.show()

