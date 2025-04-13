# Claudio Perez
"""
A patch is used to generate a number of fibers over a cross-sectional area. 
Currently there are three types of patches that fibers can be generated over: 
quadrilateral, rectangular and circular.

All patches have the following attributes:

<dl>
 <dt><code>area</code></dt>
 <dd>Total area of the patch.</dd>

 <dt><code>moic</code></dt>
 <dd>Second moment of area matrix of the patch about its centroidal axis</dd>

 <dt><code>ixc</code></dt>
 <dd>Second moment of inertia of the patch about its $x$ axis</dd>

 <dt><code>iyc</code></dt>
 <dd>Second moment of inertia of the patch about its $y$ axis</dd>

</dl>

"""
import sys
import itertools
from opensees.library.ast import *
from opensees.library.ast import Grp, Num, Ref, Int
from opensees.library.obj import LibCmd, cmd
import numpy as np

class Material: pass
class Backbone: pass


@cmd
class Fiber:
    """
    This class represents a single fiber in an
    enclosing `FiberSection` or `NDFiberSection`.
    """
    tag_space = None
    # fiber $yLoc $zLoc $A $material
    _args = [
        Grp("coord", args=[Num("x"), Num("y")], reverse=True,
            about="$x$ and $y$ coordinate of the fiber in the section "\
                  "(local coordinate system)"),
        Num("area", about="area of the fiber."),
        Ref("material", type=Material,
            about="material tag associated with this fiber (UniaxialMaterial tag"\
                  "for a FiberSection and NDMaterial tag for use in an NDFiberSection)."),
    ]
    
    @property
    def fibers(self):
        yield self

#)

_patch = LibCmd("patch")


class _Polygon:
    def __init__(self, vertices):
        self.vertices = np.asarray(vertices)
        self._moic = None
        self._moig = None
        self._area = None
        self._fibers = None

    def __contains__(self, p:tuple) -> bool:
        v = self.vertices
        inside = (1 == sum(
            _rayIntersectSeg(p, (v[i-1], v[i])) for i in range(len(v)))%2
        )
        return inside

    @property
    def area(self):
        if self._area is None:
            x,y = np.asarray(self.vertices).T
            self._area = np.sum(x*np.roll(y, -1) - np.roll(x, -1)*y)*0.5
        return self._area

    @property
    def centroid(self):
        x,y = np.asarray(self.vertices).T
        area = self.area
        x_cent = (np.sum((x + np.roll(x, -1)) *
                         (x*np.roll(y, -1) -
                          np.roll(x, -1)*y)))/(6.0*area)
        y_cent = (np.sum((y + np.roll(y, -1)) *
                         (x*np.roll(y, -1) -
                          np.roll(x, -1)*y)))/(6.0*area)
        return np.array((x_cent, y_cent))

    def moi(self, reference, vertices=None, area=None):
        x,y = (np.asarray(self.vertices) - np.asarray(reference)).T

        area = area or self.area
        alpha = x * np.roll(y, -1) - np.roll(x, -1) * y

        # planar moment of inertia wrt horizontal axis
        ixx = np.sum((y**2 + y * np.roll(y, -1) +
                      np.roll(y, -1)**2)*alpha)/12.0

        # planar moment of inertia wrt vertical axis
        iyy = np.sum((x**2 + x * np.roll(x, -1) +
                      np.roll(x, -1)**2)*alpha)/12.0

        # product of inertia
        ixy = np.sum((x*np.roll(y, -1)
                      + 2.0*x*y
                      + 2.0*np.roll(x, -1) * np.roll(y, -1)
                      + np.roll(x, -1) * y)*alpha)/24.

        return np.array([[ixx, ixy],[ixy,iyy]])

    @property
    def moic(self):
        if self._moic is None:
            self._moic = self.moi(self.centroid)
        return self._moic

    @property
    def ixc(self):
        return self.moic[0,0]

    @property
    def iyc(self):
        return self.moic[1,1]

def RegularPolygon(n, Rcol=None, diameter=None):
    Rcol = Rcol or diameter/2
    phi =  2*np.pi/n
    R = Rcol/np.cos(phi/2)
    vertices = [
        [R*np.cos(i*phi-phi/2),  R*np.sin(i*phi-phi/2)]
        for i in range(n)
    ]
    poly = _Polygon(vertices)
    return poly

@_patch
class rect(_Polygon):
    _args = [
       Ref("material", type=Material, field="material",
           about="tag of previously defined material (`UniaxialMaterial`"\
                 "tag for a `FiberSection` or `NDMaterial` tag for use "\
                 "in an `NDFiberSection`)"),
       Grp("divs", reverse=True, args=[
           Int("ij", about="number of subdivisions (fibers) in the IJ direction."),
           Int("jk", about="number of subdivisions (fibers) in the JK direction."),
       ]),
       Grp("corners", args=[
         Grp(args=[Num("yI"), Num("zI")],  reverse=True, about="$y$ & $z$-coordinates of vertex I (local coordinate system)"),
         #Grp(args=[Num("yJ"), Num("zJ")],  reverse=True, about="$y$ & $z$-coordinates of vertex J (local coordinate system)"),
         Grp(args=[Num("yK"), Num("zK")],  reverse=True, about="$y$ & $z$-coordinates of vertex K (local coordinate system)"),
         #Grp(args=[Num("yL"), Num("zL")],  reverse=True, about="$y$ & $z$-coordinates of vertex L (local coordinate system)"),
      ])
    ]
    def init(self):
        self._moic = None
        self._moig = None
        self._area = None
        self._rule = "mid"
        self._interp = None
        self._fibers = None

        if len(self.corners) == 2:
            ll, ur = self.corners
            self.vertices = [ll, [ur[0], ll[1]], ur, [ll[0], ur[1]]]

    @property
    def fibers(self):
        if self._fibers is None:
            if self.divs is None:
                return []
            from shps.gauss import iquad
            # TODO: Use shps for lq4
            interp = self._interp or lq4
            rule = self._rule

            loc, wght = zip(iquad(rule=rule, n=self.divs[0]), iquad(rule=rule, n=self.divs[1]))

            x,y= zip(*(
                    interp(r,s)@self.vertices
                            for r,s in itertools.product(*loc)
            ))

            da = 0.25*np.fromiter(
                (dx*dy*self.area  for dx,dy in itertools.product(*wght)),
                float, len(x)
            )
            self._fibers = [
                Fiber([xi,yi], dai, self.material) for yi,xi,dai in sorted(zip(y,x,da))
            ]
        return self._fibers


@_patch
class quad(_Polygon):
    """A quadrilateral shaped patch.
    The geometry of the patch is defined by four vertices: I J K L.
    The coordinates of each of the four vertices is specified in *counter clockwise* sequence
    """
    _img  = "quadPatch.svg"
    _args = [
       Ref("material", type=Material, field="material",
           about="tag of previously defined material (`UniaxialMaterial` "\
                 "tag for a `FiberSection` or `NDMaterial` tag for use in an `NDFiberSection`)"),
       Grp("divs", reverse=True, args=[
         Int("ij", about="number of subdivisions (fibers) in the IJ direction."),
         Int("jk", about="number of subdivisions (fibers) in the JK direction."),
       ]),
       Grp("vertices", args=[
         Grp("i", args=[Num("x"), Num("y")],  about="$x$ & $y$-coordinates of vertex I (local coordinate system)", reverse=True),
         Grp("j", args=[Num("x"), Num("y")],  about="$x$ & $y$-coordinates of vertex J (local coordinate system)", reverse=True),
         Grp("k", args=[Num("x"), Num("y")],  about="$x$ & $y$-coordinates of vertex K (local coordinate system)", reverse=True),
         Grp("l", args=[Num("x"), Num("y")],  about="$x$ & $y$-coordinates of vertex L (local coordinate system)", reverse=True),
      ])
    ]

    def init(self):
        self._moic = None
        self._moig = None
        self._area = None
        self._interp = None
        self._rule = "mid"
        self._fibers = None

    @property
    def fibers(self):
        if self._fibers is None:
            from shps.gauss import iquad
            interp = self._interp or lq4
            rule = self._rule

            loc, wght = zip(iquad(rule=rule, n=self.divs[0]),
                            iquad(rule=rule, n=self.divs[1]))

            x,y= zip(*(
                    interp(r,s)@self.vertices
                            for r,s in itertools.product(*loc)
            ))

            da = 0.25*np.fromiter(
                (dx*dy*self.area  for dx,dy in itertools.product(*wght)), 
                float, len(x)
            )
            #y, x, da = map(list, zip(*sorted(zip(y, x, da))))
            self._fibers = [Fiber([xi,yi], dai, self.material) for yi,xi,dai in sorted(zip(y,x,da))]
        return self._fibers

def rhom(center, height, width, slope=None, divs=(0,0)):
    vertices = [
        [center[0] - width/2 + slope*height/2, center[1] + height/2],
        [center[0] - width/2 - slope*height/2, center[1] - height/2],
        [center[0] + width/2 - slope*height/2, center[1] - height/2],
        [center[0] + width/2 + slope*height/2, center[1] + height/2]
    ]
    return quad(vertices=vertices, div=divs)

@_patch
class circ:
    """
    Create a circular patch.

    ## Examples

    >>> patch.circ()
    """
    _signature = [
        "radius", "sector", "divs", "quadrature", "material",
        "center"
    ]
    _docsigs = [
    ]
    _img  = "circPatch.svg"
    _args = [
      Ref("material", type=Material, about="tag of previously defined material ("\
                          "`UniaxialMaterial` tag for a `FiberSection` "\
                          "or `NDMaterial` tag for use in an `NDFiberSection`)"),
      Grp("divs", args=[
          Int("circ", about="number of subdivisions (fibers) in "\
                            "the circumferential direction (number of wedges)"),
          Int("rad",  about="number of subdivisions (fibers) in the radial direction (number of rings)"),
      ]),
      Grp("center",    args=[Num("y"), Num("z")],
          about="$y$ & $z$-coordinates of the center of the circle", default=[0.0, 0.0]),
      Num("intRad",   about="internal radius", default=0.0),
      Num("extRad",   about="external radius"),
      Num("startAng", about="starting angle", default=0.0),
      Num("endAng",   about="ending angle", default=np.pi*2),
    ]
    def __contains__(self, point):
        origin = np.asarray(self.center)
        vect = np.asarray(point) - origin
        size = np.linalg.norm(vect)
        angl = (np.arctan2(*reversed(vect))+2*np.pi) % 2*np.pi
        inside = [
            (size < self.extRad),
            (size > self.intRad),
            (angl > self.startAng),
            (angl < self.endAng),
        ]
        return all(inside)

    def init(self):
        self._moic = None
        self._moig = None
        self._area = None
        self._fibers = None
        if "diameter" in self.kwds:
            self.extRad = self.kwds["diameter"]/2

    @property
    def fibers(self):
        """
        - https://stackoverflow.com/questions/33510979/generator-of-evenly-spaced-points-in-a-circle-in-python
        - https://mathworld.wolfram.com/CircleLatticePoints.html
        - Sunflower algorithm: https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
        - Sunflower alpha algorithm : https://stackoverflow.com/questions/28567166/uniformly-distribute-x-points-inside-a-circle
        - Partition of the circle in cells of equal area and shape: https://orbi.uliege.be/bitstream/2268/91953/1/masset_isocell_orbi.pdf
        """
        if self._fibers is None:
            if self.divs is None:
                self._fibers = []

            elif self.kwds.get("rule", "mid") == "mid":
                ri, ro = self.intRad, self.extRad
                dr = (self.extRad - self.intRad)/self.divs[1]
                dt = (self.endAng - self.startAng)/self.divs[0]
                areas = (0.5*dt*((ri+(i+1)*dr)**2 - (ri+i*dr)**2) for i in range(self.divs[1]))
                self._fibers = [
                    Fiber([r*np.cos(theta), r*np.sin(theta)], area, self.material)
                        for r,area in zip(np.linspace(self.intRad+dr/2, self.extRad-dr/2, self.divs[1]), areas)
                        #for theta in np.linspace(self.startAng+dt/2, self.endAng-dt/2, self.divs[0])
                            for theta in np.arange(self.startAng, self.endAng, dt)
                ]

            elif self.kwds.get("rule",None) == "uniform":
                self._fibers = [
                    Fiber([r*np.cos(theta), r*np.sin(theta)], area, self.material)
                    for r,theta,area in rtuniform(n=self.divs[1], rmax=self.extRad, m=self.divs[0], rmin=self.intRad)
                ]

            elif self.kwds.get("rule",None) == "uniform-2":
                self._fibers = [
                    Fiber([r*np.cos(theta), r*np.sin(theta)], area, self.material)
                    for r, theta, area in rtuniform(n=self.divs[1], rmax=self.extRad, m=self.divs[0])
                    if (r*np.cos(theta), r*np.sin(theta)) in self
                ]
            elif self.kwds.get("rule",None) == "uniform-3":
                self._fibers = [
                    Fiber([r*np.cos(theta), r*np.sin(theta)], area, self.material)
                    for r, theta, area in rtpairs(
                        np.linspace(self.intRad, self.extRad, self.divs[2]), 
                        np.arange(self.divs[1], self.divs[2]+self.divs[1]+1)*self.divs[0]
                    )
                    if (r*np.cos(theta), r*np.sin(theta)) in self
                ]
            elif self.kwds.get("rule",None) == "uniform-4":
                self._fibers = [
                    Fiber([r*np.cos(theta), r*np.sin(theta)], area, self.material)
                    for r, theta, area in rtpairs(
                        np.linspace(0.0, self.extRad, self.divs[2]),
                        np.arange(self.divs[1], self.divs[2]+self.divs[1]+1)*self.divs[0]
                    )
                    if (r*np.cos(theta), r*np.sin(theta)) in self
                ]
            elif self.kwds.get("rule",None) == "sunflower":
                area = self.area / self.divs[0]
                self._fibers = [
                    Fiber([x,y], area, self.material)
                    for x,y in sunflower(self.divs[0], self.extRad)
                    if (x,y) in self
                ]
            else:
                raise ValueError(f"Unknown quadrature rule, '{self.kwds['rule']}'.")
        return self._fibers

    @property
    def moic(self):
        if self._moic is None:
            r2,r1 = self.extRad, self.intRad
            a2,a1 = self.endAng, self.startAng
            dsin  = np.sin(2*a2) - np.sin(2*a1)
            self._moic = np.array([
                [(a2-a1 - 0.5*dsin), 0],
                [0, (a2-a1 + 0.5*dsin)],
            ])*(r2**4 - r1**4)*0.125 + self.area*self.centroid**2
        return self._moic

    @property
    def area(self):
        r2,r1 = self.extRad, self.intRad
        a2,a1 = self.endAng, self.startAng
        return 0.5*(r2**2 - r1**2)*(a2 - a1)

    @property
    def centroid(self):
        r2,r1 = self.extRad, self.intRad
        a2,a1 = self.endAng, self.startAng
        return np.array([
            np.sin(a2) - np.sin(a1),
            - np.cos(a2) + np.cos(a1)
        ])*(r2**3 - r1**3)/(3*self.area)


    @property
    def ixc(self):
        return self.moic[0,0]

    @property
    def iyc(self):
        return self.moic[1,1]

    @property
    def J(self):
        return 0.5 * self.area * self.extRad ** 2

layer = LibCmd("layer",
        {
          "circ": [
            Ref("material",  type=Material,
                 about="material tag of previously created material "\
                       "(UniaxialMaterial tag for a FiberSection or "\
                       "NDMaterial tag for use in an NDFiberSection)"),
            Int("divs",  about="number of fibers along arc"),
            Num("area", field="fiber_area", about="area of each fiber"),
            Grp("center", args=[Num("y"), Num("z")],
                about="$y$ and $z$-coordinates of center of circular arc"),
            Num("radius", about="radius of circular arc"),
            Grp("arc", args=[
              Num("startAng",  about="starting angle (optional, default = 0.0)"),
              Num("endAng",    about="ending angle (optional, default = 360.0 - 360/$numFiber)"),
            ])
          ]
        },
        about="The layer command is used to generate a number of fibers along a line or a circular arc.",
)

class ReinforcingLayer:
    ACIBar = {
        "3" : {"area": 0.11},
        "4" : {"area": 0.20},
        "5" : {"area": 0.31},
        "6" : {"area": 0.44},
        "7" : {"area": 0.60},
        "8" : {"area": 0.79},
        "9" : {"area": 1.00},
        "10": {"area": 1.27},
        "11": {"area": 1.56},
        "14": {"area": 2.25},
        "18": {"area": 2.25},
    }
    def init(self):
        self._moic = None
        self._fibers = None
        self._area = None
        if "bar" in self.kwds:
            assert self.fiber_area is None
            assert "units" in self.kwds
            self.fiber_area = ReinforcingLayer.ACIBar[self.kwds["bar"]]["area"]*self.kwds["units"].in2

    def moi(self, reference):
        xc = np.asarray(self.centroid)
        return  np.diag(
            sum(f.area*(f.coord - xc)**2 for f in self.fibers)
        )

    @property
    def area(self):
        if self._area is None:
            self._area = sum(f.area for f in self.fibers)
        return self._area

    @property
    def centroid(self):
        return sum(np.array(fiber.coord) * fiber.area for fiber in self.fibers)/self.area

    @property
    def moic(self):
        if self._moic is None:
            self._moic = self.moi(self.centroid)
        return self._moic

    @property
    def ixc(self):
        return self.moic[0,0]

    @property
    def iyc(self):
        return self.moic[1,1]

class _layer_namespace:
    @layer
    class circ(ReinforcingLayer):
        _args = [
                Ref("material",  type=Material,
                     about="material tag of previously created material "\
                           "(UniaxialMaterial tag for a FiberSection or "\
                           "NDMaterial tag for use in an NDFiberSection)"),
                Int("divs",  about="number of fibers along arc"),
                Num("area", field="fiber_area", about="area of each fiber"),
                Grp("center", args=[Num("y"), Num("z")],
                    about="$y$ and $z$-coordinates of center of circular arc"),
                Num("radius", about="radius of circular arc"),
                Grp("arc", default=[0., 2.*np.pi], args=[
                  Num("startAng",  about="starting angle"),
                  Num("endAng",    about="ending angle"),
                ])
        ]

        @property
        def fibers(self):
            if self._fibers is None:
                r = self.radius
                self._fibers = [Fiber([r*np.cos(t), r*np.sin(t)], self.fiber_area, self.material) 
                        for t in np.linspace(*self.arc, self.divs)]
            return self._fibers


@layer
class line(ReinforcingLayer):
    _img = "straightLayer.svg"
    _args = [
      Ref("material", type=Material, about="""Reference to previously created material 
                (`UniaxialMaterial` for a `FiberSection` or `NDMaterial` 
                for use in an `NDFiberSection`)"""),
      Int("divs", about="number of fibers along line"),
      Num("area", field="fiber_area", about="area of each fiber"),
      Grp("vertices", type=Grp, args=[
          Grp("start",args=[Num("x"), Num("y")], reverse=True,
              about="""$x$ and $y$-coordinates of first fiber
                       in line (local coordinate system)"""),
          Grp("end",  args=[Num("x"), Num("y")], reverse=True,
              about="$x$ and $y$-coordinates of last fiber in line (local coordinate system)")
      ])
    ]

    @property
    def fibers(self):
        if self._fibers is None:
            if self.divs is None:
                return []
            self._fibers = [Fiber([x, y], self.fiber_area, self.material) 
                for x,y in np.linspace(*self.vertices, self.divs)]
        return self._fibers

    def __contains__(self, point):
        a,b = np.asarray(self.vertices)
        p = np.asarray(point)
        return np.isclose(_distance(a,p) + _distance(p,b), _distance(a,b))

_eps = 0.00001
_huge = sys.float_info.max
_tiny = sys.float_info.min

def _rayIntersectSeg(p, edge)->bool:
    """
    takes a point p and an edge of two endpoints a,b of a line segment 
    and return bool

    https://rosettacode.org/wiki/Ray-casting_algorithm#Python
    """
    a,b = edge
    if a[1] > b[1]:
        a,b = b,a
    if p[1] == a[1] or p[1] == b[1]:
        p = np.array((p[0], p[1] + _eps))

    intersect = False

    if (p[1] > b[1] or p[1] < a[1]) or (p[0] > max(a[0], b[0])):
        return False

    if p[0] < min(a[0], b[0]):
        return True

    else:
        if abs(a[0] - b[0]) > _tiny:
            m_red = (b[1] - a[1]) / float(b[0] - a[0])
        else:
            m_red = _huge
        if abs(a[0] - p[0]) > _tiny:
            m_blue = (p[1] - a[1]) / float(p[0] - a[0])
        else:
            m_blue = _huge
        intersect = m_blue >= m_red
    return intersect


def _distance(a,b):
    return np.linalg.norm(a-b)


#
# Interpolation
#
def _Interpolant(docs, fwd, grad):
    fwd.grad = grad
    return fwd


lq4 = _Interpolant("",
    lambda r, s: 0.25*np.array([
            (r-1)*(s-1), -(r+1)*(s-1), (r+1)*(s+1), -(r-1)*(s+1)
    ]),
    lambda r, s: np.array([
            [s-1, -(s-1), (s+1), -(s+1)],
            [r-1, -(r+1), (r+1), -(r-1)]
    ])
)

lt6 = _Interpolant(
    """
    Quadratic Lagrange polynomial interpolation over a triangle.
    """,
    lambda r,s: np.array([
             (1 - r - s) - 2*r*(1 - r - s) - 2*s*(1 - r - s),
             r - 2*r*(1 - r - s) - 2*r*s,
             s - 2*r*s - 2*s*(1-r-s),
             4*r*(1 - r - s), #6
             4*r*s,
             4*s*(1 - r - s )
        ]),

    lambda r,s: np.array([
            [4*r + 4*s - 3, 4*r - 1, 0, -8*r - 4*s + 4, 4*s, -4*s],
            [4*r + 4*s - 3, 0, 4*s - 1, -4*r, 4*r, -4*r - 8*s + 4]])
)

def rtpairs(R,N):
    """
    R - list of radii
    N - list of points per radius

    Takes two list arguments containing the desired radii
    and the number of equally spread points per radii respectively.
    The generator, when iterated, will return radius-angle polar
    coordinate pairs, in metres and radians, which can be used 
    to plot shapes, e.g. a disc in the x-y plane. 
    """
    for i in range(len(R)-1):
        theta = 0.
        dTheta = 2*np.pi/N[i]
        for j in range(N[i]):
            theta = j*dTheta
            area = 0.5*dTheta*(R[i+1]**2 - R[i]**2)
            yield (R[i+1]+R[i])/2, theta, area


def rtuniform(n,rmax,m,rmin=0.0):
    """
    n - number of radii
    rmax - maximum radius
    m - scaling of points with radius

    This generator will return a disc of radius rmax, 
    with equally spread out points within it. The number 
    of points within the disc depends on the n and m parameters.
    """
    if not isinstance(n,int):
        n0, n = n
    else:
        n0 = 1
    R = [rmin]
    N = [n0]
    rmax_f = float(rmax)
    for i in range(int(n)):
        ri = rmin + (i+1)*((rmax_f-rmin)/int(n))
        ni = int(m)*(i+1)
        R.append(ri)
        N.append(ni)
    return rtpairs(R,N)


def sunflower(n, rad, alpha=0, geodesic=False):
    from math import sqrt, sin, cos, pi
    phi = (1 + sqrt(5)) / 2  # golden ratio
    def radius(k, n, b):
        if k > n - b:
            return 1.0
        else:
            return sqrt(k - 0.5) / sqrt(n - (b + 1) / 2)
    points = []
    angle_stride = 360 * phi if geodesic else 2 * pi / phi ** 2
    b = round(alpha * sqrt(n))  # number of boundary points
    for k in range(1, n + 1):
        r = rad*radius(k, n, b)
        theta = k * angle_stride
        points.append((r * cos(theta), r * sin(theta)))
    return points

