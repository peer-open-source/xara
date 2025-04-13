from .ast import Num, Tag, Ref, Blk, Int, Lst, Grp, Str
from .obj import LibCmd, cmd
from opensees.library import Dof

_pattern = LibCmd("pattern")
_series  = LibCmd("timeSeries")


@_series
class TimeSeries:
    tag_space = "timeSeries"
    _name = "Path"
    _args = [
        Tag(),
        Lst("time",   flag="-time",   type=float),
        Lst("values", flag="-values", type=float),
        Num("scale",  flat="-fact",   reqd=False)
    ]

    def init(self, series=None, step=None, time=None):
        if hasattr(self.values, "time"):
            self.time = self.values.time
            self.values = self.values.data
        if "step" in self.kwds:
            dt = self.time_step = self.kwds["step"]
            self.time = [0.0 + i*dt for i in range(len(self.values))]
        else:
            self.time_step = None

@cmd
class GroundMotion:
    _args = [
         Tag(),
         Str("motion_type", default="Plain"),
         Ref("accel", flag="-accel", reqd=False, type=TimeSeries),
         Ref("veloc", flag="-vel",   reqd=False, type=TimeSeries),
         Ref("displ", flag="-disp",  reqd=False, type=TimeSeries)
         # <-int (IntegratorType intArgs)> 
         # <-fact $cFactor>
    ]
    # _alts = [
    #     Ref("motion")
    # ]
    _refs = ["accel", "veloc", "displ"]


    def init(self):
        kwds = {}
        if "time" in self.kwds:
            kwds["time"] = self.kwds["time"]

        m = self.kwds.get("motion", None)

        accel = getattr(m, "accel", getattr(self, "accel", None))
        if accel is not None and not isinstance(accel, TimeSeries):
            self.accel = TimeSeries(values=accel, **kwds)

        veloc = getattr(m, "veloc", getattr(self, "veloc", None))
        if veloc is not None and not isinstance(veloc, TimeSeries):
            self.veloc = TimeSeries(values=veloc, **kwds)

        displ = getattr(m, "displ", getattr(self, "displ", None))
        if displ is not None and not isinstance(displ, TimeSeries):
            self.displ = TimeSeries(values=displ, **kwds)


@cmd
class ImposedMotion:
    _args = [
        Ref("node", type="node"),
        Dof(),
        Ref("motion", type=GroundMotion)
    ]
    _refs = ["node"]



load = cmd("load", "load", args=[Ref("node"), Grp("load", min=1, type=Num)])


@_pattern
class Plain:
    _args = [
        Tag(),
        # Ref("series", type=TimeSeries, about="Reference to the time series to be used in the load pattern"),
        Str("series"),
        Blk("loads"),
        Num("scale", reqd=False, about="constant scale factor"),
    ]
    _refs = ["series"]

    def init(self):
        if isinstance(self.series, str):
            pass

        elif not isinstance(self.series, (TimeSeries,int)):
            self.series = TimeSeries(values=self.series)

        if isinstance(self.loads, dict):
            loads = []
            for k,v in self.loads.items():
                if isinstance(k, int):
                    loads.append(load(k, v))
                else:
                    loads.append(load(k[0], [v]))

            self.loads = loads


@_pattern
class UniformExcitation:
    _args = [
        Tag(),
        Dof(),
        Ref("accel", flag="-accel", type=TimeSeries, attr="name", about="acceleration time series"),
        Num("v0",    reqd=False, about="initial velocity."),
        Num("scale", flag="-fact", reqd=False)
    ]
    _refs = ["accel"]

    @property
    def series(self): return self.accel

    def init(self):
        if not isinstance(self.accel, int):
            self.accel = TimeSeries(values=self.accel)

@_pattern
class MultipleSupport:
    _args = [
        Tag(),
        Blk("motions", alt="components")
    ]
    """
    pattern.MultipleSupport(
        components = [
        #   node, dof,  history
            ( 1,   1,   ResponseComponent(accel, displ, veloc)),
            ( 2,   1,   ResponseComponent(accel, displ, veloc))
        ]
    )
    """
    def init(self):
        imposed_motions = []
        ground_motions = []

        for m in self.motions:
            if isinstance(m[2], GroundMotion):
                ground_motions.append(m[2])
            else:
                ground_motions.append(GroundMotion(motion=m[2]))
            imposed_motions.append(ImposedMotion(*m[:2], ground_motions[-1]))
        self.motions = ground_motions + imposed_motions

    @property
    def components(self):
        return []


