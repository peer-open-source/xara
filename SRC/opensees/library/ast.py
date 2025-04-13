
class Parameter:
    def __init__(self, name, typ):
        self.name = name
        self.typ = typ
    def __repr__(self):
        return f"${self.name}"

class Arg:
    __slots__ = ["name", "flag", "value", "field", "default", "type", 
                 "reqd", "ndm_reqd", "ndm_only",
                 "namespace", "kwds", "about"]#, "group"]
    def __init__(self,
        name = None,
        #help = None,
        flag: str  = None,
        reqd: bool = True,
        ndm_reqd   = (),
        ndm_only   = (),
        type       = None,
        field: str = None,
        about      = "",
        default    = None,
        # group    = None,
        **kwds
    ):
        if name and name[0] == "-":
            flag = name
            name = name[1:]
        self.name      = name
        self.flag      = flag if flag is not None else ""
        self.value     = None
        self.field     = field if field is not None  else name
        self.type      = type
        self.reqd      = reqd
        self.ndm_reqd  = ndm_reqd
        self.ndm_only  = ndm_only
        self.default = default
        self.kwds = kwds
        self.about = about
        # self.group = group
        self.init()

    def __repr__(self):
        return self.__class__.__name__ + f"<{self.name}>"

    def init(self): pass

    def get_value(self, parent, fmt):
        val = getattr(parent, self.field, None)
        if isinstance(val,bool):
            return str(val).lower()
        else:
            return val

    def _get_value(self, parent, value=None):
        value = self.value if value is None else value
        if value is None:
            value = self.default
        return value

    def as_tcl_list(self, value=None)->list:
        value = self._get_value(None, value)
        if isinstance(value, (type(None), )):
            if not self.reqd:
                return []
            else:
                value = f"${self.name}"
        elif isinstance(value, (Arg,)):
            if not self.reqd:
                return []
            else:
                value = f"${value.name}"

        return [self.flag] + [value]

    def m_src(self, value=None):
        value = self._get_value(None,value)
        if isinstance(value,bool):
            return str(value).lower()

class Flg(Arg):
    "A `Flg` is a flag-like argument"
    def init(self):
        if "enum" in self.kwds:
            self.enum = self.kwds["enum"]
            if self.type is None:
                self.type = str

        self.default = False
        if "-" in self.name:
            if self.field == self.name:
                self.field = self.name[1:]
            self.flag = self.name
            self.name = self.name[1:]

    def c_read_argv(self, struct, argidx):
        return f"{struct}->{self.field} = true;\n"

    def as_tcl_list(self, value=None):
        value = self.value if value is None else value
        return [self.flag] if value else []

class Num(Arg):
    def init(self): self.type = float

    def c_read_argv(self, struct, argidx):
        return f"{struct}->{self.field} = argv[++{argidx}];\n"

class Int(Arg):
    def init(self): self.type = int

class Str(Arg):
    def init(self): self.type = str

class Chr(Arg):
    def init(self): self.type = chr

class Ref(Arg):
    """A `Ref` is a reference to another object.

    ## Extra parameters

    `attr="name"`: referenced attribute

    """
    def init(self):
        if "attr" not in self.kwds: self.kwds["attr"] = "name"
        self.__name__ = f"Ref({self.type.__class__.__name__})"

    def __call__(self, *args,**kwds):
        return Ref(*args, type=self.type, attr=self.kwds["attr"], **kwds)

    def _get_value(self, parent, value=None):
        value = self.value if value is None else value
        if isinstance(value, dict):
            return self.kwds["typ"](**value)._get_value(self)
        elif not isinstance(value, (str,int,float,type(None))):
            return getattr(value, self.kwds["attr"])
        else:
            return value

    def as_tcl_list(self, value=None): 
        value = self._get_value(None, value)
        try:
            value = getattr(value, self.kwds["attr"])
        except:
            pass
        return [self.flag] + [value]

class Lst(Arg): pass
class Ary(Arg): pass

class Sub(Arg): pass

class Map(Arg):
    def  init(self):
        self.key_type = self.kwds["key"]
        self.val_type = self.kwds["val"]
        self.tcl_rev_kv = self.kwds.get("tcl_rev_kv", False)
        #self.__name__ = f"Map({self.key_type.__class__.__name__} : {self.val_type.__class__.__name__})"
        self.__name__ = f"{{{self.key_type.name} : {self.val_type.name} ...}}"

    def as_tcl_list(self, value=None):
        value = self._get_value(None,value)
        if value is None:
            if self.reqd:
                value = [None]*self.num
            else:
                return []
        vals = []
        i, j = (1, 0) if self.tcl_rev_kv else (0, 1)
        for arg, kv in zip(arg,value.items()):
            kv = kv[0], arg.as_tcl_list(kv[1])
            vals = vals + [kv[i], kv[j]]
        return [self.flag] + vals

class Grp(Arg):
    """Argument grouping"""
    def init(self):

        num = self.kwds.get("num", self.kwds.get("min", None))
        if isinstance(num, int):
            self.args = [
                self.type(f"{self.name}{i+1}") for i in range(num)
            ]
            self.num = len(self.args)

        elif num=="*":
            self.args = [self.type(self.name)]
            self.num = 1

        elif isinstance(self.kwds["args"],list):
            self.args = self.kwds["args"]
            self.num = len(self.args)

        else:
            raise ValueError()



    def as_tcl_list(self, value=None):
        value = self._get_value(None,value)
        if value is None:
            if self.reqd:
                value = [None]*self.num
            else:
                return []
        return [self.flag] + [a for arg,v in zip(self.args,value) for a in arg.as_tcl_list(v)]

class One(Arg): pass

class Blk(Grp):
    def init(self):
        self.args = []
    def get_value(self, value=None):
        value = self.value if value is None else value
        if "from_prop" in self.kwds:
            return [p for v in value for p in getattr(v, self.kwds["from_prop"]) ]
        else:
            return value

    def as_tcl_list(self, value=None):
        value = self.value if value is None else value
        if value is None:
            value = [None]
        return [self.flag] + [[v.get_cmd() for v in value]]

class Tag(Int):
    def init(self):
        super().init()
        self.name = "name"
        self.field = "name"
#       self.default = int(id(self))

class Alt(Arg):
    def __init__(self, name, alts, **kwds):
        Arg.__init__(self, name, **kwds)
        self._alts = alts

    def _get_value(self, value=None):
        pass



