from .writer import ModelWriter
from .emitter import Emitter, ScriptBuilder
from opensees.library.ast import Arg


class TclWriter(Emitter):
    def Arg(this, self, value=None)->list:
        if self.__class__.__name__ in dir(this):
            getattr(this, self.__class__.__name__)(self, value=value)
            return

        value = self._get_value(None, value)

        if isinstance(value, (type(None), )):
            if not self.reqd:
                return
            else:
                value = f"${self.name}"
        elif isinstance(value, (Arg,)):
            if not self.reqd:
                return
            else:
                value = f"${value.name}"

        this.write(self.flag, str(value))


    def Lst(this, self, value=None):
        this.write(self.flag)
        this.write("{")
        for v in map(self.type, value):
            this.write(v)
        this.write("}")

    def Tag(this, self, value=None):
        value = self.value if value is None else value

        if isinstance(value, int):
            force_tag = value
        else:
            force_tag = None

        ident = this.registry.register(this.current_obj, name=value, force_tag=force_tag)

        this.write("$"+ident.tclstr())

    def Flg(this, self, value=None):
        value = self.value if value is None else value
        if value: this.write(self.flag)

    def Grp(this, self, value=None):
        val = self._get_value(None,value)
        if value is None:
            if self.reqd:
                value = [None]*self.num
            else:
                return []
        if "reverse" in self.kwds and self.kwds["reverse"]:
            value = reversed(value)

        this.write(self.flag)
        if not hasattr(value,"__len__") or len(self.args) == len(value):
            for a,v in zip(self.args, value):
                this.Arg(a, value=v)
        else:
            # assert "min" in self.kwds
            assert len(self.args) == 1
            for v in value:
                this.Arg(self.args[0], value=v)


        # this.write(self.flag, " ".join(map(str,(
        #     a for arg,v in zip(self.args,value) for a in arg.as_tcl_list(v)
        # ))))
        # [this.parent.send(v) for v in value]

    def Ref(this, self, value=None):
        val = self._get_value(None, value=value)
        #value = value.name if val is None else val
        #if val is None:
        #if not isinstance(val,int):
        if not isinstance(val,int) and not isinstance(value, int):
            try:
                value = "$"+this.registry.ident(value).tclstr()
            except KeyError:
                if self.reqd:
                    raise KeyError(f"Unable to resolve reference to {value} when writing {this.current_obj}")
                else:
                    return

        else:
            value = val or value # the one thats an int, not (presumably) None
        return this.write(self.flag,value)

    def Blk(this, self, value=None):
        value = self.get_value(value=value)
        if value is None:
            value = [None]
        this.write(self.flag, "{")
        this.endln()
        this.rshift()
        for v in value:
            this.parent.send(v)
        this.lshift()
        this.write("}")

    def Map(this, self, value=None):
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

        return this.write(self.flag, *vals)

TclScriptBuilder = lambda **kwds: ScriptBuilder(TclWriter, **kwds)

class OpenSeesWriter(ModelWriter):
    def __init__(self, model=None):
        if hasattr(model,"apply"):
            model = model.apply(model.prototypes)
        self.model = model
        self.comment_char = "#"

    def dump_initialize(self, definitions={}):
        cmds = "# Parameters\n" + "\n".join(f"set {k} {v};" for k,v in definitions.items()) + "\n\n"
        ndm, ndf = self.model.ndm, self.model.ndf
        dof_keys = "dx dy dz rx ry rz"
        dofs = " ".join( str(i) for i in range(1, ndf + 1))
        cmds += f"# Create ModelBuilder (with {ndm} dimensions and {ndf} DOF/node)\n"
        cmds += f"model BasicBuilder -ndm {ndm} -ndf {ndf}\n"
        cmds += f"lassign {{{dofs}}} {dof_keys}"
        return cmds

    def dump_elements(self, *elems, definitions={}):
        transforms = set()

        cmds = "\n".join(f"set {k} {v};" for k,v in definitions.items()) + "\n"
        for i, el in enumerate(elems):
            #el.init()
            if el.name is None:
                el.name = i+1

            try:
                cmds += "\nelement " + " ".join(str(p) for p in el.ops_elem.serialize(tag=i))
            except:
                cmds += "\n" + " ".join(el.get_cmd_str()) + "\n"

            if hasattr(el, "_transform") and el._transform:
                transforms.update({el._transform})

        return "".join(" ".join(t.get_cmd_str()) + "\n" for t in transforms) + cmds

    @classmethod
    def dump_sections(self, *sections, definitions={}):
        cmds = "\n".join(f"set {k} {v};" for k,v in definitions.items()) + "\n"
        for sect in sections:
            cmds += " ".join(sect.get_cmd_str()) + "\n"
        return cmds

    def dump_materials(self, *materials, definitions={}):
        builder = TclScriptBuilder()
        writer = builder.streams[0]
        cmds = "\n".join(f"set {k} {v};" for k,v in definitions.items()) + "\n"
        writer.write(
            "\n".join(f"set {k} {v};" for k,v in definitions.items()) + "\n"
        )
        all_materials = self.model.materials
        for i,mat in enumerate(all_materials):
            if mat.name is None:
                mat.name = i+1
            cmds += " ".join(mat.get_cmd_str()) + "\n"
            writer.parent.send(mat)
        return builder.getScript()

    def dump_constraints(self, definitions={}):
        cmds = "\n".join(f"set {k} {v};" for k,v in definitions.items()) + "\n"
        for n in self.model.nodes:
            if 'boun' in n.kwds and any(n.kwds["boun"]):
                cmds += f"fix {n.name} {' '.join(str(i) for i in n.kwds['boun'])}\n"

        return cmds

    def dump_connectivity(self, definitions={}):
        cmds = "\n".join(f"set {k} {v};" for k,v in definitions.items()) + "\n"
        for i,nd in enumerate(self.model.nodes):
            nd.name = i + 1
            cmds +=  " ".join(nd.get_cmd_str()) + "\n"
        cmds += self.dump_elements(*self.model.elems)
        return cmds


