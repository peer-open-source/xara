from .writer import ModelWriter
from opensees.library.ast import *
from opensees.library.obj import Component, LibCmd, Ele

Frm = LibCmd("ElemData",
    args = [
        Flg("geom", enum=("linear")),
        # Vec("yornt")
    ],
    rels = {
        "ElasticBeamColumn3D": {
            ".transform.vecxz": {".yorn": lambda v: None},
            ".transform.name":   ".geom",
            ".transform.joint_offset":   ".JntOff",
            ".mass_density":     ".rho"
        }
    }
)

FIELDS = {
    "geom": "GeomData",
    "ixc": "Iz",
    "iyc": "Iy",
    "joint_offsets": "JntOff",
    "mass": "linear_mass",
    "materials": "MatData"
}

def fmt_val(arg,value):
    if isinstance(value, list):
        if isinstance(value[0],list):
            return "[" + \
                    "; ".join(", ".join(str(v) for v in l) for l in value) \
                    + "]"
        else:
            return "[" + \
                    ", ".join(str(v) for v in value) \
                    + "]"
    else:
        return str(value)

#import anabel.elements

import numpy as np
from datetime import datetime

FEDEASmap = {
    "2D beam": "Lin2dFrm",
    #anabel.elements.ElasticBeam: "LEFrame",
    #anabel.elements.ZeroLength: "zeroLength",
    "2D truss": "LETruss",
}

class FEDEAS_Writer(ModelWriter):
    def __init__(self, model, filename=None, simple=False):
        self.time_stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        self.model = model
        self.filename = filename
        self.simple = simple

        self.comment_char = "%"

        self.boun = []

    def dump_initialize(self, definitions={}):
        Domain = self.model
        # script += 'function [Model,ElemData,Loading] = {}()'.format(self.filename)
        script = ""
        #script += "\n" + "CleanStart\n"
        script += "\n".join(f"{k} = {v};" for k, v in definitions.items()) + "\n"
        return script

        # Node Definitions

    def dump_connectivity(self):
        Domain = self.model
        script = "% Node Definitions\n"
        for i, node in enumerate(Domain.nodes):
            node.name = i+1
            script += "XYZ({},:) = [{}];\n".format(
                i + 1, " ".join(f"{x:8.8}" for x in node.crd)
            )
            if any(node.kwds["boun"]):
                self.boun.append((node, node.kwds["boun"]))

        # Connectivity
        script += "\n% Connections\n"
        for i, elem in enumerate(Domain.elems):
            #ni = Domain.nodes.index(elem.nodes[0]) + 1
            #nj = Domain.nodes.index(elem.nodes[1]) + 1
            ni  = elem.nodes[0].name
            nj  = elem.nodes[1].name
            script += "CON({},:) = [{:8} {:8}];\n".format(i + 1, ni, nj)

        return script

    def dump_constraints(self):
        # Fixities
        script = ""
        Domain = self.model
        script += "%                " + " ".join(f"{i:4}" for i in Domain.dof_names.keys()) + "\n"
        for node, boun in self.boun:
            rx = " ".join(f"{i:4}" for i in boun)
            script += "BOUN({:3},:) = [{}];\n".format(node.name, rx)


        script += "\n" + "\n% Create model"
        if self.simple:
            script += "\n" + "Model = Create_SimpleModel(XYZ, CON, BOUN, ElemName);"
        else:
            script += "\n" + "Model = Create_Model(XYZ, CON, BOUN, ElemName);"
        return script

    def dump_materials(self,**kwds):
        script,sub_script = "", ""
        typs = set()
        ident = "    "
        for value in self.model.materials:
            script += "Components.{}('{}') = struct(...\n".format(
                value._cmd[0], value.name, ident
            )
            if value._cmd[0] not in typs:
                sub_script += "Components.{} = containers.Map('KeyType', 'char', 'ValueType', 'any');\n".format(
                    value._cmd[0]
                )
                typs = typs.union({value._cmd[0]})
            for sub_arg in value._args:
                sub_val = getattr(value,sub_arg.field)
                if sub_val is not None:
                    val_str = fmt_val(sub_arg, sub_val)
                    script += "{}'{}', {},...\n".format(
                         ident, FIELDS.get(sub_arg.name, sub_arg.name), val_str
                    )
            script = script[:script.rfind(",...")]
            script += ");\n"
        return sub_script + script
    

    def dump_elements(self):
        # Element types
        Domain = self.model
        from collections import defaultdict
        partials = defaultdict(list)
        for i, elem in enumerate(Domain.elems):
            elem.name = i+1
            partials[elem.prototype].append(elem.name)

        script = "% Specify element types\n"
        elem_script = ""
        script += f"\nElemData = cell({len(Domain.elems)},1);\n"
        for partial_name, elem in Domain.prototypes.items():
            if elem == "!fix":
                continue
            elem_typ = elem.type if hasattr(elem, "type") else elem.__class__.__name__
            elem_dic = {}
            for arg in elem._args:
                field = elem_dic.get(arg,arg.field)
                value = arg.get_value(elem,".m")
                import sys
                print(type(arg), id(type(arg)), id(Ref), file=sys.stderr)
                #print(type(arg), value, arg.name, isinstance(arg, Ref), file=sys.stderr)
                if value is None:
                    pass
                elif isinstance(arg, Ref):
                    script += "ElemType.{}.{} = Components.{}('{}');\n".format(
                        partial_name, FIELDS.get(arg.name, arg.name), value._cmd[0], value.name
                    )
                elif isinstance(arg, Grp) and Grp.type == Ref:
                    pass
                else:
                    script += "ElemType.{}.{} = {};\n".format(partial_name, FIELDS.get(arg.name, arg.name), value)
            elem_script += f"for el={partials[partial_name]!r}, ElemName{{el}}='{elem_typ}'; ElemData{{el}} = ElemType.{partial_name}; end;\n"
            script += "\n"

        return (script+elem_script).replace("$", "")

    def dump_loading(self):
        # Element Loads
        Domain = self.model
        script = "\n" + "\n%% Element loads"
        # Nodal Loads
        script += "\n" + "\n%% Nodal loads"
        script += "\n" + "Pf = zeros({});".format(Domain.nf)
        for node in Domain.nodes:
            p = node.p_vector()
            for i, dof in enumerate(node.dofs):
                if p[i] != 0.0:
                    script += "\n" + "Pf({}) = {};".format(dof, p[i])
        script += "\n" + "Loading.Pref = Pf;"
        return script

    def write(self):
        f = open(self.filename, "w+")
        f.write(self.string())
        f.close()

