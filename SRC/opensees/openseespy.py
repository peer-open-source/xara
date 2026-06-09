#===----------------------------------------------------------------------===//
#
#                                   xara
#                              https://xara.so
#
#===----------------------------------------------------------------------===//
#
# Copyright (c) 2025, OpenSees/Xara Developers
# All rights reserved.  No warranty, explicit or implicit, is provided.
#
# This source code is licensed under the BSD 2-Clause License.
# See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
#
#===----------------------------------------------------------------------===//
"""
This module implements the OpenSeesPy interface. Imports can be performed 
exactly as one would from openseespy, for example:

>>> import opensees.openseespy as ops

>>> from opensees.openseespy import node, model

>>> from opensees.openseespy import *

"""
import re
import os
import json
import pathlib 
import tempfile
import uuid
from functools import partial

from .tcl import Interpreter, _lift

# something to compare the output of model.analyze to:
successful = 0

_EXCLUDE_ECHO = {
    "nodeCoord",
    "nodeDisp",
    "nodeVel",
    "nodeReaction",
    "getNodeTags",
    "getEleTags",
    "getNDM",
    "getNDF",
    "getTime"
}


def _split_iter(source, sep=None, regex=False):
    """
    generator version of str.split()

    :param source:
        source string (unicode or bytes)

    :param sep:
        separator to split on.

    :param regex:
        if True, will treat sep as regular expression.

    :returns:
        generator yielding elements of string.
    """
    if sep is None:
        # mimic default python behavior
        source = source.strip()
        sep = "\\s+"
        if isinstance(source, bytes):
            sep = sep.encode("ascii")
        regex = True

    if regex:
        # version using re.finditer()
        if not hasattr(sep, "finditer"):
            sep = re.compile(sep)
        start = 0
        for m in sep.finditer(source):
            idx = m.start()
#           assert idx >= start
            yield source[start:idx]
            start = m.end()
        yield source[start:]

    else:
        # version using str.find(), less overhead than re.finditer()
        sepsize = len(sep)
        start = 0
        while True:
            idx = source.find(sep, start)
            if idx == -1:
                yield source[start:]
                return
            yield source[start:idx]
            start = idx + sepsize


def _obj_to_tcl(arg, name: str = None)->str:
    """
    Convert arg to a string that represents
    Tcl semantics.
    """
    import numpy as np
    if isinstance(arg, (list,np.ndarray)):
        return f"{{{' '.join(_obj_to_tcl(a) for a in arg)}}}"

    elif isinstance(arg, tuple):
        return " ".join(map(str, arg))

    # parse commands like `section Fiber {...}`
    elif isinstance(arg, dict):
        return "{\n" + "\n".join([
          f"{cmd} " + " ".join(_obj_to_tcl(a) for a in val)
              for cmd, val in arg.items()
        ]) + "}"

    elif hasattr(arg, "tag"):
        return str(arg.tag)

    else:
        return str(arg)


def _args_to_cmds(proc_name: str, *args, _final=None, **kwds):

    tcl_args = (_obj_to_tcl(i) for i in args)
    tcl_kwds = (
        (f"-{key.replace('_','-')}" if val else "") if isinstance(val, bool)
        else f"-{key} " + _obj_to_tcl(val)
            for key, val in kwds.items()
    )
    cmd = f"{proc_name} {' '.join(tcl_args)} {' '.join(tcl_kwds)}"

    if _final is not None:
        cmd += _obj_to_tcl(_final)
    return cmd


class _ModelObject:
    _tag_space: str = None 
    tag: int = None

    def _add_to_model(self, model: "OpenSeesPy", tag: int):
        raise NotImplementedError("Must be implemented by subclass")

def _is_model_object(obj):
    return hasattr(obj, "_add_to_model") and hasattr(obj, "_tag_space")

class _Surface:
    def __init__(self, nodes, cells, child, points, split):
        import shps.child
        import shps.plane
        self.nodes   = nodes
        self.cells   = cells
        self.points  = points
        self.split   = split
        self.order   = child.order
        self.outline = shps.child.IsoparametricMap(shps.plane.Q9, nodes=points)


    def walk_edge(self):
        import numpy as np
        nx, ny = self.split
        nx *= self.order
        ny *= self.order

        # setup exterior node locations in the "natural" coordinate system
        # which is on [-1,1] for x and y
        nat_exterior = [
              [ ( x, -1)  for x in np.linspace(-1, 1, nx+1)[:]],
              [ ( 1,  y)  for y in np.linspace(-1, 1, ny+1)[:]],
              [ ( x,  1)  for x in reversed(np.linspace(-1, 1, nx+1)[:])],
              [ (-1,  y)  for y in reversed(np.linspace(-1, 1, ny+1)[:])],
        ]

        def find_node(coord):
            for tag, xyz in self.nodes.items():
                if np.linalg.norm(np.array(xyz) - np.array(coord)) <= 1e-12:
                    return tag


        # get the number of nodes per element
        nen = self.order + 1
        for i,edge in enumerate(nat_exterior):
            for j in range(self.split[i%2]):
                yield tuple(find_node(self.outline.coord(xn)) for xn in edge[j*(nen-1):j*(nen-1)+nen])



class _TagRegistry:
    def __init__(self):
        self._tags =  {
            "uniaxialMaterial": set(),
            "nDMaterial": set(),
            "section": set(),
            "pattern": set(),
            "timeSeries": set(),
            # "frame_section": set(),
            # "shell_section": set(),
            # "plane_section": set(),
            # Dont do elements because they can be added in block2D and block3D without going through the element() method
        }

        self._objects = {
            k: {} for k in self._tags.keys()
        }

    def add(self, type_name: str, tag: int, obj=None):
        self._tags[type_name].add(tag)
    
    def set(self, type_name: str, obj, tag: int):
        # self._objects[type_name][id(obj)] = tag
        pass

    def get(self, type_name: str)->int:
        tag = 1
        while tag in self._tags[type_name]:
            tag += 1
        self._tags[type_name].add(tag)
        return tag


class OpenSeesPy:
    """
    This class is meant to be instantiated as a global singleton
    that is private to this Python module.

    It encapsulates an instance of Interpreter which implements an
    OpenSees state.
    """
    def __init__(self, *args,
                 save=False,
                 echo_file=None,
                 **kwds):

        error_file = pathlib.Path(tempfile.gettempdir())/f"{uuid.uuid4()}"
        self._interp  = Interpreter(*args,
                                    error_file=error_file, 
                                    **kwds)
        self._partial = partial
        self._save    = save


        # Setup an echo file to capture all Tcl commands
        if echo_file is None and "XARA_ECHO_FILE" in os.environ:
            mode = os.environ.get("XARA_ECHO_MODE", "w+")
            echo_file = open(os.environ["XARA_ECHO_FILE"], mode)

        self._echo = echo_file

        self._mesh = {"line": {}, "quad": {}}

        self._tags = _TagRegistry()

        # Enable OpenSeesPy command behaviors
        self.eval("pragma openseespy")
        self.eval("set tcl_precision 17")



    def _invoke_proc(self, proc_name: str, *args, _final=None, _return_string=False, **kwds)->object:
        """
        Invoke the Interpreter's eval method, calling
        a procedure named `proc_name` with arguments
        from args and kwds, after converting Python semantics
        to Tcl semantics (via _obj_to_tcl).

        For example, key-word arguments contained in the `kwds`
        dict are converted to a sequence of "-key" and "value"
        strings.

        """
        if False:
            do_call = True 
            for arg in args:
                if isinstance(arg, (list, dict, tuple)):
                    do_call = False
                    break
            
            if do_call and len(kwds) == 0:
                return self._call(proc_name, *args, **kwds)

        comment = ""
        if "comment" in kwds:
            comment = kwds.pop("comment")
            if isinstance(comment, str):
                comment = f"; # {comment}"
            else:
                raise TypeError(f"Invalid type for comment: {type(comment)}")

        cmd = _args_to_cmds(proc_name, *args, _final=_final, **kwds)
        cmd += comment

        #
        #
        try:
            ret = self.eval(cmd)
        except Exception as e:
            from .errors import XaraError
            raise XaraError(f"{e}").with_traceback(None) from None
        #
        #

        if ret is None or ret == "":
            return None

        if _return_string:
            return ret

        parts = ret.split()
        # Use json parse to cast return values from string. 
        # This is faster than the standard ast module.
        if len(parts) > 1:
            try:    return list(map(json.loads, parts))
            except: return ret

        elif proc_name in {"eigen"}:
            # "eigen" should always return a list
            return [float(ret)]

        else:
            try:    return json.loads(ret)
            except: return ret


    def echo(self, *args):
        print(*args)
        for arg in args:
            self.eval(f'puts "{arg}"')


    def eval(self, cmd: str, echo: bool=None) -> str:
        "Evaluate a Tcl command"
        try:
            if self._echo is not None and cmd.split(" ")[0] not in _EXCLUDE_ECHO:
                print(cmd, file=self._echo)
        except:
            pass
        return self._interp.eval(cmd)

    def _call(self, proc_name: str, *args, **kwds):
        """
        EXPERIMENTAL (2025-07-04)
        """
        if self._echo is not None and proc_name not in _EXCLUDE_ECHO:
            print(_args_to_cmds(proc_name, *args, **kwds),
                  file=self._echo)

        return self._interp._tcl.call(proc_name, *args, **kwds)

    def block3D(self, *args, **kwds):
        if isinstance(args[6], list) or isinstance(args[7], dict):
            return self._invoke_proc("block3D", *args, **kwds)

        # We have to imitate the OpenSeesPy parser, which
        # *requires* hard-coding the number of element args
        # expected by each element type. This is terribly
        # unstable and limited and should only be used when 
        # backwards compatibility with the original OpenSeesPy 
        # is absolutely necessary.
        elem_name = args[4]
        elem_argc = 7
        elem_args = args[6]

        nl  = '\n'
        ndm = self._call("getNDM")

        # loop over remaining args to form node coords
        node_args = f"""{{
            {nl.join(" ".join(map(str,args[elem_argc+i*(ndm+1):elem_argc+(i+1)*(ndm+1)])) for i in range(int(len(args[elem_argc:])/(ndm+1))))}
        }}"""

        return self._invoke_proc("block3D", *args[:6], elem_args, node_args)


    def block2D(self, *args, **kwds):
        if isinstance(args[5], list):
            return self._invoke_proc("block2D", *args, **kwds)
        from ._model.block import block2D
        return block2D(self, *args, **kwds)

    def getNodeTags(self)->list[int]:
        tags = self._call("getNodeTags")
        if tags is None:
            return []
        elif isinstance(tags, int):
            return [tags]
        else:
            return tags
        
    def getEleTags(self)->list[int]:
        tags = self._call("getEleTags")
        if tags is None:
            return []
        elif isinstance(tags, int):
            return [tags]
        else:
            return tags

    def timeSeries(self, *args, **kwds):
        """
        ['Path', 1, '-values', 0.0, 5.0, 8.0, 7.0, 5.0, 3.0, 2.0, 1.0, 0.0, '-time', 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
        ['Path', 1, '-values', [0.0, 5.0, 8.0, 7.0, 5.0, 3.0, 2.0, 1.0, 0.0], '-time', [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]]
        """

        args = list(args)
        if "-values" in args:
            iv = args.index("-values")
            # Count the number of floating-point arguments
            for nv, value in enumerate(args[iv+1:]):
                if not isinstance(value, float):
                    nv += 1
                    break
            else:
                # if we didnt break out of the for loop
                nv += 2

            values = args[iv+1:iv+nv]
            args = [a for a in args[:iv+1]] + [values] + [a for a in args[iv+nv:]]

        if "-time" in args:
            it = args.index("-time")
            for nt, value in enumerate(args[it+1:]):
                if not isinstance(value, float):
                    nt += 1
                    break
            else:
                # if we didnt break out of the for loop
                nt += 2

            time = args[it+1:it+nt]
            args = [a for a in args[:it+1]] + [time] + [a for a in args[it+nt:]]

        return self._invoke_proc("timeSeries", *args, **kwds)


    def pattern(self, *args, load=None, **kwds):

        self._current_pattern = args[1]

        if load is None and "loads" in kwds:
            load = kwds.pop("loads")

        if load is not None:
            loads = [
                    ("load", k, *v, ";\n") for k,v in load.items()
            ]
            return self._invoke_proc("pattern", *args, **kwds, _final=loads)
        else:
            return self._invoke_proc("pattern", *args, **kwds)

    def load(self, *args, pattern=None, load=None, **kwds):
        if pattern is None:
            pattern = self._current_pattern
        return self._invoke_proc("nodalLoad", *args, "-pattern", pattern, **kwds)


    def mesh(self, type, tag: int, *args, **kwds):
        if type == "line":
            return self._mesh_line(tag, 2, args[1:3], *args[3:7], args[7:])
        raise NotImplementedError(f"Mesh type '{type}' not implemented")

    def _mesh_line(self, tag, numnodes, ndtags, id, ndf:int, meshsize, eleType='', eleArgs=()):
        import numpy as np
        from itertools import count

        ndI, ndJ = ndtags
        add_node    = partial(self._invoke_proc, "node")
        add_element = partial(self._invoke_proc, "element")

        xi = np.array(self._invoke_proc("nodeCoord", ndI))
        xj = np.array(self._invoke_proc("nodeCoord", ndJ))

        L  = np.linalg.norm(xj - xi)
        nn = int(L//meshsize) + 1

        nodes = [None for _ in range(nn)]
        nodes[0]    = ndI
        nodes[nn-1] = ndJ

        node_tags = set(self._invoke_proc("getNodeTags"))
        new_node  = filter(lambda i: i not in node_tags, count(1))
        elem_tags = set(self._invoke_proc("getEleTags") or [])
        new_elem  = filter(lambda i: i not in elem_tags, count(1))

        for i,x in enumerate(np.linspace(xi, xj, nn, endpoint=True)[1:]):

            node_tag = next(new_node)
            add_node(node_tag, *x)

            nodes[i+1] = node_tag

            elem_tag = next(new_elem)

            if i < nn and eleType != '' and eleArgs:
                add_element(eleType,elem_tag,nodes[i],nodes[i+1],*eleArgs)

        self._mesh["line"][tag] = nodes


    def add_object(self, obj: _ModelObject, tag: int=None):
        if tag is None:
            if hasattr(obj, "tag") and obj.tag is not None:
                tag = obj.tag
            else:
                tag = self._tags.get(obj._tag_space)

        self._tags.set(obj._tag_space, obj, tag)
        obj._add_to_model(self, tag)


    def material(self, type, tag: int, *args, **kwds):
        if not isinstance(type, str):
            return type._add_to_model(self, tag)
        else:
            return self._invoke_proc("material", type, tag, *args, **kwds)

    def section(self, type: str, sec_tag: int, *args, **kwds):
        self._current_section = sec_tag

        # TODO: error handling

        if not isinstance(type, str):
            return type._add_to_model(self, sec_tag)
        elif len(args) > 0 and not isinstance(args[0], (int, float, str)):
            # model.section(name, tag, shape, materials, args)
            return self._add_section_shape(type, sec_tag, *args, **kwds)

        # Undocumented feature
        if "shape" in kwds:
            from opensees.section import from_shape
            ndm = int(self.eval("getNDM"))
            # kwds["shape"] looks like ("W14X90", matTag, (20,4), units?)
            shape = from_shape(type, *kwds.pop("shape"), ndm=ndm)
        else:
            shape = None
        
        try:
            ret = self._invoke_proc("section", type, sec_tag, *args, **kwds)
        except Exception as e:
            from .errors import XaraError
            raise XaraError(f"{e}").with_traceback(None) from None

        if shape is not None:
            for fiber in shape.fibers:
                self._invoke_proc("fiber", *fiber.coord, fiber.area, fiber.material, section=sec_tag)

        return ret
    
    def _add_section_shape(self, name, tag, materials, shape=None,  **kwds):
        from xara import Section, ShellSection
        if shape is None:
            shape = materials 
            materials = None

        if name in {"ShellFiber"}:
            section = ShellSection(name, materials, shape, **kwds)
        else:
            section = Section(name, shape, materials, **kwds)

        section._add_to_model(self, tag)
        return section

    def patch(self, *args, **kwds):
        if "section" not in kwds:
            kwds["section"] = self._current_section
        return self._invoke_proc("patch", *args, **kwds)

    def layer(self, *args, **kwds):
        if "section" not in kwds:
            kwds["section"] = self._current_section
        return self._invoke_proc("layer", *args, **kwds)

    def fiber(self, *args, **kwds):
        if "section" not in kwds:
            kwds["section"] = self._current_section
        return self._invoke_proc("fiber", *args, **kwds)

from collections import defaultdict

class Model:
    def __init__(self, *args, echo_file=None, **kwds):
        self._openseespy = OpenSeesPy(echo_file=echo_file)
        if len(args) > 0 or len(kwds) > 0:
            self._openseespy._invoke_proc("model", *args, **kwds)

        self._parameters = {
        }

        # Aug 2025, for xara._analysis
        self._patterns = {}

        # Apr 2026
        self._objects = defaultdict(dict)
    
    @property
    def state(self):
        return StateView(self)


    def eval(self, *args, **kwds):
        return self._openseespy.eval(*args, **kwds)

    def _call(self, proc_name: str, *args, **kwds):
        """
        EXPERIMENTAL (2025-07-04)
        """
        return self._openseespy._call(proc_name, *args, **kwds)


    def export(self, *args, **kwds):
        return self._openseespy._interp.export(*args, **kwds)

    def lift(self, type_name: str, tag: int):
        return _lift(self._openseespy._interp._tcl.interpaddr(), type_name, tag)


    # def invoke(self, *args, **kwds):
    #     if len(args) == 2:
    #         from ._invoke import _Handle
    #         return _Handle(self._openseespy, *args, **kwds)
    #     else:
    #         return self._openseespy._invoke_proc(*args, **kwds)

    
    def asdict(self):
        """April 2024"""
        return self._openseespy._interp.serialize()


    def node(self, tag: int, *args, **kwds):
        return self._openseespy._invoke_proc("node", tag, *args, **kwds)
    
    def add_object(self, obj, tag: int=None):
        r = self._openseespy.add_object(obj, tag)
        self._objects[obj._tag_space][obj.tag] = obj
        return r

    def material(self, type_or_object, *args, **kwds):
        if _is_model_object(type_or_object):
            # return self._openseespy.add_object(type_or_object)
            return self.add_object(type_or_object)
        else:
            return self._openseespy.material(type_or_object, *args, **kwds)
    
    def section(self, type_or_object, *args, **kwds):
        if _is_model_object(type_or_object):
            # return self._openseespy.add_object(type_or_object)
            return self.add_object(type_or_object)
        else:
            return self._openseespy.section(type_or_object, *args, **kwds)


    def element(self, type, tag, *args, **kwds):
        if tag is None:
            tag = 1
            ele_tags = self.getEleTags()
            if ele_tags is None:
                ele_tags = []
            elif isinstance(ele_tags, int):
                ele_tags = [ele_tags]

            for existing_tag in ele_tags:
                if tag <= existing_tag:
                    tag = existing_tag + 1
        try:
            self._openseespy._invoke_proc("element", type, tag, *args, **kwds)
        except Exception as e:
            raise e.with_traceback(None) from None
        return tag

    def symbols(self, **kwds):
        symbols = []
        for k,v in kwds.items():
            self.eval(f"set {k} {v}")
            symbols.extend((f"-{k}", f"${k}"))
        return symbols


    def surface(self, split, 
                element: str=None, 
                args=None, 
                points=None, 
                name=None, 
                kwds=None, 
                order=None, 
                shape=None):
        """
        Create a surface mesh of elements in the current model.

        Parameters
        ----------
        :param split: tuple of integers
            The number of elements in the x and y directions.
        :param element: str
            The name of the element type to use.
        :param args: tuple
            The arguments to pass to the element constructor.
        :param points: list of tuples
            The coordinates of the points in the mesh.
        :param name: str
            The name of the mesh.
        :param kwds: dict
            The keyword arguments to pass to the element constructor.
        :param order: int
            The order of the elements to use.
        :param shape: str
            The shape of the elements to use. Can be "Q" for quadrilateral or "T" for triangular.
        :return: Surface
        """
        # anchor
        # normal
        import shps.plane, shps.block

        add_node    = partial(self._openseespy._invoke_proc, "node")
        add_element = partial(self._openseespy._invoke_proc, "element")

        cell_type = None

        if shape is None:
            if element in {"ShellMITC4", "Q4", "Q8", "Q9"}:
                cell_type = {
                    "ShellMITC4": shps.plane.Q4,
                    "Q4": shps.plane.Q4,
                    "Q8": shps.plane.Q8,
                    "Q9": shps.plane.Q9,
                }[element]
            else:
                shape = "Q"

        if cell_type is None:
            if order == 1 and shape == "Q":
                cell_type = shps.plane.Q4
            elif order == 2 and shape == "Q":
                cell_type = shps.plane.Q9
            elif order == 2 and shape == "T":
                cell_type = shps.plane.T6


        if isinstance(element, str) and cell_type is None:
            # element is an element name
            cell_type = {
                    "ShellMITC4": shps.plane.Q4,
            }.get(element, shps.plane.Q4)

        elif element is None and name is None:
            cell_type = shps.plane.Q4
            element = None

        elif isinstance(name, str):
            cell_type = element
            element = name


        m_elems = self._openseespy._invoke_proc("getEleTags")
        if isinstance(m_elems, int):
            m_elems = {m_elems}

        elif m_elems is not None:
            m_elems = {tag for tag in m_elems}

        m_nodes = self._openseespy._invoke_proc("getNodeTags")
        if isinstance(m_nodes, int):
            m_nodes = {m_nodes}
        if m_nodes is not None:
            m_nodes = {
                int(tag): self._openseespy._invoke_proc("nodeCoord", f"{tag}")
                for tag in m_nodes
            }

        if m_nodes is not None and len(m_nodes) > 0:
            join = dict(nodes=m_nodes, cells=m_elems)
        else:
            join = None

        if kwds is None:
            kwds = {}

        nodes, elems = shps.block.block(split, cell_type, points=points,
                                        append=False, join=join, **kwds)


#       anchor_point, anchor_coord = next(iter(anchor.items()))
#       if isinstance(anchor_coord, int):
#           anchor_coord = self._openseespy._str_call("nodeCoord", f"{anchor_coord}")

#       anchor_point = np.array([*nodes[anchor_point], 0.0])
        for tag, coord in nodes.items():
            add_node(tag, *coord)

        if isinstance(args, dict):
            ekwds = args
            args = ()
        else:
            ekwds = {}

        if element is not None:
            for tag, elem_nodes in elems.items():
                add_element(element, tag, list(map(int,elem_nodes)), *args, **ekwds)

        return _Surface(nodes=nodes,
                        cells=elems,
                        child=cell_type,
                        points=points,
                        split=split)
    #
    # Loading
    #

    def pattern(self, *args, **kwds):
        if _is_model_object(args[0]):
            return self.add_object(args[0])
        else:
            return self._openseespy.pattern(*args, **kwds)

    # 
    # Analysis 
    #
    def getIterationCount(self):
        return self._openseespy._invoke_proc("numIter")
    
    def testNorms(self, *args):
        return self._call("testNorms", *args)

    def getResidual(self):
        import numpy as np
        residual_string = self._openseespy._invoke_proc("printB", "-ret", _return_string=True)
        n = sum(1 for _ in _split_iter(residual_string))
        return np.fromiter(map(float, _split_iter(residual_string)), count=n, dtype=float)


    def getTangent(self, **kwds):
        import numpy as np

        tangent_string = self._openseespy._invoke_proc("printA", "-ret", _return_string=True, **kwds)

        nn = sum(1 for _ in _split_iter(tangent_string))

        A  = np.fromiter(

                map(float, _split_iter(
                    tangent_string
                )),

                count=nn,
                dtype=float
        )

        # Assigning to .shape as opposed to calling .reshape()
        # should enforce no copying
        A.shape = tuple([int(np.sqrt(len(A)))]*2)

        # For large systems, avoid clogging memory
        if nn > 100:
            import gc
            del tangent_string
            gc.collect()
        return A

    def __getattr__(self, name: str):
        if name in _OVERWRITTEN:
            return getattr(self._openseespy, name)
        elif name in __all__ or name in {"print"}:
            return self._openseespy._partial(self._openseespy._invoke_proc, name)
        else:
            raise AttributeError(f"class Model has no attribute '{name}'")




# A list of symbol names that are importable
# from this module. All of these are dynamically
# resolved by the function __getattr__ below.
__all__ = [
# 
    "tcl",
    "OpenSeesError",
    "invoke",

# OpenSeesPy attributes

    "uniaxialMaterial",
    "testUniaxialMaterial",
    "setStrain",
    "getStrain",
    "getStress",
    "getTangent",
    "getDampTangent",
    "wipe",
    "model",
    "node",
    "element",
    "section",
    "fiber", "patch", "layer",
    "geomTransf",
    "transform",
    "beamIntegration",
    "nDMaterial",
    "material",
    "frictionModel",
    "limitCurve",
    # Constraints
    "fix",
    "sp",
    "fixX",
    "fixY",
    "fixZ",
    "mass",
    "constrain",
    "equalDOF",
    "equalDOF_Mixed",
    "rigidLink",
    "rigidDiaphragm",
    # Damping
    "rayleigh",
    "modalDamping",
    "modalDampingQ",
    "setElementRayleighDampingFactors",
    # Meshing
    "block2D",
    "block3D",
    "mesh",
    "remesh",

    # Loading
    "timeSeries",
    "pattern",
    "load",
    "eleLoad",
    "imposedMotion",
    "imposedSupportMotion",
    "groundMotion",
    # Analysis
    "loadConst",
    "system",
    "numberer",
    "constraints",
    "integrator",
    "algorithm",
    "analysis",
    "analyze",
    "test",
    "modalProperties",
    "responseSpectrumAnalysis",
    "eigen",
    #
    "reactions",
    "setTime",
    "remove",
    "setCreep",
    "reset",
    "initialize",
    "getLoadFactor",
    "build",
    "printModel",
    "printA",
    "printB",
    "printGID",
    "testNorm",
    "testNorms",
    "testIter",
    "recorder",
    "database",
    "save",
    "restore",

    "getTime",
    # node
    "nodeUnbalance",
    "nodeDisp",
    "nodeRotation",
    "nodeEigenvector",
    "nodeVel",
    "nodeAccel",
    "nodeReaction",
    "nodeResponse",
    "nodeCoord",
    # element
    "eleResponse",
    "eleForce",
    "eleDynamicalForce",

    "updateElementDomain",

    "getNDM",
    "getNDF",
    "eleNodes",
    "eleType",
    "nodeDOFs",
    "nodeMass",
    "nodePressure",
    "setNodePressure",
    "nodeBounds",
    "getPatterns",
    "getFixedNodes",
    "getFixedDOFs",
    "getConstrainedNodes",
    "getConstrainedDOFs",
    "getRetainedNodes",
    "getRetainedDOFs",
    "getNodeTags",
    # Setters
    "setNodeCoord",
    "setNodeVel",
    "setNodeAccel",
    "setNodeDisp",

    "start",
    "stop",
    "region",
    "setPrecision",
    "searchPeerNGA",
    "domainChange",
    "record",
    "metaData",
    "defaultUnits",
    "stripXML",
    "convertBinaryToText",
    "convertTextToBinary",
    "getEleTags",
    "getCrdTransfTags",
    "getParamTags",
    "getParamValue",
    "sectionForce",
    "sectionDeformation",
    "sectionStiffness",
    "sectionFlexibility",
    "sectionLocation",
    "sectionWeight",
    "sectionTag",
    "sectionDisplacement",
    "cbdiDisplacement",
    "basicDeformation",
    "basicForce",
    "basicStiffness",
    "InitialStateAnalysis",
    #
    "totalCPU",
    "solveCPU",
    "accelCPU",
    "numFact",
    "numIter",
    "wipeAnalysis",
    "systemSize",
    "version",
    "setMaxOpenFiles",
    "ShallowFoundationGen",
    "setElementRayleighFactors",
    # Parameters
    "parameter",
    "addToParameter",
    "updateParameter",
    "setParameter",
    # Sensitivity
    "computeGradients",
    "sensitivityAlgorithm",
    "sensNodeDisp",
    "sensNodeVel",
    "sensNodeAccel",
    "sensLambda",
    "sensSectionForce",
    "sensNodePressure",
    # Query
    "getNumElements",
    "getEleClassTags",
    "getEleLoadClassTags",
    "getEleLoadTags",
    "getEleLoadData",
    "getNodeLoadTags",
    "getNodeLoadData",
    "IGA",
    "NDTest",
    # Parallel
    "getPID",
    "getNP",
    "barrier",
    "send",
    "recv",
    "Bcast",
    # Reliability
    "randomVariable",
    "getRVTags",
    "getRVParamTag",
    "getRVValue",
    "getMean",
    "getStdv",
    "getPDF",
    "getCDF",
    "getInverseCDF",
    "correlate",
    "performanceFunction",
    "gradPerformanceFunction",
    "transformUtoX",
    "wipeReliability",
    "updateMaterialStage",
    "sdfResponse",
    "probabilityTransformation",
    "startPoint",
    "randomNumberGenerator",
    "reliabilityConvergenceCheck",
    "searchDirection",
    "meritFunctionCheck",
    "stepSizeRule",
    "rootFinding",
    "functionEvaluator",
    "gradientEvaluator",
    "getNumThreads",
    "setNumThreads",
    "logFile",
    "setStartNodeTag",
    "hystereticBackbone",
    "stiffnessDegradation",
    "strengthDegradation",
    "strengthControl",
    "unloadingRule",
    "partition",
    "pressureConstraint",
    "domainCommitTag",
#   "runFOSMAnalysis",
    "findDesignPoint",
    "runFORMAnalysis",
    "getLSFTags",
    "runImportanceSamplingAnalysis",
]


# Commands that are pre-processed in Python
# before forwarding to the Tcl interpreter
_OVERWRITTEN = {
    "timeSeries",
    "pattern", "load",
    "eval",
    "material",
    "section", "patch", "layer", "fiber",
    "block2D",
    "block3D",
    "mesh",
    "getNodeTags",
    "getEleTags",
}



class State:
    class NodalVector:
        def __init__(self, dofs, values=None):
            self._dofs   = dofs
            self._values = values

        def __call__(self, node=None, dof=None):
            if node is None:
                return self._values
    
            elif dof is None:
                ndofs = self._dofs.get(node)
                if ndofs is None:
                    return None
                return self._values[node]

            else:
                ndofs = self._dofs.get(node)
                if ndofs is None or dof > len(ndofs):
                    return None
                return self._values[node][dof-1]


    def __init__(self, model,
                 u=None,
                 v=None,
                 a=None,
                 reactions=None,
                 time=None):
        self._nodes = {
            node: model._call("nodeDOFs", node)
                for node in model._call("getNodeTags") or []
        }
        self._u         = u
        self._v         = v
        self._a         = a
        self._reactions = reactions
        self._time      = time

    @property
    def u(self):
        return self.NodalVector(self._nodes, self._u)

    def v(self):
        return self.NodalVector(self._nodes, self._v)

    def a(self):
        return self.NodalVector(self._nodes, self._a)

    def reactions(self):
        return self.NodalVector(self._nodes, self._reactions)

    def time(self):
        return self._time



class StateView:
    class _NodalVector:
        def __init__(self, model, response: str):
            self._model = model
            self._response = response

        def __call__(self, node=None, dof=None, element=None):
            if dof is not None:
                return self._model._call(f"node{self._response.capitalize()}", node, dof)
            else:
                return self._model._call(f"node{self._response.capitalize()}", node)
        
        def __getitem__(self, index):
            if isinstance(index, tuple):
                node, dof = index
                return self._model._call(f"node{self._response.capitalize()}", node, dof)
            else:
                node = index
                return self._model._call(f"node{self._response.capitalize()}", node)


    def __init__(self, model):
        self._model = model

    def store(self, fields=None, file=None):
        if fields is None:
            fields = ["u", "v", "a", "reactions"]

        data = {
            key: getattr(self, key).values for key in fields
        }

        state = State(self._model, time=self.time, **data)
        if file is None:
            return state
        else:
            # TODO
            pass

    @property
    def time(self):
        return self._model._call("getTime")

    @property
    def u(self):
        return self._NodalVector(self._model, "disp")
    
    @property
    def v(self):
        return self._NodalVector(self._model, "vel")
    
    @property
    def a(self):
        return self._NodalVector(self._model, "accel")

    # @property
    # def tangent(self):
    #     return self._model.getTangent()
    
    # @property 
    # def residual(self):
    #     return self._model.getResidual()
    
    @property
    def reactions(self):
        return self._NodalVector(self._model, "reaction")
    
    def element(self, tag: int):
        class ElementStateView:
            def __init__(self, model, tag: int):
                self._model = model
                self._tag   = tag

            def force(self, dof: int=None):
                if dof is None:
                    return self._model._call("eleForce", self._tag)
                else:
                    return self._model._call("eleForce", self._tag, dof)

            def deformation(self, dof: int=None):
                if dof is None:
                    return self._model._call("eleDeformation", self._tag)
                else:
                    return self._model._call("eleDeformation", self._tag, dof)
            
            def section(self, tag: int):
                class SectionStateView:
                    def __init__(self, model, ele_tag: int, sec_tag: int):
                        self._model   = model
                        self._ele_tag = ele_tag
                        self._sec_tag = sec_tag

                    def tangent(self, expand=False):
                        import numpy as np
                        if not expand:
                            K =  self._model._openseespy._invoke_proc(
                                "eleResponse", 
                                self._ele_tag,
                                "section", self._sec_tag,
                                "tangent"
                            )
                        else:
                            K = self._model._openseespy._invoke_proc("sectionStiffness", 
                                                                    self._ele_tag, self._sec_tag)
                        Ka = np.array(K)
                        Ka.shape = tuple([int(np.sqrt(len(Ka)))]*2)
                        return Ka

                    def stress(self):
                        return self._model._openseespy._invoke_proc(
                            "eleResponse", 
                            self._ele_tag,
                            "section", self._sec_tag,
                            "resultant"
                        )
                    def strain(self):
                        return self._model._openseespy._invoke_proc(
                            "eleResponse", 
                            self._ele_tag,
                            "section", self._sec_tag,
                            "deformation"
                        )

                return SectionStateView(self._model, self._tag, tag)
        return ElementStateView(self._model, tag)


    def print(self):
        for node in self._model.getNodeTags():
            u = self.u(node)
            print(f"  {node}: " + "  ".join(f"{val:14.6g}" for val in u))


# The global singleton, for backwards compatibility
try:
    _openseespy = OpenSeesPy()
except:
    _openseespy = None


def __getattr__(name: str):
    # For reference:
    #   https://peps.python.org/pep-0562/#id4
    if name in _OVERWRITTEN:
        return getattr(_openseespy, name)
    else:
        return _openseespy._partial(_openseespy._invoke_proc, name)

