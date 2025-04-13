import dataclasses
from io import StringIO

import numpy as np
from .registry import Registry

class ObjectSerializationError(Exception): pass

class Emitter:
    def __init__(self, strm, parent):
        self.idnt = ""
        self.strm = strm
        self.parent = parent
        self.refs = set()
        self.newline = True
        self.registry = Registry()
        # self.tags = {}

    def write(self, *args, end=" "):
        if self.newline:
            self.newline=False
            print(self.idnt, end="", file=self.strm)
        for arg in args:
            if isinstance(arg, (int,float,str, bool, np.number)):
                print(f"{arg}", end=end, file=self.strm)
            else:
                raise ValueError(f"{type(arg)}")

    def endln(self):
        self.newline = True
        print("", file=self.strm)

    def rshift(self):
        self.idnt += "\t"

    def lshift(self):
        self.idnt = self.idnt[:-1]

class ScriptBuilder:
    # TAB = object()
    # RET = object()

    def __init__(self, emitter, skip_int_refs: bool = False):
        self.emitter = emitter
        self.sent = set()
        self.python_objects = {}
        self._skip_int_refs = skip_int_refs

        self.streams = [emitter(StringIO(), self)]
        self.stream_names = {"root": self.streams[0]}

    def new_stream(self):
        self.streams.append(self.emitter(StringIO(), self))
        return self.streams[-1]

    def get_stream(self, name):
        if name not in self.stream_names:
            self.stream_names[name] = self.new_stream()
        return self.stream_names[name]

    def __repr__(self):
        return self.getScript()

    def getIndex(self):
        return "\n".join(
            (f"set {i.tclstr()} {tag}"
                for i,tag in self.registry.index().items())
        )

    def getScript(self, indexed=True)->str:
        index = self.getIndex() if indexed else ""
        return "\n\n".join((index, *[s.strm.getvalue() for s in self.streams]))

    @property
    def registry(self):
        return self.streams[0].registry

    def send(self, obj, idnt=None):
        w = self.streams[0]

        if isinstance(obj, (list,tuple)):
            #import sys
            for i in obj:
                self.send(i)
            return self

        if self.registry.registered(obj):
            return self

        if dataclasses.is_dataclass(obj):
            obj._args = dataclasses.fields(obj)
            for field in obj._args:
                value = getattr(obj, field.name)
                self.registry.register(value)
                if dataclasses.is_dataclass(value):
                    self.send(value)
                else:
                    #typ = val.__class__.__name__
                    #getattr(w, typ)(arg, value=value)
                    w.Arg(value, value=value)
            return self

        elif not hasattr(obj,"_args"):
            raise ObjectSerializationError(f"object {obj} of type {type(obj)}")

        # recurse into references
        for ref,tag_space in obj.get_refs():

            if isinstance(ref, int) and self._skip_int_refs:
                continue

            try:
                self.send(ref)

            except ObjectSerializationError:
                ident = self.registry.register(ref, tag_space=tag_space)
                self.python_objects[ident] = ref

        w.write(" ".join(obj._cmd))
        w.current_obj = obj # TODO: Clean this up

        for arg in obj._args:
            typ = arg.__class__.__name__
            value = getattr(obj, arg.field)
            if value is None and not arg.reqd:
                continue

            try:
                getattr(w, typ)(arg, value=value)
            except AttributeError:
                w.Arg(arg, value=value)

        w.endln();

        if not self.registry.registered(obj):
            self.registry.register(obj)

        return self

