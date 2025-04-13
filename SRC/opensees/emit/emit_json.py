import sys
import dataclasses
import numpy as np
from .emitter import Emitter, ScriptBuilder, ObjectSerializationError

SEQUENTIAL, ASSOCIATIVE, TERMINAL = map(lambda x: 1<<x, range(3))

def _is_primitive(obj):
    return isinstance(obj, (int, float, bool, type(None), str, np.ndarray))
    return not _is_collection(obj)

def _is_collection(obj, mapping:int = SEQUENTIAL|ASSOCIATIVE|TERMINAL):
    test =  (
        (dataclasses.is_dataclass(obj) or isinstance(obj, dict)) and mapping&ASSOCIATIVE
    ) or (
        isinstance(obj, (list, set, tuple)) and mapping&SEQUENTIAL
    )


    if test and (mapping == SEQUENTIAL):
        assert isinstance(obj, (tuple, set, list, np.ndarray))

    elif test and (mapping == ASSOCIATIVE):
        assert isinstance(obj, (dict)) or dataclasses.is_dataclass(obj)

    return test

class SerializationStream:
    pass


class Partial: pass

class _Builder(ScriptBuilder):
    def serialize_partial(self, obj):
        if _is_primitive(obj):
            w.Arg(obj)
        elif _is_collection(obj, SEQUENTIAL):
            return list()
        elif _is_collection(obj, ASSOCIATIVE):
            return dict()

    def serialize_references(self, obj):
        if isinstance(obj, dict):
            w.write("{")
            w.rshift()
            for k,v in obj.items():
                w.write(f'"{k}": ')
                self.send(v)
            w.lshift()
            w.endln()
            return self

    def serialize(obj):
        if self.registry.registered(obj):
            return self

        partial = self.serialize_partial(obj)
        self.registry.register(obj)
        return self.serialize_references(partial, obj)


    def send(self, obj, writer = None, idnt=None):
        w0 = w = writer or self.streams[-1]

        if _is_primitive(obj):
            #w.Arg(obj)
            if isinstance(obj, np.ndarray):
                w.write(str(obj.tolist()))
                return self
            try:
                w.write(str(obj))
                return self
            except Exception as e:
                print(e, file=sys.stderr)
                return self

        if self.registry.registered(obj):
            w.write(self.registry.ident(obj).tclstr())
            return self
        else:
            ident = self.registry.register(obj)


        w.write(ident.tclstr())

        w = self.get_stream(ident.nm[0])
        w.write(f'"{ident.nm[1]}": ')
        # w.rshift()

        if _is_collection(obj, SEQUENTIAL):
            try:
                rep = "[" + ", ".join(str(i) for i in obj) + "]"
                w.write(rep)
                return self
            except Exception as e:
                print(e, file=sys.stderr)
                pass
            w.write("[")
            w.endln()
            w.rshift()
            for i in obj:
                self.send(i)
                #w.write(self.registry.ident(i).tclstr())
                w.write(", ")
            w.lshift()
            w.endln()
            w.write("],")
            w.endln()
            return self

        elif isinstance(obj, dict):
            w.write("{")
            w.endln()
            w.rshift()
            for k,v in obj.items():
                w.write(f'"{k}": ')
                self.send(v, writer=w)
                w.write(",")
                w.endln()
            w.lshift()
            w.write("},")
            w.endln()
            return self
        # if self.registry.registered(obj):
        #     return self

        elif dataclasses.is_dataclass(obj):
            obj._args = dataclasses.fields(obj)
            w.write("{")
            w.endln()
            w.rshift()
            for field in obj._args:
                value = getattr(obj, field.name)
                w.write(f'"{field.name}": ')
                self.send(value, writer=w)
                w.write(",")
                w.endln()

            w.lshift()
            w.write("},")
            w.endln()
            self.registry.register(obj)
            return self
        else:
            try:
                pass
                print(obj, file=sys.stderr)
                #w.Arg(obj)
                #return self
            except Exception as e:
                print(e, file=sys.stderr)
                # pass
                # if not hasattr(obj,"_args"):
                #     raise ObjectSerializationError(f"object {obj} of type {type(obj)}")


class JSON_Emitter(Emitter):

    #def endln(self):
    #    if not hasattr(self, "eol") or not self.eol:
    #        self.eol = True

    def Arg(this, self, value=None)->list:
        try:
            if self.__class__.__name__ in dir(this):
                getattr(this, self.__class__.__name__)(self, value=value)
                return

            try:
                value = self._get_value(None, value)
            except Exception as e:
                value = self
                # this.write(self)
                # return
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
            this.write(self.flag, str(value))
            this.endln()
        except:
            pass


    def Lst(this, self, value=None):
        this.write(self.flag)
        this.write("[")
        for v in map(self.type, value):
            this.write(v, end=", ")
        this.write("]")
        this.endln()

    def Tag(this, self, value=None):
        value = self.value if value is None else value

        if isinstance(value, int):
            force_tag = value
        else:
            force_tag = None

        ident = this.registry.register(this.current_obj, name=value, force_tag=force_tag)

        this.write(f"${ident}")
        this.endln()

    def Flg(this, self, value=None): 
        value = self.value if value is None else value
        if value: this.write(self.flag)
        this.endln()

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
            assert "min" in self.kwds
            assert len(self.args) == 1
            for v in value:
                this.Arg(self.args[0], value=v)
    
    def Ref(this, self, value=None): 
        try:
            val = self._get_value(None, value=value)
        except:
            val = value = self
        if not isinstance(val,int):
            try:
                value = "$"+this.registry.ident(value).tclstr()
            except KeyError:
                raise KeyError(f"Unable to resolve reference to {value} when writing {this.current_obj}") 
        else:
            value = val
        return this.write(value)

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

def dumps(obj, **kwds):
    stream = _Builder(JSON_Emitter).send(obj)
    return ",\n".join(f'"{k}": {{\n  {v.strm.getvalue()}}}\n  ' for k,v in stream.stream_names.items())
#.getScript(indexed=False)

def load(model):
    nodes = {
    }

