import re
import sys
import inspect
import textwrap
import importlib

from opensees.emit.emitter import Emitter, ScriptBuilder
from opensees.library import ast
from opensees.library.ast import Arg

def write_grp(a):
    args = (write_grp(i) if isinstance(i, ast.Grp) else (i.name if i.name else "") for i in a.args)
    return "["+",".join(args)+"]"

def write_obj(v, w, qual=None, ndm=None):

    args = v._args if ndm is None else filter(lambda arg: arg.reqd or ndm in arg.ndm_reqd, v._args)

    name = v.__name__
    if qual is not None:
        #qual = qual + "."
        qual = qual + "(\""
    else:
        qual = ""
#   s = "(" + ", ".join(arg.name for arg in v._args if arg.reqd) + ", **kwds)"
    s = "\", " + ", ".join(arg.name for arg in args) + ", **extra)"

    #s = str(inspect.signature(v)).replace('=None','')

    # if the signature is too long, add in line breaks
    if len(s) > 45:
        s = s.replace(", ", ",<br>&emsp;&emsp;&emsp;&emsp;")

    w.write(textwrap.dedent(f"""
    <span style="font-feature-settings: kern; color: var(--md-code-fg-color) !important; font-family: var(--md-code-font-family);">
        {qual}<span style="color:#900">{name}</span>{s}
    </span>
    """))
    try:
        w.write(v.about)
    except:
        w.write(v.__doc__ or "")
    if hasattr(v, "_img"):
        w.write(f"![](/figures/{v._img})")


    w.write(textwrap.dedent("""
    <div>
    <table>
    <!-- <colgroup>
      <col style="width: 10%" ><col style="width: 30%" ><col style="width: 60%" >
    </colgroup> -->
    <tbody>
    """))

class ApiEmitter(Emitter):
    def Arg(this, a, value=None)->list:

        if a.__class__.__name__ in dir(this):
            getattr(this, a.__class__.__name__)(a)
            return
        name = a.name or ""
        default = " = "+str(a.default) if a.default is not None else ""
        about = re.sub('[\s+]', ' ', a.about.replace('\n',' '))
        try:
            typ = f'<code style="white-space: nowrap;">{a.type.__name__}</code>'
        except:
            typ = f"<code>{a.__class__.__name__}</code>"
        if "enum" in a.kwds:
            about += "<table>" + "".join(f"<tr><td><code>\"{k}\"</code></td><td>{v}</td></tr>" for k,v in a.kwds["enum"].items()) + "</table>"

        this.write(f"<td><tt>{name+default}</tt></td><td>{typ}</td><td>{about}")


    def Lst(self, arg, value=None):
        name = arg.name or ""
        about = re.sub('[\s+]', ' ', arg.about.replace('\n',' '))
        self.write(f"<td><tt>{name}</tt></td><td><code>Array[{arg.type.__name__}]</td></code><td>{about}</td>")

    def Tag(this, self, value=None):
        pass


    def Flg(this, a, value=None):
        name = a.name or ""
        default = " = "+str(a.default) if a.default is not None else ""
        about = re.sub('[\s+]', ' ', a.about.replace('\n',' '))
        typ = f"<code>bool</code>"
        this.write(f"<td><tt>{name+default}</tt></td><td>{typ}</td><td>{about}")

    def Grp(this, a, value=None):
        name = (a and a.name) or ""
        default = " = "+str(a.default) if a.default is not None else ""
        about = re.sub('[\s+]', ' ', a.about.replace('\n',' '))
        try:
            typ = f"<code>{a.__class__.__name__}({a.type.__name__})</code>"
        except:
            typ = f"<code>{a.__class__.__name__}</code>"

        if name == "":
            name = f"<code>{write_grp(a)}</code>"
        else:
            typ = f"<code>{write_grp(a)}</code>"
        #this.write(f"<td>{name+default}</td><td>{typ}</td><td>{about}")
        this.write(f'<td colspan="2">{name+default}</td><td>{typ}, {about}')
        this.endln()
        #this.rshift()
        this.write("<table>")
        for i in a.args:
            this.endln()
            this.write("<tr>")
            this.Arg(i)
            this.write("</tr>")
        this.write("</td></table>")
        #this.lshift()

    def Ref(this, a, value=None):
        name = a.name or ""
        default = " = "+str(a.default) if a.default else ""
        about = re.sub('[\s+]', ' ', a.about.replace('\n',' '))
        if isinstance(a.type, str):
            typ = f"<code>{a.__class__.__name__}({a.type})</code>"
        else:
            typ = f"<code>{a.__class__.__name__}({a.type.__name__})</code>"
        this.write(f"<td>{name+default}</td><td>{typ}</td><td>{about}")

    # def Blk(this, self, value=None):
    #     value = self.get_value(value=value)
    #     if value is None:
    #         value = [None]
    #     this.write(self.flag, "{")
    #     this.endln()
    #     this.rshift()
    #     for v in value:
    #         this.parent.send(v)
    #     this.lshift()
    #     this.write("}")

    def Map(this, a, value=None):
        typ = f"<code>{a.__class__.__name__}</code>"
        # about = a.about.replace('\n',' ').strip()
        name = a.name or ""
        default = " = "+str(a.default) if a.default else ""
        about = re.sub('[\s+]', ' ', a.about.replace('\n',' '))
        typ = f"<code>{a.__name__}</code>"
        this.write(f"<td>{name+default}</td><td>{typ}</td><td>{about}", end="")
        this.write("<table>")
        this.Arg(a.key_type)
        this.Arg(a.val_type)
        this.write("</table>")


class ApiDocWriter(ScriptBuilder):
    def __init__(self):
        ScriptBuilder.__init__(self, ApiEmitter)

    def send(self, obj, idnt=None, qual=None, ndm=None):
        w = self.streams[0]

        if ndm is not None:
            for dm in ndm:
                write_obj(obj, w, qual=qual, ndm=dm)
                w.endln();

        else:
            write_obj(obj, w, qual=qual)
            w.endln();


        for arg in obj._args:
            w.write("<tr>")
            typ = arg.__class__.__name__

            w.Arg(arg)

            w.write("</tr>")
            w.endln();

        w.write(textwrap.dedent("""
        </tbody>
        </table>
        </div>
        """))

        w.endln();

        return self


if __name__ == "__main__":
    import opensees

    _, *module, name = sys.argv[1].split(".")
    module = ".".join(module)

    obj = getattr(getattr(opensees, module), name)

    ndm = getattr(obj, "_ndms", None)

    print(ApiDocWriter().send(obj, qual="Model."+module, ndm=ndm))





