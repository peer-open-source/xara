

import os
from pathlib import Path
import opensees.tcl

import platformdirs
from prompt_toolkit import prompt
from prompt_toolkit.styles import Style
from prompt_toolkit.completion import NestedCompleter
from prompt_toolkit import PromptSession
from prompt_toolkit.history import FileHistory

from prompt_toolkit.lexers import PygmentsLexer
from prompt_toolkit.styles.pygments import style_from_pygments_cls
from pygments.lexers.tcl import TclLexer
import pygments.styles
from pygments.styles import get_style_by_name

from .unix import file_util_commands

try:
    style = style_from_pygments_cls(get_style_by_name('pastie'))
    # style = style_from_pygments_cls(get_style_by_name('nord'))
except:
    style = None


_banner = """\033[1m\
                 ┌─┐┌─┐┌─┐┌─┐  ┌──┌─┐┌─┐ ┌──
                 └─┘├─┘└──┘ │ ─┘  └──└───┘
 ───────────────────┘Berkeley, California ──────────────────────
                         © UC Regents
\033[0m"""

if os.name == "nt" and "PSModulePath" not in os.environ:
    # on Windows cmd, dont add ANSI escapes
    _banner = """
                 ┌─┐┌─┐┌─┐┌─┐  ┌──┌─┐┌─┐ ┌──
                 └─┘├─┘└──┘ │ ─┘  └──└───┘
 ───────────────────┘Berkeley, California ──────────────────────
                         © UC Regents
"""


completions = {
        "model": {
            "basic": {
                "1": {1},
                "2": {2, 3},
                "3": {3, 6},
            },
        },
        "node": None,
        "element": {
            "ForceBeamColumn"
        },
        "uniaxialMaterial": {
            "Elastic": None,
            "ElasticPP": None,
            "Steel02": None
        },
        "exit": None,
}

def nested_prompt(session=None, _prompt=None, inputs="", depth=1, ret=True, **kwds):
    if _prompt is None: _prompt = "  "*depth

    if session is None:
        inputs = "\n".join((inputs, prompt(_prompt, **kwds)))
    else:
        inputs = "\n".join((inputs, session.prompt(_prompt, **kwds)))

    if inputs[-1] == "{":
        while True:
            inputs = nested_prompt(session=None, inputs=inputs, depth=depth+1, ret=False, **kwds)
            if inputs[-2:] == "\n}":
                inputs = inputs + "\n"
                break

    return inputs

class OpenSeesREPL:
    """
    Create an OpenSees read-eval-print loop (REPL).
    """
    prompt = "opensees \N{WHITE PARALLELOGRAM} "

    def __init__(self, interp=None, banner=True):

        self.banner = banner

        if interp is None:
            self.interp = opensees.tcl.TclRuntime(verbose=False)
        else:
            self.interp = interp

        self.interp.eval(file_util_commands)

        try:
            history_file = Path(platformdirs.PlatformDirs("OpenSeesRT", "OpenSees").user_data_dir)
            history_file.touch()
            self.session = PromptSession(history=FileHistory(str(history_file)))

        except:
            self.session = PromptSession()


    def repl(self):
        use_vi = True
        prompt = self.prompt
        cwd_files = {}
        tcl_locals = {}

        completions.update({"source": cwd_files, "cd": cwd_files})


        interp = self.interp
        lexer  = PygmentsLexer(TclLexer)

        print(_banner)

        while True:
            cwd_files.clear()
            for file in Path(interp.eval("pwd")).glob("*.tcl"):
                cwd_files[str(file.name)] = None

            #tcl_locals.clear()
            completions.update({
                k: None for k in self.interp.eval("info commands").split() if k not in completions
            })

            inputs = nested_prompt(self.session, [('class:prompt',prompt)],
                        vi_mode=use_vi,
                        style=style,
                        lexer=lexer,
                        completer = NestedCompleter.from_nested_dict(completions),
                        complete_while_typing=False
            )

            try:
                value = interp.eval(inputs)
                if value is not None and value != "" and value != 0:
                    print(value)

            except opensees.tcl.tkinter._tkinter.TclError as e:
                pass

            except Exception as e:
                # interp.eval(f'error {{{e}}}')
                print(e)
                # pass

