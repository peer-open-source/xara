#!/usr/bin/env python3
#===----------------------------------------------------------------------===#
#
#         STAIRLab -- STructural Artificial Intelligence Laboratory
#
#===----------------------------------------------------------------------===#
#
"""
This file implements the primary command line interface for
the package which is invoked by running:

    python -m opensees

"""
import os
import sys

import opensees
import opensees.tcl

HELP = """\
usage: opensees <file> [args..]
   or  opensees [options..] <file> [args..]

Execute an OpenSees Tcl script given in <file> with optional 
script arguments [args].

Options
  -i                          After evaluating <file>, drop into
                              an interactive Tcl shell.
  -c <command>                Run <command> before executing <file>,
                              but after executing stdin, if opened.
                              This option can be given several times.
  --enable-tk                 Enable the use of the Tk library for GUIs

  -h/--help                   Print this message and exit.

  -v/--verbose


  --no-load                   Do not pre-load the OpenSeesRT library.
"""

"""
                -trans        Transient analysis
                -eigen        Eigenvalue analysis returning frequencies
                -modes        Eigenvalue analysis returning modes
                -displ        Displacement-controlled static analysis
                -force        Force-controled static analysis
                -modal        Modal response history

       opensees -print
       opensees -json
       opensees -emit
"""


def parse_args(args):
    opts = {
        "subproc":   False,
        "preload":   True,
        "verbose":   False,
        "interact":  False,
        "enable_tk": False,
        "commands":  []
    }
    file = None
    argi = iter(args[1:])
    for arg in argi:
        if arg[0] == "-":
            if arg == "-":
                file = "-"
            if arg == "-h" or arg == "--help":
                print(HELP)
                sys.exit()

            elif arg == "-modes":
                import opensees.repl.eigen
                opensees.repl.eigen.modes(*argi)
                sys.exit()

            elif arg == "-eigen":
                import opensees.repl.eigen
                opensees.repl.eigen.eigen(*argi)
                sys.exit()

            elif arg == "--enable-tk":
                opts["enable_tk"] = True

            elif arg == "--subproc":
                opts["subproc"] = True

            elif arg == "--no-load":
                # dont load OpenSees Tcl package
                opts["preload"] = False

            elif arg == "--version" or arg == "-version":
                import opensees.tcl
                print(opensees.tcl.eval("version"))
                sys.exit()

            elif arg == "-v":
                opts["verbose"] = True

            elif arg == "-c":
                opts["commands"].append(next(argi))

            elif arg == "-i":
                opts["interact"] = True


        else:
            file = arg if file is None else file
            break


    return file, opts, argi


def main():

    file, opts, argi = parse_args(sys.argv)

    tcl = opensees.tcl.TclRuntime(verbose=opts["verbose"],
                                  preload=opts["preload"],
                                  enable_tk=opts["enable_tk"])

    from_pipe = not sys.stdin.isatty()

    if file is None and len(opts["commands"]) == 0 and not from_pipe:

        if opts["subproc"]:
            OpenSeesShell().cmdloop()
            sys.exit()

        try:
            # Try full-featured REPL
            from opensees.repl.ptkshell import OpenSeesREPL
            tcl.eval("set tcl_interactive 1")
            OpenSeesREPL(interp=tcl).repl()

        except ImportError:
            from opensees.repl.cmdshell import TclShell
            TclShell().cmdloop()
        except EOFError:
            pass
        sys.exit()

    #
    #
    #
    argv = list(argi)
    tcl.eval(f"set argc {len(argv)+1}")
    tcl.eval(f"set argv {{{file} {' '.join(argv)}}}")

    script = None
    run_cmds = "before"
    if from_pipe and (file == "-" or (file is None)): # and len(opts["commands"]) == 0):
        script = sys.stdin.read()
        run_cmds = "after"

    elif file is not None:

        try:
            script = open(file).read()
        except UnicodeDecodeError:
            script = open(file, encoding="latin-1").read()

        except FileNotFoundError as e:
            print(e, file=sys.stderr)
            sys.exit(1)

        tcl.eval(f"info script {file}")
        run_cmds = "before"


    if run_cmds == "before":
        for cmd in opts["commands"]:
            tcl.eval(cmd)

    # Store the exit code of the Tcl script
    code = 0

    #
    if script is not None:
        try:
            code = tcl.eval(script)

        except opensees.tcl.InterpreterError as e:
            print(e, file=sys.stderr)
            if not opts["interact"]:
                sys.exit(-1)

    if run_cmds == "after":
        for cmd in opts["commands"]:
            tcl.eval(cmd)

    if not opts["interact"]:
        tcl.eval("wipe")

    else:
        from opensees.repl.ptkshell import OpenSeesREPL
        # TODO: do something for windows
        if os.name != "nt":
            sys.stdin = open("/dev/tty")
        tcl.eval("set tcl_interactive 1")
        try:
            OpenSeesREPL(interp=tcl).repl()
        except EOFError:
            pass

    try:
        code = int(code)
    except:
        code = 0

    sys.exit(code)


if __name__ == "__main__":
    main()

