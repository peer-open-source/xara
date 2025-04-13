import cmd
import opensees.tcl

class TclShell(cmd.Cmd):
#   intro = """\

#   OpenSees -- Open System For Earthquake Engineering Simulation
#  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
#                       Berkeley, California
#   """

    intro = """\
    OpenSees -- Open System For Earthquake Engineering Simulation
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
            Pacific Earthquake Engineering Research Center
    """

    prompt = PROMPT
    file = None
    def __init__(self, *args, interp=None, **kwds):
        if interp is None:
            interp = opensees.tcl.TclInterpreter()
        self.tcl_interp = interp
        super().__init__(*args, **kwds)

    def default(self, arg):
        try:
            value = self.tcl_interp.eval(arg)
            if value:
                print(value)
            return None
        except Exception as e:
            print(e)

    def precmd(self, line):
        if self.file and 'playback' not in line:
            print(line, file=self.file)
        return line

    def completedefault(self, text, line, begidx, endidx):
        print(text,line,begidx,endidx)
        return ["hi"]

    def close(self):
        if self.file:
            self.file.close()
            self.file = None

