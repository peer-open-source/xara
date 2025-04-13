from opensees.library.ast import *
from opensees.library.obj import LibCmd

recorder = LibCmd("recorder")

class TimeSeries: pass
class Node: pass
class DOF: pass

node_responses = Str("recorder", about="a string indicating response required. Response types are given in table below.", enum={
                "disp":            "displacement*",
                "vel":             "velocity*",
                "accel":           "acceleration*",
                "incrDisp":        "incremental displacement",
                "eigen {i}":       "eigenvector for mode i",
                "reaction":        "nodal reaction",
                "rayleighForces":  "damping forces",
            }
        )

common = [
        Str("destination"),
#        Alt("destination", [
#            Str("dest_txt", flag="-file"), 
#            Str("dest_xml", flag="-xml"),
#            Str("dest_bin", flag="-binary"),
#            Grp("dest_tcp", flag="-tcp", args=[ 
#                    Str("inetAddr", about='ip address, "xx.xx.xx.xx", of remote machine to which data is sent'),
#                    Str("port",     about='port on remote machine awaiting tcp'),
#                ]
#            )
#          ],
#          reqd=True,
#          about="name of file to which output is sent."\
#                "file output is either in xml format (`-xml` option), textual (`-file` option) or binary (`-binary` option)"
#        ),

        Int("precision", reqd=False, flag="-precision",
            about="number of significant digits (default is 6)"\
                  "(optional, default: records at every time step)"
        ),

        Flg("-time", reqd=False,
            about="using this option places domain time in first entry of each data line, default is to have time omitted"
        ),

#   <-closeOnWrite> 
        Flg("-closeOnWrite", reqd=False, about="""
            using this option will instruct the recorder to invoke a
            close on the data handler after every timestep. If this is a file
            it will close the file on every step and then re-open it for the
            next step. Note, this greatly slows the execution time, but is
            useful if you need to monitor the data during the analysis."""),

        Num("time_step", flag="-dT", reqd=False,
            about="time interval for recording. will record when next step is `deltaT` greater than last recorder step."
        ),
#   <-region $regionTag> 
#       $regionTag       a region tag; to specify all nodes in the previously defined region. (optional)

]


class Recorder:
    def init(self, format=None):
        format = format or self.kwds.get("format", None)
        if format is None:
            format = self.destination.split(".")[-1]

        if format not in ["txt", "bin", "xml", "binary", "tcp"]:
            raise ValueError("Unable to deduce format")

        format = {"txt": "file", "bin": "binary"}.get(format, format)

        self._args[0].flag = "-" + format

        self._runtime = None

    def link(self, *args, runtime):
        self._runtime = runtime
        return runtime.send(self)

    def parse(self):
        if self.dest_txt:
            with open(self.dest_txt, "r") as f:
                return self.parse_txt(f)

    def update(self):
        self._data = self.parse()

    def __getitem__(self, indices):
        return self._data[self.index(indices)]



#recorder Node 
@recorder
class Node(Recorder):
    """
    Record the response of a number of nodes at every converged step.

    - In case you want to remove a recorder, you need to know the tag for that
      recorder. Here is an example on how to get the tag of a recorder:

          set tagRc [recorder Node -file nodesD.out -time -node 1 2 3 4 -dof 1 2 disp]
    """

    def parse_xml(self, stream):
        pass

    def parse_txt(self, stream):
        import numpy as np
        return np.loadtxt(stream).reshape((-1, len(self.nodes), len(self.dofs)))

    def parse_bin(self, stream):
        pass

    def _index(self, indices):
        pass

    @property
    def node(self):
        nn = len(self.nodes)
        ndf = len(self.dofs)
        return np.moveaxis(self._data, 1, 0)

    _pysn = [0, -1]
    _args = [
        *common,

        Ref("series", flag="-timeSeries", type=TimeSeries, reqd=False,
            about="the tag of a previously constructed `TimeSeries`, results from node at each time step are added to load factor from series"
        ),

#   <-node $node1 $node2...> 
       Grp("nodes", type=Ref(type="node"), reqd=False, min=1, flag="-node", about="tags of nodes whose response is being recorded (optional, default: omitted)"),

#   <-nodeRange $startNode $endNode> 
#       $startNode $endNode..    tag for start and end nodes whose response is being recorded (optional, default: omitted)

#   -dof ($dof1 $dof2 ...) 
        Grp("dofs", type=Int, reqd=False, num="*", flag="-dof", about="the specified dof at the nodes whose response is requested."),

        node_responses
    ]



@recorder
class Element:
    _args = [
            *common
    ]

# recorder Element 
#     """
#     <-file $fileName> <-xml $fileName> <-binary $fileName> 
#         $fileName	name of file to which output is sent.
#                     file output is either in xml format (-xml option), textual (-file option) or binary (-binary option)
#     <-precision $nSD> 
#         $nSD	    number of significant digits (optional, default is 6)
#     <-time> 
#         -time	    (optional using this option places domain time in first entry of each data line, default is to have time ommitted)
#     <-closeOnWrite> 
#         -closeOnWrite	optional. using this option will instruct the recorder to invoke a close on the data handler after every timestep. If this is a file it will close the file on every step and then re-open it for the next step. Note, this greatly slows the execution time, but is useful if you need to monitor the data during the analysis.
#     <-dT $deltaT> 
#         $deltaT	time interval for recording. will record when next step is $deltaT greater than last recorder step. (optional, default: records at every time step)
#     <-ele ($ele1 $ele2 ...)> 
#         $ele1 $ele2 ..	tags of elements whose response is being recorded -- selected elements in domain (optional, default: omitted)
#     <-eleRange $startEle $endEle> 
#         $startEle $endEle ..	tag for start and end elements whose response is being recorded -- range of selected elements in domain (optional, default: omitted)
#     <-region $regTag> 
#         $regTag	previously-defined tag of region of elements whose response is being recorded -- region of elements in domain (optional)
#     $arg1 $arg2 ...
#         $arg1 $arg2 ...	arguments which are passed to the `setResponse()` element method
# 

