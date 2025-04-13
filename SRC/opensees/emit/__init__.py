from .opensees import OpenSeesWriter
from .writer import ModelWriter
from .fedeas import FEDEAS_Writer

def serialize(ast):
    pass

class JSON(ModelWriter):
    def dump_from_tcl(self, definitions=None):
        import tempfile, os, json
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'temp.json')
            import opensees.tcl
            opensees.tcl.eval(f"""
            {OpenSeesWriter(self.model).dump(definitions=definitions)}

            print -JSON -file {path}
            """)
            with open(path,"r") as f:
                json_string = f.read()
        return json_string




    def dump(self, definitions=None):
        return self.dump_from_tcl(definitions=definitions)
        # out = {
        #     "StructuralAnalysisModel":{
        #         "properties": {"crdTransformations": [
        #             
        #         ]},
        #         "geometry": {
        #             "nodes": [
        #                 {**n, "crd": [float(i) for i in n.crd]} for n in self.model.nodes
        #             ],
        #             "elements": [
        #             ]
        #         }
        #     }
        # }
        # return json.dumps(out, indent=2)
