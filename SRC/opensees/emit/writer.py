
class ModelWriter:
    def __init__(self, model):
        if not hasattr(model, "header"):
            model.header = ""

        self.model = model

    def dump(self, definitions={}, **kwds):
        model = self.model
        c = self.comment_char
        ndm, ndf = model.ndm, model.ndf
        try:
            cmds  = model.header.replace("\n", f"\n{c} ") + "\n"
        except:
            cmds = ""


        head = "\n" + self.heading(1, "Initializations")
        body = self.dump_initialize(definitions=definitions, **kwds)
        if body and not body.isspace():
            cmds += "\n" + head + body

        head  = self.heading(1, "Materials")
        body  = self.dump_materials(**kwds)
        if body and not body.isspace():
            cmds += "\n" + head + body

        head  = "\n" + self.heading(1, "Elements")
        body  = self.dump_elements(**kwds)
        if body and not body.isspace():
            cmds += "\n" + head + body

        cmds += "\n" + self.heading(1, "Assemblage")
        cmds += self.dump_connectivity(**kwds)

        cmds += "\n" + self.heading(1, "Constraints")
        cmds += self.dump_constraints(**kwds)

        return "\n\n".join([cmds])

    def dump_materials(self, definitions={})->str:
        return ""

    def dump_sections(self, definitions={})->str:
        raise Exception("This model is not yet implemented by the dispatched writer")

    def dump_elements(self, *elems)->str:
        return ""

    def dump_assembly(self)->str:
        pass

    def heading(self, level, text):
        c = self.comment_char
        if level == 1:
            return f"""{c}{'='*60}\n{c} {text}\n{c}{'='*60}\n"""
        elif level == 2:
            return f"\n{c} {text}\n{c}{'-'*50}"


