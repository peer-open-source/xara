from collections import defaultdict

class Identifier: 
    def __init__(self, tag_space, space_name):
        self.nm = (tag_space, space_name)

    def __hash__(self):
        return hash(self.nm)

    def tclstr(self):
        return f"{self.nm[0]}({self.nm[1]})"

class TagSpace:
    def __init__(self):
        self.current_tag = 1
        self.obj_by_tag = {}
        self.tag_by_obj = {}
        self.forced_tags = set()
        self.forced_names = set()
    
    def __getitem__(self, name):
        return self.tag_by_obj[name]

    def new_tag(self):
        t = self.current_tag
        while t in self.forced_tags: t += 1
        self.current_tag = t + 1
        return t

    def add(self, name, force_tag=None)->bool:
        if force_tag is not None:
            tag = force_tag
            if force_tag in self.forced_tags:
                raise ValueError("Duplicate forced tag: '" +
                        f"{name}' and '{self.obj_by_tag[force_tag]}'"
                )
            elif tag in self.obj_by_tag:
                new_tag = self.new_tag()
                old_obj = self.obj_by_tag[tag]
                self.tag_by_obj[old_obj] = new_tag
                self.obj_by_tag[new_tag] = old_obj
            self.forced_tags.add(force_tag)
        else:
            tag = self.new_tag()

        self.tag_by_obj[name] = tag
        self.obj_by_tag[tag] = name
        return True

class Registry:
    def __init__(self):
        self.tag_spaces = defaultdict(TagSpace)
        self.identifiers = {}
        self.objects = {}
        self.anonid = 1

    def __getitem__(self, tag_space: str)->TagSpace:
        return self.tag_spaces[tag_space]

    def registered(self, obj)->bool:
        return id(obj) in self.objects

    def register(self, obj, name=None, tag_space=None, force_tag=None) -> Identifier:
        try:
            ts = tag_space or obj.tag_space
        except AttributeError:
            # references where type is given as string
            if isinstance(obj, str):
                ts = str(obj)
            else:
                ts = str(type(obj))

        if force_tag is not None:
            id2 = str(force_tag)
        elif name is None:
            # Anonymous
            assert obj is not None
            id2 = f"<{self.anonid}>"
            self.anonid += 1
        else:
            id2 = str(name)

        self.tag_spaces[ts].add(id2, force_tag=force_tag)

        ident = Identifier(ts, id2)
        self.objects[id(obj)] = ident
        self.identifiers[ident] = id(obj)
        return ident

    def index(self)->dict:
        return {
            i: self.tag_spaces[i.nm[0]][i.nm[1]]
            for i in self.identifiers if i.nm[0] is not None
        }

    def ident(self, obj)->Identifier:
        return self.objects[id(obj)]

