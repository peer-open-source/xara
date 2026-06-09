from dataclasses import dataclass

@dataclass
class Element:
    body_type: str
    cell_type: str
    external_nodes: list = []
    internal_nodes: list = []

    @property
    def total_nodes(self)->int:
        match self.cell_type:
            case "Q4":
                return 4
            case "Q8":
                return 8
            case "Q9":
                return 9
            case "T3":
                return 3
            case "T6":
                return 6
            case "H8":
                return 8


def cell_by_name(name, ndm):
    match name.lower():
        case "quad" | "stdquad" | "Q4":
            return Element("Plane", "Q4")
        case "tri31":
            return Element("Plane", "T3")
        case "tri6":
            return Element("Plane", "T6")
        
        case "truss":
            return Element("Truss", "L2")

        case "*Frame":
            return Element("Frame", "L2")
