import meshio

def load(interp, file):
    mesh = meshio.read(file)


def dump(model, file, format="vtk"):
    if "StructuralAnalysisModel" in model:
        model = model["StructuralAnalysisModel"]

    nodes = {
        int(n["name"]): i for i,n in enumerate(model["geometry"]["nodes"])
    }

    points = [
        n["crd"] for n in model["geometry"]["nodes"]
    ]

    cells = [
        ("quad", [
            [nodes[int(n)] for n in e["nodes"]]
              for e in model["geometry"]["elements"]
                if ("quad" in e["type"].lower() or ("shell" in e["type"].lower() and len(e["nodes"]))) and len(e["nodes"]) == 4
            ]),
        ("quad8", [
            [nodes[int(n)] for n in e["nodes"]]
              for e in model["geometry"]["elements"]
                if ("quad" in e["type"].lower() or ("shell" in e["type"].lower() and len(e["nodes"]))) and len(e["nodes"]) == 8
            ]),
        ("quad9", [
            [nodes[int(n)] for n in e["nodes"]]
              for e in model["geometry"]["elements"]
                if ("quad" in e["type"].lower() or ("shell" in e["type"].lower() and len(e["nodes"]))) and len(e["nodes"]) == 9
            ]),
        ("triangle", [
            [nodes[int(n)] for n in e["nodes"]]
              for e in model["geometry"]["elements"]
                if "tri" in e["type"] or ("shell" in e["type"].lower() and len(e["nodes"]) == 3)
            ]),
        ("tetra", [
            [nodes[int(n)] for n in e["nodes"]]
                for e in model["geometry"]["elements"]
                if "tet" in e["type"] or ("tetrahedron" in e["type"].lower())
            ]),
        ("hexahedron", [
            [nodes[int(n)] for n in e["nodes"]]
                for e in model["geometry"]["elements"]
                if "brick" in e["type"].lower()
            ])
    ]

    cells = [type for type in cells if len(type[1]) > 0]

    mesh = meshio.Mesh(
        points,
        cells,
    )


    return mesh

if __name__ == "__main__":
    import sys, json

    with open(sys.argv[1]) as f:
        data = json.load(f)

    dump(data["StructuralAnalysisModel"], None).write(sys.argv[2])

