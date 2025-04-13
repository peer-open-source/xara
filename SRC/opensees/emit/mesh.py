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
                if ("quad" in e["type"].lower() or ("shell" in e["type"].lower() and len(e["nodes"]) != 3))
            ]),
        ("triangle", [
            [nodes[int(n)] for n in e["nodes"]]
                for e in model["geometry"]["elements"]
                if "tri" in e["type"] or ("shell" in e["type"].lower() and len(e["nodes"]) == 3)
            ]),
        ("tetra", [
            [nodes[int(n)] for n in e["nodes"]]
                for e in model["geometry"]["elements"]
                if "tri" in e["type"] or ("tetrahedron" in e["type"].lower())
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
        # Optionally provide extra data on points, cells, etc.
        # point_data={"T": [0.3, -1.2, 0.5, 0.7, 0.0, -3.0]},
        # Each item in cell data must match the cells array
        # cell_data={"a": [[0.1, 0.2], [0.4]]},
    )


    return mesh

    # Alternative with the same options
    # meshio.write_points_cells("foo.vtk", points, cells)

if __name__ == "__main__":
    import sys, json

    with open(sys.argv[1]) as f:
        data = json.load(f)

    dump(data["StructuralAnalysisModel"], None).write(sys.argv[2])

