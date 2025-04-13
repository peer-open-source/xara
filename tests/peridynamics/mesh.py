
import pygmsh
import opensees.openseespy as ops
import numpy as np

Lx = 1.0
Ly = 1.0
meshsz = 0.02
with pygmsh.geo.Geometry() as geom:
    geom.add_polygon(
        [
            [0.0, 0.0],
            [Lx, 0.0],
            [Lx, Ly],
            [0.0, Ly],
        ], mesh_size=meshsz)
    mesh = geom.generate_mesh()

# mesh.write("./tests/peridynamics/plate.vtk")

coord = mesh.points[:, :2]

ndim = 3
plane_type = 'e'
totnode = coord.shape[0]
delta = 4.01*meshsz
maxfam = int((2*delta/meshsz)**ndim)
print(totnode, maxfam, delta)

model = ops.Model('basic', '-ndm', ndim)
if ndim == 2:
    model.eval(f"peri init {ndim} {totnode} {maxfam} {plane_type}")
else:
    model.eval(f"peri init {ndim} {totnode} {maxfam}")

for i, node in enumerate(coord):
    x, y = node
    model.eval(f"peri node {i} {x} {y} 0.0")

model.eval(f"peri fam {delta}")
model.eval(f"peri prin node 0 1 2 3")
model.eval(f"peri prin fam 0 1 2 3")

conn = mesh.cells_dict['triangle'];
vol = np.zeros(totnode)
for i in range(conn.shape[0]):
    x1, y1 = coord[conn[i, 0]]
    x2, y2 = coord[conn[i, 1]]
    x3, y3 = coord[conn[i, 2]]
    v_elmt = 0.5*(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
    vol[conn[i, :]] += v_elmt / 3.0
# print(np.sum(vol))

for i in range(totnode):
    model.eval(f"peri svol {i} {vol[i]}")

model.eval(f"peri cvol {meshsz}")
model.eval(f"peri prin vol 0 1 2 3")

print(model.eval("peri form"))

#print(model.eval("peri form 8"))


# model.eval(f"peri prin node {5} {int(totnode/2)} {totnode-1}")
# model.eval(f"peri prin node 0")
# model.eval{f"peri prin node {int(totnode/2)}"}
# model.eval(f"peri prin node {totnode-1}")
# ndim = 2
# totnode = 100
# maxfam = 100

# model = ops.Model()

# model.eval(f"peri init {ndim} {totnode} {maxfam}")

# mesh = meshio.read()


# for i,node in enumerate(mesh.points):
#     x, y, z = ...
#     model.eval(f"peri node {i} {x} {y} {z}")
