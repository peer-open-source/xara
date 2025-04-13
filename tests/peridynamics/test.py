import numpy as np
import opensees.openseespy as ops

ndim = 2
plane_type = 'e'
E = 38.4e3
nu = 0.2
dens = 2.4e-9
vel = 0.1
# ====================================================
# BUILD GEOMETRIC MODEL
# ====================================================
Lx = 10.0
Ly = 10.0
space = Lx / 40.0
delta = 3.01 * space
geom_err = 0.01 * space
totnode = (int(Lx/space)+1) * (int(Ly/space)+1)
maxfam = int((2*delta/space)**ndim)


def add_box_quad(x1, y1, x2, y2):
    xmin = min(x1, x2)
    xmax = max(x1, x2)
    ymin = min(y1, y2)
    ymax = max(y1, y2)
    ndivx = int((xmax-xmin)/space)
    ndivy = int((ymax-ymin)/space)
    x = np.linspace(xmin, xmax, ndivx+1)
    y = np.linspace(ymin, ymax, ndivy+1)
    x, y = np.meshgrid(x, y, indexing='ij')
    x = x.flatten()
    y = y.flatten()
    return x, y

x, y = add_box_quad(-0.5*Lx, -0.5*Ly, 0.5*Lx, 0.5*Ly)

model = ops.Model('basic', '-ndm', ndim)
model.eval(f"peri init {ndim} {totnode} {maxfam} {plane_type}")

for i in range(totnode):
    model.eval(f"peri node {i} {x[i]} {y[i]} 0.0")

# ====================================================
# TEST 1: COORDINATES
# model.eval(f"peri prin node {totnode-2} {totnode-1}")
# PASSED...
# ====================================================

model.eval(f"peri fam {delta}")

# ====================================================
# TEST 2: FAMILIES
# model.eval(f"peri prin fam 0 1 2 3")
# PASSED...
# ====================================================

vol = np.zeros(totnode)
for i in range(totnode):
    vol[i] = space**ndim
    model.eval(f"peri svol {i} {vol[i]}")
model.eval(f"peri cvol {space}")

# ====================================================
# TEST 3: VOLUMES
# model.eval(f"peri prin vol 0 1 2 3")
# PASSED...
# ====================================================

model.eval("peri suco")

# ====================================================
# TEST 4: SURFACE CORRECTION
# model.eval("peri prin corr 0 1 2 3")
# PASSED...
# =================================================

# kbulk = E / (3.0 * (1.0 - 2.0*nu))
# dt_max = space / np.pi / np.sqrt(kbulk / dens)
dt_max = 2.0e-8

xmin = -0.5*Lx-geom_err
xmax = -0.5*Lx+0.1*space
ymin = -0.5*Ly-geom_err
ymax = 0.5*Ly+geom_err
model.eval(f"peri boun {xmin} {ymin} {xmax} {ymax} {-vel*dt_max} 0 d")
xmin = 0.5*Lx-0.1*space
xmax = 0.5*Lx+geom_err
ymin = -0.5*Ly-geom_err
ymax = 0.5*Ly+geom_err
model.eval(f"peri boun {xmin} {ymin} {xmax} {ymax} {vel*dt_max} 0 d")

# ====================================================
# TEST 5: BOUNDARY CONDITIONS
# PASSED...
# ====================================================

tt = 0;
part = 0;
model.eval(f"peri step {tt} {part} {dt_max} {dens}")

model.eval("peri form")
# ====================================================
# TEST 6: SHAPE TENSOR FORMULATION
# PASSED...
# ====================================================
