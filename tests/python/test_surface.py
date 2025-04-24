import xara
import veux

mat_1 = ...
mat_2 = ...

block_1  = {
    1: (0.0, 0.0),
    2: (1.1, 0.0),
    3: (1.0, 1.0),
    4: (0.0, 1.0),
    5: (0.5,-0.1),
    6: (1.1, 0.5)
}


block_2  = {
    1: (1.1, 0.0),
    2: (2.0, 0.0),
    3: (2.0, 1.0),
    4: (1.0, 1.0),
    5: (1.5,-0.1),
#   7: (2.1, 0.5),
    8: (1.1, 0.5)
}

model = xara.Model(ndm=2, ndf=2)

model.nDMaterial("ElasticIsotropic", 1, E=29e3, nu=0.3)
model.surface((2,2), "Quad", ["1", "PlaneStress", "1"], block_1)

model.surface((2,2), "Quad", ["2", "PlaneStress", "1"], block_2)


#print(opensees.tcl.dumps(model))


veux.serve(veux.render(model))
