
# Basic syntactic elements
from opensees.library.ast import Tag, Grp, Ref, Flg, Num, Str, Int
from opensees.library.obj import _LineElement, Sec, Ele, Trf
# Framework objects
from opensees.library import BeamInt, Node, LinearTransform
# Frame parameters
from opensees.library import Area, Yng

Iyc = lambda **kwds: Num("iyc", field="iyc",  alt="section", about=r"Centroidal moment of inertia $I_{yy}$", **kwds)
Ixc = lambda **kwds: Num("ixc", field="ixc",  alt="section", about=r"Centroidal moment of inertia $I_{zz}$", **kwds)
Avy = lambda **kwds: Num("avy", field="avy",  alt="section", about=r"Centroidal moment of inertia $A_{yy}$", **kwds)
Avx = lambda **kwds: Num("avx", field="avx",  alt="section", about=r"Centroidal moment of inertia $A_{zz}$", **kwds)

@Ele
class PrismFrame:
    "Create a prismatic frame element."
    _call = "material, geometry, section"
    _ndms = (2, 3)
    _args = [
        Tag(),
        Grp("nodes", args=[
          Ref("iNode", type=Node,  attr="name", about=""),
          Ref("jNode", type=Node,  attr="name", about=""),
        ]),
        Area(alt="section"),
        Yng( alt="material"),
        # 3D arguments
        Num("G",    field="shear_modulus",   about="", alt="material", ndm_reqd=(3,), reqd=False),
        Num("J",    field="torsion_modulus", about="", alt="section", ndm_reqd=(3,), reqd=False),
        Iyc(ndm_reqd=(3,), reqd=False),
        # 2D/3D arguments again
        Ixc(),
        Ref("geom",  field="transform",    type=Trf, attr="name", default=LinearTransform()),
        Flg("-cMass", field="consistent_mass",
            about="Flag indicating whether to use consistent mass matrix."),
        Num("mass",field="mass_density", flag="-mass", default=0.0, reqd=False,
            about="element mass per unit length"),
        Int("delta",field="geom_flag", flag="-delta", default=0, reqd=False,
            about="Higher-order local geometry option.", enum={
                0: "No higher order effects; linear",
                1: r"include first-order $P-\delta$ effects"
            }),
        Avy(reqd=False, default=0.0),
        Avx(reqd=False, default=0.0),
    ]
    _refs=["transform"]

if False:
    ElasticBeam2D = Ele("ElasticBeam2D",
        "elasticBeamColumn",
        args = [
            Tag(),
            Grp("nodes", args=[
              Ref("iNode", type=Node,  attr="name", about=""),
              Ref("jNode", type=Node,  attr="name", about=""),
            ]),
            Area(alt="section"),
            Yng( alt="material"),
            Iyc(),
            Ixc(),
            Ref("geom",  field="transform",    type=Trf, attr="name", default=LinearTransform()),
            Num("mass",field="mass_density", flag="-mass", default=0.0, reqd=False,
                about="element mass per unit length"),
            Flg("-cMass", field="consistent_mass",
                about="Flag indicating whether to use consistent mass matrix.")
        ],
        refs=["transform"],
        alts=[
            Ref("material", type="Material"),
            Ref("section",  type=Sec)
        ],
        inherit=[_LineElement],
    )

    ElasticBeamColumn3D = Ele("ElasticBeamColumn3D",
        "elasticBeamColumn",
        args = [
            Tag(),
            Grp("nodes", args=[
              Ref("iNode", type=Node,  attr="name", about=""),
              Ref("jNode", type=Node,  attr="name", about=""),
            ]),
            Area(alt="section"),
            Yng( alt="material"),
            Num("G",    field="shear_modulus",   about="", alt="material"),
            Num("J",    field="torsion_modulus", about="", alt="section"),
            #Grp("moi", ctype="struct", args=[
              #Num("iyc", field="iyc",  about="Centroidal moment of inertia", alt="section"),
              #Num("ixc", field="ixc",  about="", alt="section"),
              Iyc(),
              Ixc(),
            #]),
            Ref("geom",  field="transform",    type=Trf, attr="name"),
            Num("mass",field="mass_density", flag="-mass", default=0.0, reqd=False,
                about="element mass per unit length"),
            Flg("-cMass", field="consistent_mass",
                about="Flag indicating whether to use consistent mass matrix.")
        ],
        refs=["transform"],
        alts=[
            Ref("material", type="Material"),
            Ref("section",  type=Sec)
        ],
        inherit=[_LineElement],
    )


