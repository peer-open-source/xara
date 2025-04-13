"""
This module provides constructors for `SectionForceDeformation` objects
which represent force-deformation (or resultant stress-strain)
relationships at beam-column and plate sample points.

"""
from math import pi, sin, cos
from opensees.library import LibCmd, Cmd, Component
from opensees.library import uniaxial
from opensees.library.frame import Area, Iyc, Ixc, Yng
from opensees.library.ast import Tag, Num, Blk, Ref, Flg, Map
from . import patch


_section = LibCmd("section")

class _FiberCollection:
    def __init__(self, shapes):
        self.shapes = shapes
    def __contains__(self, point):
        return any(point in area for area in self.shapes)

@_section
class ElasticSection:
# section Elastic $secTag $E $A $Iz <$G $alphaY>
# section Elastic $secTag $E $A $Iz $Iy $G $J <$alphaY $alphaZ>
    _ndms = (2, 3)
    _args = [
        Tag(),
        Yng(),
        Area(),
        Ixc(),
        Iyc(),
        Num("G", about="Shear Modulus (optional for 2D analysis, required for 3D analysis)", ndm_reqd=(3,), reqd=False),
        Num("J", about="torsional moment of inertia of section (required for 3D analysis)", ndm_reqd=(3,), reqd=False),
        Num("alphaY", about="shear shape factor along the local $y$-axis", reqd=False),
        Num("alphaZ", about="shear shape factor along the local $z$-axis", reqd=False),
    ]

@_section
class FiberSection(_FiberCollection):
    _ndms = (2, 3)

    _args = [
        Tag(),
        Num("GJ", flag="-GJ", field="torsional_stiffness", ndm_reqd=(3,), optional=True,
            about="linear-elastic torsional stiffness $GJ$"),
        Blk("shapes", from_prop="fibers", default=[], type=Cmd, defn=dict(
           fiber=LibCmd("fiber"),
           layer=LibCmd("layer")
          )
        )
    ]
    #_refs = ["materials"]

    def __enter__(self):
        return Component.__enter__(self)

    def init(self):
        self._fibers = None
        self._area   = None
        if "material" in self.kwds:
            mat = self.kwds["material"]
            for i in self.shapes:
                if i.material is None:
                    i.material = mat

    @property
    def materials(self):
        return (f.material for f in self.fibers)

    def get_refs(self):
        return ((f.material,"uniaxialMaterial") for f in self.fibers)

    def add_patch(self, patch):
        self._fibers = None
        self.shapes.append(patch)

    def add_patches(self, patch):
        self._fibers = None
        self.shapes.extend(patch)

    def __repr__(self):
        import textwrap
        return textwrap.dedent(f"""\
        SectionGeometry
            area: {self.area}
            ixc:  {self.ixc}
            iyc:  {self.iyc}
        """)

    @property
    def patches(self):
        return [p for p in self.shapes if p.get_cmd()[0] == "patch"]

    @property
    def layers(self):
        return [p for p in self.shapes if p.get_cmd()[0] == "layer"]

    @property
    def fibers(self):
        if self._fibers is None:
            self._fibers = [
             f for a in (a.fibers if hasattr(a,"fibers") else [a] for a in self.shapes)
                for f in a
            ]
        return self._fibers

    @property
    def area(self):
        if self._area is None:
            self._area = sum(i.area for i in self.shapes)
        return self._area

    @property
    def centroid(self):
        # TODO: cache
        return sum(i.centroid * i.area for i in self.shapes) / self.area

    @property
    def ixc(self):
        # TODO: cache
        yc = self.centroid[1]
        return sum(
            p.ixc + (p.centroid[1]-yc)**2*p.area for p in self.shapes
        )

    @property
    def iyc(self):
        # TODO: cache
        xc = self.centroid[0]
        return sum(
            p.iyc + (p.centroid[0]-xc)**2*p.area for p in self.shapes
        )


    @property
    def moic(self):
        # TODO: cache
        return [
            [p.moi[i] + p.centroid[i]**2*p.area for i in range(2)] + [p.moi[-1]]
            for p in self.shapes
        ]

    def print_properties(self):
        import textwrap
        print(textwrap.dedent(f"""
        Ixx (centroid)    |   {self.ixc:9.0f}
        Iyy (centroid)    |   {self.iyc:9.0f}
        """))

@_section
class SectionAggregator:
    """
    This class groups previously-defined `UniaxialMaterial`
    objects into a single section force-deformation model.

    Each `UniaxialMaterial` object represents the section force-deformation response for a particular section degree-of-freedom (dof). There is no interaction between responses in different dof directions. The aggregation can include one previously defined section.
    """
    _img="SectionAggregator.gif"
    _args=[

        # section Aggregator $secTag $matTag1 $dof1 $matTag2 $dof2 ....... <-section $sectionTag>

        Tag(),# ,$secTag unique section tag
        Map("materials",
            about = "the force-deformation quantity to be modeled by this section object.",
            val = Ref("material", type=uniaxial, attr="name", about="tags of previously-defined `UniaxialMaterial` objects"),
            key = Flg(
                    name = "dof",
                    enum = {
                        "P":  "Axial force-deformation",
                        "Mz": "Moment-curvature about section local z-axis",
                        "Vy": "Shear force-deformation along section local y-axis",
                        "My": "Moment-curvature about section local y-axis",
                        "Vz": "Shear force-deformation along section local z-axis",
                        "T":  "Torsion Force-Deformation"
                    }
            )
        ),
        Ref("section", type=_section, flag="-section",
            about="tag of previously-defined Section object to which the UniaxialMaterial objects are aggregated as additional force-deformation relationships")
    ]
    example="""
    create new section with IDtag 2, taking the existing material tag 2 to
    represent the shear and adding it to the existing section tag 4, which
    may be a fiber section where the interaction betweeen axial force and
    flexure is already considered.

    section Aggregator 2 2 Vy -section 4;
    """

    reference="http://earthquakespectra.org/doi/abs/10.1193/1.4000136"

    authors=["Micheal H. Scott"]

