# Claudio Perez
from .ast import *
from .obj import *
from .obj import _LineElement
from .model import Node

def Dof(about=""):
    return Int("dof", about=about, enum={
        "1": "corresponds to translation along the global 1 axis",
        "2": "corresponds to translation along the global 2 axis",
        "3": "corresponds to translation along the global 3 axis",
        "4": "corresponds to rotation about the global 1 axis",
        "5": "corresponds to rotation about the global 2 axis",
        "6": "corresponds to rotation about the global 3 axis",
    })

def Yng(about=None, **kwds):
    about = about or "Young's modulus of elasticity"
    return Num("E", field="elastic_modulus", about=about, **kwds)

def Yld(nature="strength", about=None, **kwds):
    about = about or f"Yield {nature}."
    return Num("fy", field=f"yield_{nature}", about=about, **kwds)

def Area(**kwds):
    return Num("A", field="area", about="cross-sectional area", **kwds)

def Izz(**kwds):
    return Num("Iz", field="izz", about="area moment of inertia", **kwds)

Rho = lambda: None
class ConstitutiveData:
    _args = [Yng(), Yld(), Rho()]


class BeamInt(Arg):
    def as_tcl_list(self, value):
        args = []
        for kk in ["rule", "section", "num"]:
            k = value[kk]
            if isinstance(k,Arg):
                args.append(k._get_value())
            else:
                args.append(k)
        #return ['"' + ' '.join(args) + '"']
        return args
        #return f"{self.kwds['rule']} {self.kwds['section']} {self.kwds['num']}"

LinearTransform = Trf("LinearTransform", "Linear",
  args = [
    Tag("name", about="Tag used to identify the transform object."),
    Grp("vecxz", type=Num, num=3, ndm=3),
    Grp("joint_offsets", reqd=False, flag="-jntOffset", type=Grp, args=[
        Grp(type=Num, num=3, default=[0.0]*3),
        Grp(type=Num, num=3, default=[0.0]*3)
    ])
  ],
)

CorotTransform = Trf("CorotTransform", "Corotational",
  args = [
    Tag("name", about="Tag used to identify the transform object."),
    Grp("vecxz", type=Num, num=3, ndm=3),
    Grp("joint_offsets", reqd=False, flag="-jntOffset", type=Grp, args=[
        Grp(type=Num, num=3, default=[0.0]*3),
        Grp(type=Num, num=3, default=[0.0]*3)
    ])
  ],
)

class uniaxial:
    """
    `UniaxialMaterial` object library.

    A `UniaxialMaterial` object typically represents a pair of work conjugate
    scalars such as axial stress/strain, moment/cuvature, or force/deformation.
    """
    class UniaxialWrapper:
        def __init__(self, wrapped):
            self._wrapped = wrapped
            self._wrapping_several = isinstance(self._wrapped, (tuple, list))

        def __repr__(self):
            return f"{self.__class__.__name__}<{self._wrapped.__class__.__name__}>"

        def __enter__(self):
            if self._wrapping_several:
                self.wrapped = type(self._wrapped)(
                    w.__enter__() for w in self._wrapped
                )
            else:
                self.wrapped = self._wrapped.__enter__()

            return self

        def __exit__(self, *args, **kwds):
            if self._wrapping_several:
                for i,w in enumerate(self._wrapped):
                    w.__exit__(*args, **kwds)
                    del self.wrapped[i]
            else:
                self._wrapped.__exit__(*args, **kwds)
            del self.wrapped

    @classmethod
    def create(cls, name, tag, *args):
        obj = Uni(name, name, args=[Tag()]+[Str(f"arg{i}") for i,arg in enumerate(args)])
        return obj(tag, *map(str,args))

    Series = Uni("SeriesMaterial", "Series", args = [
            Tag(),
            Grp("materials", min=2, type=Ref(type="uniaxial", attr="name"))
        ],
        refs = ["materials"]
    )

    Parallel = Uni("ParallelMaterial", "Parallel", args = [
            Tag(),
            Grp("materials", min=2, type=Ref(type="uniaxial", attr="name"))
        ],
        refs = ["materials"]
    )


    Elastic = ElasticSpring = Uni("ElasticUniaxialMaterial", "Elastic",
        args = [
            Tag(),
            Yng(),
            Num("eta",  field="damp_tangent",     reqd=False, default=0.0, about="damping tangent"),
            Num("Eneg", field="negative_modulus", reqd=False),
            Num("G",    field="shear_modulus",    reqd=False),
        ]
    )

    ElasticBilin = Uni("ElasticBilin", "ElasticBilin", args=[
        Tag(),
        Yng(about=r"tangent in tension for stains: 0 <= strains $\le$ `epsP2`"),
        Num("E2", about="tangent when material in tension with strains > `epsP2`"),
        Num("eps2", about="strain at which material changes tangent in tension."),
        Num("EN1", reqd=False, about=r"optional, default = `EP1`. tangent in compression for stains: 0 < strains $\le$ `epsN2`"),
        Num("EN2", reqd=False, about="optional, default = `EP2`. tangent in compression with strains < `epsN2`"),
        Num("epsN2", reqd=False, about="optional, default = -epsP2. strain at which material changes tangent in compression.")
    ])

    ElasticMultilinear = Uni("ElasticMultilinear", "ElasticMultilinear", args=[
            Tag(),
            Grp("points", type=Grp, args=[Grp(type=Num, args=[Num("-strain"), Num("-stress")])],
                about="Points defining multilinear backbone."
            ),
        ]
    )

    Hardening = Uni("HardeningMaterial", "Hardening", args=[
        Tag(),
        Yng(),
        Yld("stress"),
        Num("H_iso", about="isotropic hardening Modulus"),
        Num("H_kin", about="kinematic hardening Modulus"),
        Num("eta",   reqd=False, about="visco-plastic coefficient (optional, default=0.0)"),
    ])

    ElasticPP = Uni("ElasticPP", "ElasticPP",
        about="This command is used to construct an elastic perfectly-plastic uniaxial material object.",
        args=[
          Tag(),
          Yng(),
          Num("epsyP",  about="strain or deformation at which material reaches plastic state in tension"),
          Num("epsyN", reqd=False, about="strain or deformation at which material reaches plastic state in compression. (optional, default is tension value)"),
          Num("eps0",  reqd=False, about="initial strain (optional, default: zero)")
    ])

    Hysteretic = Uni("Hysteretic", "Hysteretic", args = [
        Tag(),
        Grp("points", args = [
          Grp(args=[Num("s1p"), Num("e1p")], about="stress and strain (or force & deformation) at first point of the envelope in the positive direction"),
          Grp(args=[Num("s2p"), Num("e2p")], about="stress and strain (or force & deformation) at second point of the envelope in the positive direction"),
          Grp(args=[Num("s3p"), Num("e3p")], about="stress and strain (or force & deformation) at third point of the envelope in the positive direction (optional)"),
          Grp(args=[Num("s1n"), Num("e1n")], about="stress and strain (or force & deformation) at first point of the envelope in the negative direction"),
          Grp(args=[Num("s2n"), Num("e2n")], about="stress and strain (or force & deformation) at second point of the envelope in the negative direction"),
          Grp(args=[Num("s3n"), Num("e3n")], about="stress and strain (or force & deformation) at third point of the envelope in the negative direction (optional)"),
        ]),
        Num("pinchx", about="pinching factor for strain (or deformation) during reloading"),
        Num("pinchy", about="pinching factor for stress (or force) during reloading"),
        Num("damage1", reqd=False, about="damage due to ductility: D1(mu-1)"),
        Num("damage2", reqd=False, about="damage due to energy: D2(Eii/Eult)"),
        Num("beta", about="power used to determine the degraded unloading stiffness based on ductility, mu-beta (optional, default=0.0)")
        ])

    RambergOsgoodSteel = Uni("RambergOsgoodSteel", "RambergOsgoodSteel", args=[
        Tag(), Yng(), Yld("stress"), Num("a", about="yield offset parameter"),
        Num("n", about="Parameter to control the transition from elastic to plastic branches. Additionally controls the hardening of the material: by increasing $n$, hardening ratio will be decreased.")
    ])

    DoddRestrepo = Uni("Dodd_Restrepo", "Dodd_Restrepo", args=[
        Tag(),
        Yld("stress"),
        Num("fu",  about="Ultimate tensile strength (UTS)"),
        Num("esh", about="Tensile strain at initiation of strain hardening"),
        Num("esu", about="Tensile strain at the UTS"),
        Yng(),
        Num("eshi", about="Tensile strain for a point on strain hardening curve, recommended range of values for eshI: " \
                r"$\left[ \frac{\texttt{esu} + 5 \texttt{esh}}{6}, \frac{\texttt{esu} + 3 \texttt{esh}}{4}\right]$"),
        Num("fshi", about="Tensile stress at point on strain hardening curve corresponding to eshI"),
        Num("OmegaFac", default=1.0, reqd=False, about=r"Roundedness factor for Bauschinger curve in cycle reversals from the strain hardening curve. Range: $\left[0.75, 1.15\right]$. Largest value tends to near a bilinear Bauschinger curve. Default = 1.0."),
    ])


    Steel02 = Uni("Steel02",
        "Steel02",
        about="""
        This command is used to construct a uniaxial
        Giuffre-Menegotto-Pinto steel material object
        with isotropic strain hardening.""",
        args = [
          Tag(),
          Num("fy", about="yield strength"),
          Yng(),
          Num("b", default=0.0, about="strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent"),

          Num("R0", reqd=False, default=18),
          Num("cR1", reqd=False, default=0.925, about="Bauschinger evolution parameter"), 
          Num("cR2", reqd=False, default=0.150, about="Bauschinger evolution parameter"),

          # Grp("a", reqd = False, about="isotropic hardening parameters", args=[
          Num("a1", reqd=False, default=0.0, about="""
                    increase of compression yield envelope as proportion
                    of yield strength after a plastic strain of `a2∗(Fy/E0)`"""
          ),
          Num("a2", reqd=False, default=1.0, about="see explanation under `a1`."),
          Num("a3", reqd=False, default=0.0, about="""
                    increase of tension yield envelope as proportion
                    of yield strength after a plastic strain of `a4∗(Fy/E0)`"""
          ),
          Num("a4", reqd=False, default=1.0, about="see explanation under `a3`."),

          Num("sigInit", reqd=False, default=0.0, about="initial stress")
          #]),
        ],
        # params (list (float)) parameters to control the transition from 
        # elastic to plastic branches. params=[R0,cR1,cR2].
    )

    SteelMPF = Uni("SteelMPF",
        "SteelMPF",
        about="""
        This command is used to construct a uniaxial
        Giuffre-Menegotto-Pinto steel material object
        with isotropic strain hardening.""",
        args = [
          Tag(),
          Num("fyp", about="yield strength"),
          Num("fyn", about="yield strength"),
          Yng(),
          Num("bp", default=0.0, about="strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent"),
          Num("bn", default=0.0, about="strain-hardening ratio (ratio between post-yield tangent and initial elastic tangent"),

          Num("R0", reqd=False, default=18),
          Num("cR1", reqd=False, default=0.925, about="Bauschinger evolution parameter"),
          Num("cR2", reqd=False, default=0.150, about="Bauschinger evolution parameter"),
          # Grp("a", reqd = False, about="isotropic hardening parameters", args=[
          Num("a1", reqd=False, default=0.0, about="""
                    increase of compression yield envelope as proportion
                    of yield strength after a plastic strain of `a2∗(Fy/E0)`"""
          ),
          Num("a2", reqd=False, default=1.0, about="see explanation under `a1`."),
          Num("a3", reqd=False, default=0.0, about="""
                    increase of tension yield envelope as proportion
                    of yield strength after a plastic strain of `a4∗(Fy/E0)`"""
          ),
          Num("a4", reqd=False, default=1.0, about="see explanation under `a3`."),

          Num("sigInit", reqd=False, default=0.0, about="initial stress")
          #]),
        ],
        # params (list (float)) parameters to control the transition from 
        # elastic to plastic branches. params=[R0,cR1,cR2].
        # Recommended values: R0=between 10 and 20, cR1=0.925, cR2=0.15"),
    )

    ReinforcingSteel = Uni("ReinforcingSteel", "ReinforcingSteel", args=[
        Tag(),
        Yld("stress", about="Yield stress in tension (see Figure 1)"),
        Num("fu",     about="Ultimate stress in tension"),
        Yng(          about="Initial elastic tangent"),
        Num("Esh",    about="Tangent at initial strain hardening"),
        Num("esh",    about="Strain corresponding to initial strain hardening"),
        Num("eult",   about="Strain at peak stress"),
        Grp("-GABuck",reqd=False, about="Buckling Model Based on Gomes and Appleton (1997)", args=[
            Num("lsr", about="Slenderness Ratio (see Figure 2)"),
            Num("beta", about="Amplification factor for the buckled stress strain curve. (see Figure 3)"),
            Num("r", about="""Buckling reduction factor
                r can be a real number between [0.0 and 1.0]
                r=1.0 full reduction (no buckling)
                r=0.0 no reduction
                0.0<r<1.0 linear interpolation between buckled and unbuckled curves
                """),
            Num("gamma", about="Buckling constant (see Figures 3 and 4)")
            ]
          ),


        Grp("-DMBuck", reqd=False, about="Buckling model based on Dhakal and Maekawa (2002)", args=[
            Num("lsr", about="Slenderness Ratio (see Figure 2)"),
            Num("alpha", about="Adjustment Constant usually between 0.75 and 1.0", default=1.0),
          ]
        ),

        Grp("-CMFatigue", reqd=False, about="Coffin-Manson Fatigue and Strength Reduction", args=[
            Num("Cf", about="Coffin-Manson constant C (see Figure 5)"),
            Num("alpha", about="Coffin-Manson constant a (see Figure 5)"),
            Num("Cd", about="Cyclic strength reduction constant (see Figure 6 and Equation 3)"),
          ]
        ),

        Grp("-IsoHard", reqd=False, about="Isotropic Hardening / Diminishing Yield Plateau", args=[
            Num("a1", about="Hardening constant", default=4.3),
            Num("limit", about="Limit for the reduction of the yield plateau. "\
                               "% of original plateau length to remain (0.01 < limit < 1.0 )."\
                               "If `limit` = 1.0, then no reduction takes place", default=0.01),
          ]
        ),

        Grp("-MPCurveParams", about="Menegotto and Pinto Curve Parameters see Fig 6b", args=[
            Num("R1", default = 0.333),
            Num("R2", default = 18),
            Num("R3", default = 4)
          ]
        )
    ])

    Steel04 = Uni("Steel04", "Steel4", args=[
        Tag(),
        Yld("stress"),
        Yng(),
        Grp("-kin", reqd=False, type=Num, about="""apply kinematic hardening
            Kinematic hardening is based on the Menegotto-Pinto model. The parameters and their use is identical to those of the Steel02 material.

            Steel4 param kin.png

            recommended values: `R_0 = 20`, `r_1 = 0.90`, `r_2 = 0.15`""", args=[

            Num("b_k", about="hardening ratio (E_k/E_0)"),
            Num("R_0", about="control the exponential transition from linear elastic to hardening asymptote"),
            Num("r_1"),
            Num("r_2"),
        ]),

        Grp("-iso", reqd=False, about="""apply isotropic hardening
            Isotropic hardening increases the yield strength of the material.
            The applied increase is calculated as a function of the accumulated
            plastic strain. The following parameters control that function.

            Steel4 param iso.png""", args=[

            Num("b_i", about="initial hardening ratio (E_i/E_0)"),
            Num("b_l", about="saturated hardening ratio (E_is/E_0)"),
            Num("rho_i", about="specifies the position of the intersection point between initial and saturated hardening asymptotes"),
            Num("R_i", about="control the exponential transition from initial to saturated asymptote"),
            Num("l_yp", about=r"length of the yield plateau in $\varepsilon_{y0} = f_y / E_0$ units"),
        ]),

        Grp("-ult", type=Num, reqd=False, about="""apply an ultimate strength limit
            The ultimate strength limit serves as an upper limit of material
            resistance. After the limit is reached the material behaves in a
            perfectly plastic manner. Exponential transition is provided from
            the kinematic hardening to the perfectly plastic asymptote. Note
            that isotropic hardening is also limited by the ultimate strength,
            but the transition from the isotropic hardening to the perfectly
            plastic asymptote is instantaneous.

            Steel4 param ult.png
            """, args=[
                Num("f_u", about="ultimate strength"),
                Num("R_u", about="control the exponential transition from kinematic hardening to perfectly plastic asymptote"),
        ]),

        Grp("-asym", reqd=False, type=Grp, about="""assume non-symmetric behavior
                If non-symmetric behavior is assumed, material response under tension and compression will be controlled by two different parameter sets. The normal parameters control behavior under tension. Additional parameters shall be specified to describe behavior under compression. The following parameters are expected after the normal parameters when the options below are used.
                Steel4 param asymi.png
                Steel4 param asymk.png""", args=[
                Grp("-kin", type=Num, num=4, args=[Num("b_kc"), Num("R_0c"), Num("r_1c"), Num("r_2c")]),
                Grp("-iso", type=Num, num=4, args=[Num("b_ic"), Num("rho_ic"), Num("b_lc"), Num("R_ic")]),
                Grp("-ult", type=Num, num=2, args=[Num("f_uc"), Num("R_uc")]),
            ]
        ),

        Num("sig_init", flag="-init", reqd=False, about="""apply initial stress
            Initial stress is assumed at 0 strain at the beginning of the
            loading process. The absolute value of the initial stress is
            assumed to be less than the yield strength of the material.

            Steel4 param init.png"""
        ),

        Int("cycNum", flag="-mem", reqd=False, about="""expected number of half-cycles during the loading process
            Efficiency of the material can be slightly increased by correctly setting this value. The default value is `cycNum = 50`
            Load history memory can be turned off by setting `cycNum = 0`."""
        ),

    ])

    UVCuniaxial = UVC = Uni("UVCuniaxial", "UVCuniaxial", args = [
        Tag(),
        Yng(),
        Yld("stress"),
        Num("QInf",   about="Maximum increase in yield stress due to cyclic hardening (isotropic hardening)."),
        Num("b",      about="Saturation rate of QInf, b > 0."),
        Num("DInf",   about="Decrease in the initial yield stress, to neglect the model updates set DInf = 0."),
        Num("a",      about="Saturation rate of DInf, $a > 0$. If $D_{\\infty} == 0$, then a is arbitrary (but still a > 0)."),
        Int("N",      about="Number of backstresses to define, $N \\ge 1$."),
        Grp("backstress", type=Grp, args=[Grp(type=Num, args=[Num("C"), Num("gamma")])],
            about="Backstress parameters, up to 9 pairs may be specified. If `C` is specified, then the corresponding `gamma` must also be specified. Note that only the first N backstresses will be read by the parser."
        ),
    ])

    ConfinedConcrete01  = Uni("ConfinedConcrete01", "ConfinedConcrete01", args=[
         Tag(),
         Str("secType", about="label for the transverse reinforcement configuration.", enum={
             'S1' : "square section with S1 type of transverse reinforcement with or without external FRP wrapping",
             'S2' : "square section with S2 type of transverse reinforcement with or without external FRP wrapping",
             'S3' : "square section with S3 type of transverse reinforcement with or without external FRP wrapping",
             'S4a': "square section with S4a type of transverse reinforcement with or without external FRP wrapping",
             'S4b': "square section with S4b type of transverse reinforcement with or without external FRP wrapping",
             'S5' : "square section with S5 type of transverse reinforcement with or without external FRP wrapping",
             'C'  : "circular section with or without external FRP wrapping",
             'R'  : "rectangular section with or without external FRP wrapping."
         }),
         Num("fpc"),
         Yng(),
         Alt("eps_ult", [
               Num("epscu", flag="-epscu", about="confined concrete ultimate strain"),
               Num("gamma", flag="-gamma",
                     about="the ratio of the strength corresponding to "\
                         "ultimate strain to the peak strength of the "\
                         "confined concrete stress-strain curve. "\
                         "If gamma cannot be achieved in the range "\
                         "`[0, epscuLimit]` then epscuLimit "\
                         "(optional, default: 0.05) will be assumed as ultimate strain.")
            ]
          ),
          Alt("poisson", [
                Num("nu",      flag="-nu"),
                Flg("varub",   flag="-varub"),
                Flg("varnoub", flag="-varnoub"),
            ]
          ),
          Num("L1",
              about="length/diameter of square/circular core section measured respect to the hoop center line."),
          Num("L2",  reqd=False,  about="additional dimensions when multiple hoops are being used."),
          Num("L3",  reqd=False,  about="additional dimensions when multiple hoops are being used."),
          Num("phis",
              about="hoop diameter. If section arrangement has multiple hoops it refers to the external hoop."),
          Num("S",         about="hoop spacing."),
          Num("fyh",       about="yielding strength of the hoop steel."),
          Num("Es0",       about="elastic modulus of the hoop steel."),
          Num("haRatio",   about="hardening ratio of the hoop steel."),
          Num("mu",        about="ductility factor of the hoop steel."),
          Num("phiLon",    about="diameter of longitudinal bars."),
          Grp("internals", flag="-internal", reqd=False,
              about="optional parameters for defining the internal transverse reinforcement."\
                    "If they are not specified they will be assumed equal to the external ones "\
                    "(for S2, S3, S4a, S4b and S5 typed).",
              args = [
                  Num("phisi"), Num("Si"), Num("fyhi"), Num("Es0i"), Num("haRatioi"), Num("mui"),
              ]
          ),
          Grp("wrap", flag="-wrap", reqd=False, args=[
              Num("cover", about="cover thickness measured from the outer line of hoop."),
              Num("Am",    about="total area of FRP wraps (number of layers x wrap thickness x wrap width)."),
              Num("Sw",    about="spacing of FRP wraps (if continuous wraps are used the spacing is equal to the wrap width)."),
              Num("ful",   about="ultimate strength of FRP wraps."),
              Num("Es0w",  about="elastic modulus of FRP wraps."),
          ]),
          Flg("gravel", reqd=False),
          Flg("silica", reqd=False),
          Num("tol"       , flag="-tol",         reqd=False),
          Num("maxNumIter", flag="-maxNumIter",  reqd=False),
          Num("epscuLimit", flag="-epscuLimit",  reqd=False),
          Num("stRatio"   , flag="-stRatio",     reqd=False),
    ])

    Concrete02 = Uni("Concrete02", "Concrete02",
        args=[
          Tag(about="integer tag identifying material"),
          Num("fpc",   about="concrete compressive strength at 28 days (compression is negative)"),
          Num("epsc0", about="concrete strain at maximum strength"),
          Num("fpcu",  about="concrete crushing strength"),
          Num("epsU",  about="concrete strain at crushing strength"),
          Num("lamda", about="ratio between unloading slope at `epscu` and initial slope"),
          Num("ft",    about="tensile strength"),
          Num("Ets",   about="tension softening stiffness (absolute value) " +
                             "(slope of the linear tension softening branch)"),
        ],
        about="""This command is used to construct a uniaxial Kent-Scott-Park
            concrete material object with degraded linear unloading/reloading
            stiffness according to the work of Karsan-Jirsa and no tensile
            strength. (REF: Fedeas)."""
        )

    Concrete04 =  Uni("Concrete04", "Concrete04", args=[
        Tag(),
        Num("fc"  ,
            about="floating point values defining concrete compressive strength "\
                  "at 28 days (compression is negative)*"),
        Num("ec"  ,  about="floating point values defining concrete strain at maximum strength*"),
        Num("ecu" ,  about="floating point values defining concrete strain at crushing strength*"),
        Num("Ec"  ,  about="floating point values defining initial stiffness**"),
        Grp("tension", reqd=False, type=Num, args=[
            Num("fct" ,  about="floating point value defining the maximum tensile strength of concrete"),
            Num("et"  ,  about="floating point value defining ultimate tensile strain of concrete"),
        ]),
        Num("beta",  reqd=False, about="floating point value defining the exponential curve parameter to define the residual stress (as a factor of $ft) at $etu"),
    ])

    Concrete02IS = Uni("Concrete02IS", "Concrete02IS", [
         Tag(),
         Yng(),
         Num("fpc", about="peak compressive stress"),
         Num("epsc0"),
         Num("fpcu"),
         Num("epscu"),
         Grp("tension", reqd=False, type=Num, args=[
            Num("rat"), Num("ft"), Num("Ets", about="Tensile initial modulus")
         ])
    ])

    Concrete06 = Uni("Concrete06", "Concrete06", args=[
        Tag("matTag", about="integer tag identifying material"),
        Num("fc", about="concrete compressive strength (compression is negative)*"),
        Num("e0", about="strain at compressive strength*"),
        Num("n", about="compressive shape factor"),
        Num("k", about="post-peak compressive shape factor"),
        Num("alpha1", about=r"$\alpha_1$ parameter for compressive plastic strain definition"),
        Num("fcr", about="tensile strength"),
        Num("ecr", about="tensile strain at peak stress (`fcr`)"),
        Num("b", about="exponent of the tension stiffening curve"),
        Num("alpha2", about=r"$\alpha_2$ parameter for tensile plastic strain definition"),
    ])


    ConcreteCM = Uni("ConcreteCM", "ConcreteCM", args=[
        Tag(), #"mattag  Unique uniaxialMaterial tag"),
        Num("fpcc" , about=r"Compressive strength (f'c)"),
        Num("epcc" , about=r"Strain at compressive strength ($\epsilon^\prime_c$)"),
        Num("Ec"   , about=r"Initial tangent modulus ($E_c$)"),
        Num("rc"   , about=r"Shape parameter in Tsai’s equation defined for compression ($r_c$)"),
        Num("xcrn" , about=r"Non-dimensional critical strain on compression envelope ($\epsilon^-_{cr}$, where the envelope curve starts following a straight line)"),
        Num("ft"   , about=r"Tensile strength ($f_t$)"),
        Num("et"   , about=r"Strain at tensile strength ($\epsilon$t)"),
        Num("rt"   , about=r"Shape parameter in Tsai’s equation defined for tension ($r_t$)"),
        Num("xcrp" , about=r"Non-dimensional critical strain on tension envelope ($\epsilon^+_{cr}$, where the envelope curve starts following a straight line – large value [e.g., 10000] recommended when tension stiffening is considered)"),
        Int("gap",   reqd=False, flag="-GapClose",  about="gap = 0, less gradual gap closure (default); gap = 1, more gradual gap closure"),
    ])


    #
    # MOMENT-CURVATURE
    #

    BoucWen = Uni("BoucWen", "BoucWen", args=[
        Tag(),
        Num("alpha", about=r"ratio of post-yield stiffness to the initial elastic stiffenss ($0< \alpha < 1$)"),
        Num("ko",    about="initial elastic stiffness"),
        Int("n",     about="parameter that controls transition from linear to nonlinear range (as n increases the transition becomes sharper; n is usually grater or equal to 1)"),
        Grp("a",     args=[Num("gamma"), Num("beta")], about="parameters that control shape of hysteresis loop; depending on the values of γ and β softening, hardening or quasi-linearity can be simulated (look at the NOTES)"),
        Grp("b",     args=[Num("Ao"), Num("deltaA")], about="parameters that control tangent stiffness"),
        Grp("c",     args=[Num("deltaNu"), Num("deltaEta")], about="parameters that control material degradation"),
    ])


    #
    # Soil
    #
    PySimple1 = Uni("PySimple1", "PySimple1", args=[
        Tag(),
        Int("soilType", enum={
             "1": "Backbone of $p-y$ curve approximates Matlock (1970) soft clay relation.",
             "2": "Backbone of $p-y$ curve approximates API (1993) sand relation."
            }
        ),
        Num("pult", about="Ultimate capacity of the p-y material. Note that $p$ or $p_{ult}$ are distributed loads [force per length of pile] in common design equations, but are both loads for this uniaxialMaterial [i.e., distributed load times the tributary length of the pile]."),
        Num("y50", about="Displacement at which 50% of pult is mobilized in monotonic loading."),
        Num("Cd", about="Variable that sets the drag resistance within a fully-mobilized gap as $C_d p_{ult}$."),
        Num("c", reqd=False, about="The viscous damping term (dashpot) on the far-field (elastic) component of the displacement rate (velocity). (optional Default = 0.0). Nonzero c values are used to represent radiation damping effects"),
    ])

    QzSimple1 = Uni("QzSimple1", "QzSimple1", args=[
        Tag(),
        Str("qzType", enum={
                  "1": "Backbone of $q-z$ curve approximates Reese and O'Neill's (1987) relation for drilled shafts in clay.",
                  "2": "Backbone of $q-z$ curve approximates Vijayvergiya's (1977) relation for piles in sand."
                }
        ),

        Num("qult", about="Ultimate capacity of the $q-z$ material. SEE NOTE 1."),
        Num("Z50", about="Displacement at which 50% of qult is mobilized in monotonic loading. SEE NOTE 2."),
        Num("suction", reqd=False, about="Uplift resistance is equal to suction*qult. Default = 0.0. The value of suction must be 0.0 to 0.1.*"),
        Num("c", reqd=False, about="The viscous damping term (dashpot) on the far-field (elastic) component of the displacement rate (velocity). Default = 0.0. Nonzero c values are used to represent radiation damping effects.*"),
    ])


    # Clough = Uni("Clough", "SNAP::Clough", args=[
    #     Tag(),
    #     #Num("elstk",     about="initial elastic stiffness"),
    #     Num("fyieldPos", about="positive yield strength"),
    #     Num("fyieldNeg", about="yield strength in compression"),
    #     Num("alpha",     about="strain hardening ratio (fraction of elstk)"),
    #     Num("elstk",     about="initial elastic stiffness"),
    #     Num("Resfac	",   about="residual stress after collapse"),
    #     Num("dyieldPos", about="positive yield displacement"),
    #     Num("dyieldNeg", about="negative yield displacement"),
    #         
    #     Grp("a", args=[Num("ecaps"), Num("ecapd"), Num("ecapk"), Num("ecapa")], about="parameter expressing the hystetic energy dissipacion capacity."),
    #     Grp("b", args=[Num("Enrgts"), Num("Enrgtd"), Num("Enrgtk"), Num("Enrgta")], about="hysteretic energy dissipation capacity")
    # ])

class element:

    ZeroLength = ZeroLength3D = Ele("ZeroLength3D",
        "zeroLength",
        args = [
            Tag(),
            Grp("nodes", args=[
                Ref("iNode", type=Node,  attr="name", about=""),
                Ref("jNode", type=Node,  attr="name", about=""),
            ]),
            Grp("materials", flag="-mat", type=Ref(type=uniaxial, attr="name"), num=6),
            Grp("dofs",      flag="-dir", type=Dof, num=6, default=[
                  "$dx", "$dy", "$dz", "$rx", "$ry", "$rz"
            ]),
            Grp("orientation", reqd=False, flag="-orient", args=[
                Grp("x",  type=Num, num=3),
                Grp("yp", type=Num, num=3)
            ]),
            Int("do_rayleigh", flag="-doRayleigh", default=0)
        ],
        refs=["materials"]
    )

    @Ele
    class Truss:
        _args = [
            Tag(),
            Grp("nodes", type=Node,
                args=[Ref("iNode", attr="name"), Ref("jNode", attr="name")], about="end nodes"),
            Num("A", about="cross-sectional area of element"),
            Ref("material", attr="name", type=Uni, about="tag associated with previously-defined UniaxialMaterial"),
            #Ref("section", type="Section", about="tag associated with previously-defined Section"),
            Num("rho",   reqd=False, about="mass per unit length, optional, default = 0.0"),
            Int("cFlag", reqd=False, about="consistent mass flag, optional, default = 0"+
                "0: lumped mass matrix (default)"+
                "1: consistent mass matrix"),
            Int("rFlag", reqd=False, about="Rayleigh damping flag, optional, default = 0", enum={
                    "0": "NO RAYLEIGH DAMPING (default)",
                    "1": "include Rayleigh damping",
                }
            )
        ]
        _refs=["material"]

    @Ele
    class forceBeamColumn:
        "Create a forceBeamColumn element."
        _call = "material, geometry, section"
        _args=[
            Tag(),
            Grp("nodes", args=[
              Ref("iNode", type=Node,  attr="name", about=""),
              Ref("jNode", type=Node,  attr="name", about=""),
            ]),
            Ref("geom",  field="transform", type=Trf, attr="name"),
            BeamInt("integration"),
            Flg("-cMass", field="consistent_mass",
                about="Flag indicating whether to use consistent mass matrix."),
            Num("mass",field="mass_density", flag="-mass", default=0.0, reqd=False, 
                about="element mass per unit length"),
        ]
        _refs=["transform", "section"]

        @property
        def section(self):
            return self.integration["section"]


    @Ele
    class forceBeamColumn:
        "Create a forceBeamColumn element."
        _call = "material, geometry, section"
        _args=[
            Tag(),
            Grp("nodes", args=[
              Ref("iNode", type=Node,  attr="name", about=""),
              Ref("jNode", type=Node,  attr="name", about=""),
            ]),
            Ref("geom",  field="transform", type=Trf, attr="name"),
            BeamInt("integration"),
            Flg("-cMass", field="consistent_mass",
                about="Flag indicating whether to use consistent mass matrix."),
            Num("mass",field="mass_density", flag="-mass", default=0.0, reqd=False, 
                about="element mass per unit length"),
        ]
        _refs=["transform", "section"]

        @property
        def section(self):
            return self.integration["section"]

        def init(self):
            pass

    DisplBeamColumn = Ele("DispBeamColumn", "dispBeamColumn",
        #, eleTag, *eleNodes, transfTag, integrationTag, '-cMass', '-mass', mass=0.0)
        about="Create a dispBeamColumn element.",
        args=[
            Tag(),
            Grp("nodes", args=[
              Ref("iNode", type=Node,  attr="name", about=""),
              Ref("jNode", type=Node,  attr="name", about=""),
            ]),
            Ref("geom",  field="transform", type=Trf, attr="name"),
            Ref("integration", type="", alt="quadrature"),
            Flg("-cMass", field="consistent_mass",
                about="Flag indicating whether to use consistent mass matrix."),
            Num("mass",field="mass_density", flag="-mass", default=0.0, reqd=False, 
                about="element mass per unit length"),
        ],
        refs=["transform"],
    )


class constraint:
    RigidBeamLink = Lnk("RigidBeamLink",
        "beam", [Tag(), Grp("nodes", type=Ref(type=uniaxial), num=2)]
    )

    RigidDiaphragm = cmd("rigidDiaphragm", args=[
        Int("perpDirn", about="direction perpendicular to the rigid plane (i.e. direction 3 corresponds to the 1-2 plane)"),
        Ref("rNodeTag", about="integer tag identifying the retained node"),
        Grp("nodes", type=Ref(type="node"), min=2, about="integer tags identifying the constrained nodes"),
    ])



Backbone = LibCmd("backbone")

class backbone:

    Popovics = Mander = Backbone("Mander",
        args = [Tag(), Num("fc"), Num("epsc"), Yng()]
    )

    ReeseSoftClay = Backbone("ReeseSoftClay",
        args=[
          Tag(), Num("pu"), Num("y50"), Num("n")
        ]
    )

    ReeseSand = Backbone("ReeseSand",
        args=[
          Tag(), Num("kx"), Num("ym"), Num("pm"), Num("yu"), Num("pu")
        ]
    )

    ReeseStiffClayBelowWS = Backbone("ReeseStiffClayBelowWS",
        args=[
          Tag(), Num("Esi"), Num("y50"), Num("As"), Num("Pc")
        ]
    )

    Raynor = Backbone("Raynor",
        args=[
          Tag(), Yng(), Yld("stress"), Num("fsu"), Num("Epsilonsh"), Num("Epsilonsm"), Num("C1"), Num("Ey")
        ]
    )

    Backbone("Capped",
        args=[Tag(), Ref("backbone"), Num("capTag")]
    )

    LinearCapped = Backbone("LinearCapped", 
        args=[
          Tag(), Ref("backbone",type=Backbone), Num("eCap"), Num("E"), Num("sRes")
        ]
    )

    Multilinear = Backbone("Multilinear", "Multilinear", args=[
            Grp("points", type=Grp, args=[Grp(type=Num, args=[Num("e"), Num("s")])],
                about="Points defining multilinear backbone."
            ),
        ]
    )


# region = Cmd("region", about="""
#     The region command is used to label a group of nodes and elements. This command
#     is also used to assign rayleigh damping parameters to the nodes and elements in
#     this region. The region is specified by either elements or nodes, not both. If
#     elements are defined, the region includes these elements and the all connected
#     nodes, unless the `-eleOnly` option is used in which case only elements are
#     included. If nodes are specified, the region includes these nodes and all
#     elements of which all nodes are prescribed to be in the region, unless the
#     `-NodeOnly` option is used in which case only the nodes are included. 
#     """,
#     #region $regTag <-ele ($ele1 $ele2 ...)> 
#     #               <-eleOnly ($ele1 $ele2 ...)> 
#     #               <-eleRange $startEle $endEle> 
#     #               <-eleOnlyRange $startEle $endEle> 
#     #
#     #               <-node ($Node1 $Node2...)> 
#     #               <-nodeOnly ($Node1 $Node2 ...)> 
#     #               <-nodeRange $startNode $endNode>
#     #               <-nodeOnlyRange $startNode $endNode> 
#     #               <-node all> 
#     #               <-rayleigh $alphaM $betaK $betaKinit $betaKcomm>
# 
#     # region 1 -ele 1 5 -eleRange 10 15
#     # region 2 -Node 2 4 6 -NodeRange 9 12 
#     args=[
#     Tag(),
#     Grp(flag="-ele",
#         about="tags of selected elements in domain to be included in region (optional, default: omitted)"),
#     Grp(flag="-eleRange",
#         args=[Ref("startEle"), Ref("endEle")],
#         about="tag for start and end elements -- range of selected elements in domain"),
#     Grp(flag="-Node",
#         optional=True,
#         about="tags of selected nodes in domain to be included in region (optional, default: omitted)"),
#     Grp(flag="-NodeRange",
#         optional=True,
#         args=[Ref("startNode"), Ref("endNode")], 
#         about="tag for start and end nodes -- range of nodes in domain"),
#     Grp(flag="-rayleigh", 
#         optional=True,
#         args=["alphaM", "betaK", "betaKinit", "betaKcomm"],
#         about="Arguments to define Rayleigh damping matrix (optional, default: zero)"),
#     ])


redirect = Cmd("redirect",[
      Blk("commands"),
      One(optional=True, enum=[
         Str("filename", flag=">"),
         Str("filename", flag=">>"),
         Grp("shell_pipe", flag="|", args=[
            Str("shell_name", optional=True),
            Blk("shell_block", flag="|"),
         ]),
         Str("variable")
      ])
    ])
