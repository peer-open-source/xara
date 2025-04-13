from .obj import LibCmd
from .ast import Num, Tag, Str
from opensees.library import Yld, Yng

def BulkModulus():
    return Num("K", about="Bulk modulus")

def ShearModulus():
    return Num("G", about="Shear modulus")

def SaturatedYieldStress():
    return Num("ssat", about="Saturated yield stress.")


_MaterialLibrary = LibCmd("nDMaterial", subs={
    "J2Plasticity": [
        Tag(),
        BulkModulus(),
        ShearModulus(),
        Yld(),
        SaturatedYieldStress(),
        Num("delta", about="exponential hardening parameter"),
        Num("H", about="linear hardening parameter")
    ]
})


def Material(type: str, args, name=None):
    if not hasattr(_MaterialLibrary, type):
        subclass = _MaterialLibrary.subclass(type,
                        args = [Tag()]+[
                            Str(f"arg{i+1}") for i in range(len(args))
                        ]
                 )
    else:
        subclass = getattr(_MaterialLibrary, type)

    return subclass(*args, name=name)

def MultiaxialMaterial(type, parameters, properties=None):
    return getattr(_MaterialLibrary, type)(*parameters)


