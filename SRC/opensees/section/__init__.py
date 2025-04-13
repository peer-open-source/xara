from .section import (
      patch,
      FiberSection,
      ElasticSection
)
layer = patch.layer

from .shapes import (
      ConfinedPolygon,
      ConfiningPolygon,
      PolygonRing,
      sect2gmsh
)

SectionGeometry = FiberSection
PlaneShape      = FiberSection

def _WideFlange(aisc_data, mesh_data, material, tag=None, ndm=None)->SectionGeometry:

    if isinstance(aisc_data, dict):
        d  = aisc_data['d']
        bf = aisc_data['bf']
        tf = aisc_data['tf']
        tw = aisc_data['tw']

    if isinstance(mesh_data, dict):
        nft = mesh_data['nft']
        nwl = mesh_data['nwl']
        nfl = mesh_data.get('nfl', 1, ) #mesh_data['nft'])
        nwt = mesh_data.get('nwt', 1, ) #mesh_data['nwl'])

        if ndm is None:
            ndm = mesh_data.get("ndm", 3)

        int_typ = mesh_data.get("IntTyp", None)
        flg_opt = mesh_data.get('FlgOpt', True)
#       web_opt = mesh_data.get('WebOpt', False)

    else:
        assert isinstance(mesh_data, tuple)
        nft, nwl = mesh_data
        nfl, nwt = 1, 1
        int_typ  = None
        flg_opt  = True

    yoff = ( d - tf) / 2
    zoff = (bf + tw) / 4

    dw = d - 2 * tf
    bi = bf - tw

    if ndm == 2:
        GJ =  1.0

    else:
        J = aisc_data.get("J")

        # return SectionGeometry(name=tag, GJ=GJ, areas=[
        #     patch.rect(vertices=[[-],[]], material=material)
        #     # patch.Fiber([x,y], A, material) for x,y,A in zip(xfib, yfib, wfib)
        # ])

    if flg_opt:
        return SectionGeometry(name=tag, GJ=GJ, shapes=[
            patch.rect(corners=[[-bf/2, yoff-tf/2],[bf/2,  yoff+tf/2]], material=material, divs=(nfl, nft), rule=int_typ),
            patch.rect(corners=[[-tw/2,-yoff+tf/2],[tw/2,  yoff-tf/2]], material=material, divs=(nwt, nwl), rule=int_typ),
            patch.rect(corners=[[-bf/2,-yoff-tf/2],[bf/2, -yoff+tf/2]], material=material, divs=(nfl, nft), rule=int_typ),
        ])

    else:
        return SectionGeometry(name=tag, GJ=GJ, shapes=[
            patch.rect(corners=[[-zoff-bi/4, yoff-tf/2],[-zoff+bi/4,  yoff+tf/2]], material=material, divs=(nfl, nft), rule=int_typ),
            patch.rect(corners=[[ zoff-bi/4, yoff-tf/2],[ zoff+bi/4,  yoff+tf/2]], material=material, divs=(nfl, nft), rule=int_typ),
            patch.rect(corners=[[     -tw/2,-yoff-tf/2],[      tw/2,  yoff+tf/2]], material=material, divs=(nwt, nwl), rule=int_typ),
            patch.rect(corners=[[-zoff-bi/4,-yoff-tf/2],[-zoff+bi/4, -yoff+tf/2]], material=material, divs=(nfl, nft), rule=int_typ),
            patch.rect(corners=[[ zoff-bi/4,-yoff-tf/2],[ zoff+bi/4, -yoff+tf/2]], material=material, divs=(nfl, nft), rule=int_typ),
        ])


def from_shape(type, identifier: str, material=None, mesh=None, units=None, ndm=None, tag=None, **kwds):
    if identifier == "WF":
        # TODO
        return _WideFlange

    else:
        return from_aisc(type, identifier, material, mesh, units, ndm, tag, **kwds)


def from_aisc(type, identifier: str, material = None, mesh:dict=None, units=None, ndm=None, tag=None, **kwds):
    if mesh is None:
        mesh = {}
    if units is None:
        import opensees.units.english as units

    aisc_data = load_aisc(identifier, units=units)
    if aisc_data is None:
        raise ValueError(f"Cannot find section with identifier {identifier}")


    if identifier[0] == "W":
        geom = _WideFlange(aisc_data, mesh, tag=tag, ndm=ndm, material=material)

    if type == "Fiber":
        return geom
    if type == "Resultant":
        pass
    if type == "Elastic":
        pass

def load_aisc(SectionName, props="", units=None)->dict:
    """Load cross section properties from AISC database.

    props:
        A list of AISC properties, or one of the following:
        - 'simple': `A`, `Ix`, `Zx`

    """

    from .aisc_imperial import imperial
    SectData = imperial[SectionName.upper()]

    if props == "simple":
        props = ""
        return

    elif props:
        props = props.replace(" ", "").split(",")
        sectData = {k: v for k, v in SectData.items() if k in props}
        if "I" in props:
            sectData.update({"I": SectData["Ix"]})
        return sectData

    for k,v in list(SectData.items()):
        try:
            SectData[k] = float(v)
        except:
            continue

    if units is not None:
        SectData["d"]  *= units.inch
        SectData["bf"] *= units.inch
        SectData["tw"] *= units.inch
        SectData["tf"] *= units.inch
        SectData["A"]  *= units.inch**2
        SectData["Ix"] *= units.inch**4
        SectData["Iy"] *= units.inch**4

    return SectData
