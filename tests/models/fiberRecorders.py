import json, sys, fnmatch

import numpy as np
from opensees import patch, section
HELP = """
python -m fiberRecorders model.json output.txt -d Dcol -e E... -l STATE

Options
-file | -xml            output file name
-logging <int>          logging flag
-d <float>              section diameter
-s <int>,[<int>]        section tag range

To install run:
    python fiberRecorders.py --setup develop

"""

# --8<--------------------------------------------------------
def damage_states(Dcol):
    cover = 2.0
    Rcol = Dcol/2
    coverl = cover + Rcol*(1-np.cos(np.pi/8))
    return {
        "dsr1" : {
            "regions": [
                 # external radius    internal radius
                section.PolygonRing(8, Rcol,         Rcol-coverl/4)
            ]
        },
        "dsr2" : {
            "regions": [
                section.PolygonRing(8, Rcol,         Rcol-coverl/4)
            ]
        },
        "dsr3" : {
            "regions": [
                section.PolygonRing(8, Rcol-cover/2, Rcol-3*cover/4)
            ],
            "material": "*concr*"
        },
        "dsr4" : {
            "regions": [
                section.PolygonRing(8, Rcol-3*coverl/4, Rcol-coverl)
            ],
            "material": "*concr*"
        },
        "dsr5": {
            "regions": [
                section.FiberSection(areas=[
                    patch.circ(intRad=Rcol - cover - 2, extRad=Rcol - cover)
                ])
            ],
            "material": "*concr*"
        },
        "dsr6": {
            "regions": [
                section.FiberSection(areas=[
                    patch.circ(intRad=Rcol - cover - 2, extRad=Rcol - cover)
                ])
            ],
            "material": "*steel*"
        },
        "all": {
            "regions": [
                section.ConfinedPolygon(8, Rcol)
            ]
        }
    }

# --8<--------------------------------------------------------

def iter_elem_fibers(model:dict, elements:list, sections: list=(0,-1), filt:dict=None):
    sam = model["StructuralAnalysisModel"]
    model["sections"] = {
        str(s["name"]): s for s in sam["properties"]["sections"]
    }
    model["materials"] = {
        str(m["name"]): m for m in sam["properties"]["uniaxialMaterials"]
    }
    for el in sam["geometry"]["elements"]:
        if el["name"] in elements and "sections" in el:
            #for tag in el["sections"]:
            for i in sections:
                idx = len(el["sections"]) - 1 if i==-1 else i
                tag = el["sections"][idx]
                s = model["sections"][tag]
                if "section" in s:
                    s = model["sections"][s["section"]]
                    for s,f in iter_section_fibers(model, s, filt):
                        yield el,idx+1,f

                        
                        
def iter_elem_strains(model, el, s, recorder_data):
    sam = model["StructuralAnalysisModel"]
    model["sections"] = {
        str(s["name"]): s for s in sam["properties"]["sections"]
    }
    model["materials"] = {
        str(m["name"]): m for m in sam["properties"]["uniaxialMaterials"]
    }
    for el in sam["geometry"]["elements"]:
        if el["name"] in elements and "sections" in el:
            #for tag in el["sections"]:
            for i in sections:
                idx = len(el["sections"]) - 1 if i==-1 else i
                tag = el["sections"][idx]
                s = model["sections"][tag]
                if "section" in s:
                    s = model["sections"][s["section"]]
                    for s,f in iter_section_fibers(model, s, filt):
                        yield el,idx+1,f

    
def iter_section_fibers(model, s, filt=None):
    if filt is not None:
        if "material" not in filt:
            filt["material"] = "*"
        for fiber in s["fibers"]:
            if any(
                fiber["coord"] in region
                for region in filt["regions"]
            ) and fnmatch.fnmatch(
                model["materials"][fiber["material"]]["type"].lower(),
                filt["material"]
            ):
                yield s,fiber
    else:
        for fiber in s["fibers"]:
            yield s,fiber

    
def print_fiber(c, s, base_cmd, options):
    out_file = options["record_file"] +f"_{s}" + "_" + str(c[0]) + "_" + str(c[1]) + ".txt"
    base_cmd = base_cmd.format(fmt=options["format"],out_file=out_file)
    fiber_cmd = base_cmd + f"fiber {c[0]} {c[1]} stressStrain;\n"
    print(fiber_cmd)
    if options["logging"]:
        print("puts \""+fiber_cmd+"\"")

        
def getDictData(allData, curDict):
    "Arpit Nema"
    if type(curDict) == dict:
        for key in curDict.keys():
            if type(curDict[key]) == dict:
                getDictData(allData, curDict[key])
            else:
                if type(curDict[key]) == int:
                    curDict[key] = allData[:, curDict[key]]
                    
def read_sect_xml(xml_file):
    "Arpit Nema"
    import xml.etree.ElementTree as ET
    with open(xml_file) as f:
        xml_data = f.read()
    root = ET.XML(xml_data)  # Parse XML

    hdrs = []
    dataDict = {}
    colCtr = 0
    for i, child in enumerate(root):
        if child.tag == "TimeOutput":
            hdrs.append(child[0].text)
            dataDict[child[0].text] = colCtr
            colCtr += 1
        if child.tag == "ElementOutput":
            eleKey = child.attrib["eleTag"]
            secKey = child[0].attrib["number"]
            hdrPre = eleKey + "_" + secKey + "_" + child[0][0].attrib["secTag"]

            dataDict[eleKey] = {secKey: {}}
            for respCtr in range(len(child[0][0])):
                hdrs.append(hdrPre + "_" + child[0][0][respCtr].text)
                respKey = child[0][0][respCtr].text
                if respKey in dataDict[eleKey][secKey].keys():
                    respKey = respKey + "_"
                dataDict[eleKey][secKey][respKey] = colCtr
                colCtr += 1
        if child.tag == "Data":
            tmp = child.text
    data = np.array(tmp.replace("\n", "").split(), dtype=float).reshape((-1, len(hdrs)))
    getDictData(data, dataDict)
    return dataDict

def fiber_strain(recorder_data, el, s, f, t=None):
    if t is not None:
        eps = recorder_data[str(el)][str(s)]["eps"][t]
        kz =  recorder_data[str(el)][str(s)]["kappaZ"][t]
        ky =  recorder_data[str(el)][str(s)]["kappaY"][t]
    else:
        eps = recorder_data[str(el)][str(s)]["eps"]
        kz =  recorder_data[str(el)][str(s)]["kappaZ"]
        ky =  recorder_data[str(el)][str(s)]["kappaY"]

    return eps - kz * f["coord"][1] + ky * f["coord"][0]

def print_help():
    print(HELP)

REQUIREMENTS = """
numpy
opensees
"""

def install_me():
    import os
    import subprocess
    try:
        from setuptools import setup
    except ImportError:
        from distutils.core import setup

    sys.argv = sys.argv[:1] + sys.argv[2:]

    setup(name = "fiberRecorders",
          version = "0.0.1",
          description = "Utilities for creating fiber recorders",
          long_description = HELP,
          author = "",
          author_email = "",
          url = "",
          py_modules = ["fiberRecorders"],
          scripts = ["fiberRecorders.py"],
          license = "",
          install_requires = [*REQUIREMENTS.strip().split("\n")],
    )


def parse_args(args)->dict:
    opts = {
        "model_file": None,
        "record_file": None,
        "logging": 0,
        "format": "file"
    }
    argi = iter(args)
    for arg in argi:
        if arg == "--setup":
            install_me()
            sys.exit()

        if arg in ["--help","-h"]:
            print_help()
            sys.exit()

        if arg[:2] == "-d":
            opts["Dcol"] = float(next(argi))

        elif arg == "-e":
            opts["elements"] =  [
                int(e) for e in next(argi).split(",")
            ]

        elif arg == "-s":
            opts["sections"] = [int(s) for s in next(argi).split(",")]

        elif arg == "-l":
            opts["state"] = next(argi)
        
        # Formats
        elif arg == "-file":
            opts["format"] = "file"
        elif arg == "-xml":
            opts["format"] = "xml"

        elif arg == "-logging":
            opts["logging"] = next(argi)

        elif opts["model_file"] is None:
            opts["model_file"] = arg


        else:
            opts["record_file"] = arg

    return opts

base_cmd = "recorder Element -{fmt} {out_file} -time "


if __name__=="__main__":

    opts = parse_args(sys.argv[1:])

    damage_state = damage_states(opts["Dcol"])[opts["state"]]

    elements = opts["elements"]

    sections = opts["sections"]

    with open(opts["model_file"], "r") as f:
        model = json.load(f)

    for e,s,f in iter_elem_fibers(model, elements, sections, damage_state):
        elem_cmd = base_cmd + f"-ele {e['name']} "
        print_fiber(f["coord"], s, elem_cmd + f"section {s} ", opts)