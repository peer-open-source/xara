#!/usr/bin/env python
import os
import sys
from glob import glob
from pathlib import Path
from os.path import basename, splitext

try:
    import amoeba
except:
    amoeba = None
import setuptools

#--------------------------------------------------

version    = "0.0.2"
build_type = "no-build"

#--------------------------------------------------

options = {
        "release": [
            "-DCMAKE_BUILD_TYPE=RELEASE",
            "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION:BOOL=TRUE",
        ],
        # 
        "local": [
            "-DCMAKE_BUILD_TYPE=RELEASE",
            "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION:BOOL=FALSE",
#           "-DProfileBuild:BOOL=TRUE",
        ],
        "debug": [
            "-DCMAKE_BUILD_TYPE=DEBUG",
            "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION:BOOL=FALSE",
        ],
        "no-build": []
}

if os.name == "nt":
    EnvArgs = []

elif "CONDA_PREFIX" in os.environ:
    # Ensure that conda libraries and compilers are used
    EnvArgs = [
        "-DDependencies=Conda",
        f"-DCMAKE_PREFIX_PATH:FILEPATH={os.environ['CONDA_PREFIX']}",
        f"-DCMAKE_IGNORE_PATH:FILEPATH=/usr/lib/;/lib",
    ]

else:
    EnvArgs = [
        "-DDependencies=Unix",
    ]


try:
    assert os.name != "nt"
    import pybind11
    OpenSeesPyRT_Target = ["--target", "OpenSeesPyRT"]
    OpenSeesPyRT_Config = [
        f"-Dpybind11_DIR:FILEPATH={pybind11.get_cmake_dir()}",
        f"-DPYTHON_EXECUTABLE:FILEPATH={sys.executable}"
    ]

except (AssertionError,ImportError):
    OpenSeesPyRT_Config = ["-DNoOpenSeesPyRT=True"]
    OpenSeesPyRT_Target = []


if __name__ == "__main__":
    setuptools.setup(
       data_files=[('bin', [*map(str,Path("win32/").glob("*.*"))]),
       ] if os.name == "nt" else [],
       cmdclass = {"build_ext": amoeba.BuildExtension,
                   "cmake": amoeba.CMakeCommand} if build_type != "no-build" else {},
       ext_modules = [
           amoeba.CMakeExtension(
               name = build_type,
               install_prefix="opensees",
               cmake_build_type=options[build_type][0].split("=")[-1],
               cmake_configure_options = [
#                  "-G", "Unix Makefiles",
                   *EnvArgs,
                   *options[build_type],
                   f"-DOPENSEESRT_VERSION={version}",
                   *OpenSeesPyRT_Config,

               ],
               cmake_build_options=["-j15",
                   "--target", "OpenSeesRT",
                   *OpenSeesPyRT_Target
               ]
           )
       ] if build_type != "no-build" else []
    )

