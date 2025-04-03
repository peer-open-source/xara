#!/usr/bin/env python
import os
import sys
from glob import glob
from pathlib import Path
import subprocess
from os.path import basename, splitext

import amoeba
import setuptools

#--------------------------------------------------

version    = "0.1.13"
build_type = os.environ.get("OPENSEESRT_BUILD", "local")

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
        "light": [
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

use_conan = False

if os.name == "nt":
    EnvArgs = []
    use_conan = True

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
    assert False
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


if use_conan:
    EnvArgs + ["-DCMAKE_TOOLCHAIN_FILE=conan/conan_toolchain.cmake"]


class BuildOpenSeesRT(amoeba.BuildExtension):
    def build_extension(self, ext):
        # Ensure Conan dependencies are installed using Conan 2.0 commands
        if use_conan:
            self.run_conan(ext)

        super().build_extension(ext)

    def run_conan(self, ext):
        toolchain = str(Path(".").absolute()/"build"/"generators"/"conan_toolchain.cmake")
        ext.cmake_configure_options.append(
                f"-DCMAKE_TOOLCHAIN_FILE={toolchain}"
        )

        # Create the Conan profile and run the Conan install command
        subprocess.run([
            "conan", "install", ".",
            "--build=missing",
        ])


# BuildOpenSeesRT = amoeba.BuildExtension

if __name__ == "__main__":
    setuptools.setup(
       data_files=[
           ('bin', [*map(str,Path("win32/").glob("*.*"))]),
       ] if os.name == "nt" else [],
       cmdclass = {
            "build_ext": BuildOpenSeesRT, # amoeba.BuildExtension,
            "cmake": amoeba.CMakeCommand
       } if build_type != "no-build" else {},
       ext_modules = [
           amoeba.CMakeExtension(
               name = build_type,
               install_prefix="opensees",
               cmake_build_type=options[build_type][0].split("=")[-1],
               cmake_configure_options = [
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

