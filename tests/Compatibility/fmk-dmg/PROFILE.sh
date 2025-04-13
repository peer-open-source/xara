#!/usr/bin/bash
#
build="/home/claudio/opensees/OpenSeesRT/build"
lib="/src/libg3/SRC/runtime/libOpenSeesRT.so"
opensees="python -m opensees"
#opensees="$build/temp.linux-x86_64-cpython-39_stack/src/libg3/SRC/executable/OpenSees"
perf="perf stat -B -r 1 -e cache-references,cache-misses,cycles,instructions,faults,page-faults,branches,branch-misses"

export OPENBLAS_NUM_THREADS=12
#script="20Story_LignoImplicit.tcl"
script="fmk2.tcl"
export OPENSEESRT_LIB="$build/temp.linux-x86_64-cpython-39_local/$lib"
#unset OPENSEESRT_LIB
$perf $opensees $script

export OPENSEESRT_LIB="$build/temp.linux-x86_64-cpython-39_stack/$lib"
$perf $opensees $script

