#!/bin/bash
#

#OPENSEES=OpenSees
OPENSEES="python -m opensees"

Results="README.md"
#Results="Memory.md"
export OPENSEESRT_LIB=/home/claudio/packages/opensees-pypi/build/temp.linux-x86_64-cpython-39_local/src/libg3/SRC/runtime/libOpenSeesRT.so
# export LD_PRELOAD=/home/claudio/mambaforge/envs/py39/x86_64-conda-linux-gnu/lib/libasan.so
# export OPENSEESRT_LIB=/home/claudio/packages/opensees-pypi/build/temp.linux-x86_64-cpython-39_debug/src/libg3/SRC/runtime/libOpenSeesRT.so

cat - > $Results <<EOF
# Examples

| Status  |     File     |
|---------|--------------|
EOF

# -e Shell_Dynamic
#
ls $@ | grep -v -e "Rocking" -e "Brick27"  | while read i; do
  echo "FILE $i"; 
  $OPENSEES $i; 
  printf "| $?\t|"' `'"$i"'` |\n' | tee -a $Results;
done

