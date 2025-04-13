library="element"

shopt -s globstar
#files=`echo SRC/material/uniaxial/**/*.cpp | sed 's/ /\n/g'`

for i in SRC/element/**/*.cpp; do
  case $i in
    Tcl*)
      continue
      ;;
    *)
      filename="${i##*$library/}"
      classname="${filename##*/}"
      classname="${classname%.*}"
      echo ${filename/$classname/IO_$classname} $classname
      mkdir -p "SRC/api/tclCommandPackage/libraries/$library/${filename/$classname/}"
      {
      printf "\n\n#include <g3_api.h>\n"
      printf "\n\n#include <${i/cpp/h}>\n"
      cat BuildTools/scripts/cut_function.sed \
        | grep -v '^#' \
        | while read -r line; do
        #cmd="$(sed -n "${l}p" )"
            cat $i | clang-format | sed -n "$line"
          done
      } > SRC/api/tclCommandPackage/libraries/$library/${filename/$classname/IO_$classname}
      ;;
    esac
done

