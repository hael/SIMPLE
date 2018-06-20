#!/bin/bash
# Michael eager 2016

set +e #disable fail on error
SIMPLE_ROOT_PATH=/media/seagate1/Projects/Simple-beta
SIMPLE_SRC_PATH=${SIMPLE_ROOT_PATH}/src
SIMPLE_BUILD_PATH=${SIMPLE_ROOT_PATH}/build/lib/simple
print_connection () {
    local inmod=$1
    local modfile=''
    local -a deps
    if [[ $# -ne 2 ]];then
        modfile=$(find $SIMPLE_SRC_PATH $SIMPLE_BUILD_PATH -name "$inmod.f90");
    else
        modfile=$(find $SIMPLE_ROOT_PATH/production  $SIMPLE_SRC_PATH $SIMPLE_BUILD_PATH -name "$inmod.f90");
    fi
    if [ -f $modfile ]; then
#       echo " Could not find $modfile";
#       return 0;
#    fi

        if grep -q "include 'simple_lib.f08'" ${modfile}; then
            echo "\"$inmod\" -> \"simple_lib.f08\"";
        fi
#    echo "printing $modfile";
    deps=($(awk '/^[ ]*use [^i]/ {print $2}'  ${modfile} | tr -d ",'"| tr '\n' ' '))
#    echo "printing $modfile";
#    echo $deps;
    if [ $? -eq 0 ] && [ ${#deps[@]} -ne 0 ]; then
        for name in ${deps[@]}; do
            echo "\"$inmod\" -> \"$name\"";
            print_connection ${name} | sort -u
        done
    fi
    fi
}
print_simple_lib() {
    local -a deps
    local libfile=$SIMPLE_SRC_PATH/inc/simple_lib.f08
    deps=($(awk '/^[ ]*use [^i]/ {print $2}'  $libfile | tr -d ",'") )
    for name in ${deps[@]}; do
        echo "\"simple_lib.f08\" -> \"$name\"";
        print_connection ${name}
    done
}


modtop=$1
outfile="${modtop/\//_}_graph.dot"

modtopfile=$(find $SIMPLE_ROOT_PATH/production  $SIMPLE_SRC_PATH $SIMPLE_BUILD_PATH  -name "$1.f90")
if [ ! -f $modtopfile ]; then
    echo "no module file"; exit 0
else


(
echo "digraph g{
  rankdir=LR;
"
print_simple_lib | sort -u
print_connection $modtop 1
echo "
}
"
) # > $outfile
fi

# if [ -x dot ]; then
#     dot -Tpng $outfile -o ${outfile/dot/png}
#     display ${outfile/dot/png}
# else
#     cat $outfile
#     echo "Go to webgraphviz.com or copy the text above and run: dot -Tpng $outfile -o ${outfile/dot/png} "
# fi
