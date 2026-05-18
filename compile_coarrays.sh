#!/bin/bash
module load mpich opencoarrays
rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=debug -D USE_COARRAYS=on
make -j install
#exit
# run example
#cafrun -n 2 single_exec prg=tseries_motion_correct nparts=1 nthr=14
