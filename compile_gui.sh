#!/bin/bash
export TMPDIR=$HOME/tmp
mkdir ~/tmp
chmod u+rwx -R ~/tmp
rm -rf build
mkdir build
cd build
cmake .. -D NICE=YES -D USE_LIBTIFF=ON -D GUI=OFF
make -j install
