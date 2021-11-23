#!/bin/bash
PATH=$SIMPLE_PATH/bin/nodejs/bin/:$PATH npm --prefix=$SIMPLE_PATH/gui_data/server/single_user/ start "${@:1}"

