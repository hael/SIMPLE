#!/bin/bash
PATH=$SIMPLE_PATH/bin/nodejs/bin/:$PATH node --max-old-space-size=8192 $SIMPLE_PATH/gui_data/server/multi_user/simple_multiuser.js "${@:1}"

