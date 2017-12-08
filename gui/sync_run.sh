#!/bin/bash

LOOPTIME=120
WAITTIME=60

while true; do
	echo "	Finding files ..."
	filesarray=(`find $1 -iname "*$3*"`)
	echo "	Found ${#filesarray[@]} files. Waiting for writes to finish"
	sleep $WAITTIME
	for file in "${filesarray[@]}"; do 
		FILEBASE=`basename $file`
		echo "  " `date` " : " $FILEBASE

		if [ ! -f $2/movies/$FILEBASE ]; then
			if [ "$4" = TRUE ]; then
				mv $file $2/movies/
			else
				cp $file $2/movies/
			fi
		fi
	done
	sleep $LOOPTIME
done

