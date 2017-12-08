#!/bin/bash

REMOVE=FALSE
while getopts s:d:i:r opts; do
   case ${opts} in
      s) SOURCE=${OPTARG} ;;
      d) DESTINATION=${OPTARG} ;;
      i) IDENTIFIER=${OPTARG} ;;
      r) REMOVE=TRUE ;;
      :) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
   esac
done

if [ ! -n "$SOURCE" ]; then
	echo "ERROR: Source not specified with -s flag"
	exit 1
fi

if [ ! -n "$DESTINATION" ]; then
	echo "ERROR: Destination not specified with -d flag"
	exit 1
fi

if [ ! -n "$IDENTIFIER" ]; then
	echo "ERROR: Identifier not specified with -i flag"
	exit 1
fi

if [ ! -e "$SOURCE" ]; then
	echo "ERROR: Source does not exist"
	exit 1
fi

if [ ! -e "$DESTINATION" ]; then
	echo "ERROR: Destination does not exist"
	exit 1
fi

DESTINATION="$DESTINATION/`basename $SOURCE`"

mkdir -p  $DESTINATION

if [ ! -d $DESTINATION ]; then
	echo "ERROR: Error creating destination directory"
	exit 1
fi

mkdir -p  $DESTINATION/movies

if [ ! -d $DESTINATION/movies ]; then
	echo "ERROR: Error creating destination movies directory"
	exit 1
fi

mkdir -p  $DESTINATION/xml

if [ ! -d $DESTINATION/xml ]; then
	echo "ERROR: Error creating destination xml directory"
	exit 1
fi

echo ""
echo "Syncing movies :"
echo "	Source             : $SOURCE"
echo "	Destination        : $DESTINATION"
echo "	Movie Identifier   : $IDENTIFIER"

if [ "$REMOVE" = TRUE ]; then
	echo "	Remove Source      : True"
else
	echo "	Remove Source      : False"
fi

echo ""

srun -n 1 -c 1 -p eclipse_krios_sync ./sync_run.sh $SOURCE $DESTINATION $IDENTIFIER $REMOVE



