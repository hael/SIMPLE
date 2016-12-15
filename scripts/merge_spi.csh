#!/bin/bash

echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!                   Make sure that you have added the                   !"
echo "!               sinmple library to your path for this script to work    !"
echo "!               or that your are pointing to the correct path           !"
echo "!                                                                       !"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"

set sptr    = "/"            # set the name for the separator.
set all     = "*"            # unix commands for all files eg: file*

set rundir = /home/frederic/Monash/SourceCode/Simple/Restructured/HansVersion/Simple_Restruct.projet/bin/
set exeName = $rundir"simple_stackops"

#/media/frederic/LaCie/Data/Simple/Titan_Krios/run2-ribos_3262/temp/
set local = ./

set ptcl_dir = "Ribos_eman2/particles"
set ptcl_spi_dir = $local"MrgStck_spi"
set ptcl_spi_file = $ptcl_spi_dir$sptr"MrgStck_spi.spi"

mkdir $local$ptcl_spi_dir

set ctffind_dir = $local$ptcl_dir
set ptcl_spi_dir = $local$ptcl_spi_dir

cd $local
pwd

set file_ptclspiname = "target_merge_spi.txt"
ls $local$ptcl_dir/FoilHole_*frames_als_rot_ptcls.spi > $local$file_ptclspiname

set FILENAME = $local$file_ptclspiname


foreach LINE (`cat $FILENAME`)
@ count ++
tput setaf 2;
echo "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
tput bold;
echo $LINE
tput setaf 4; echo "Processing image(spi) number: $count"
              printf "The merged spi file: ";tput setaf 1;
              printf "$ptcl_spi_file\n"; tput sgr0;
tput setaf 2;
echo "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
tput sgr0

#initialising the copying procedure
if ( $count == 1 ) then
    tput bold; tput setaf 2;
    echo "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~"
    tput setaf 7;
    printf "Initialising the top of the stack\n"   ; tput setaf 4;
    printf "Merging the first file by copying: "   ; tput setaf 2;
    printf "$LINE\n"                               ; tput setaf 4;
    printf "into the list on top of the stack: "   ; tput setaf 1;
    printf "$ptcl_spi_file\n"                      ; tput setaf 7;
    printf "To initialise the start of the stack\n"; tput sgr0;
    tput bold; tput setaf 2;
    echo "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~"
    tput sgr0;
    cp $LINE $ptcl_spi_file
endif

#excuting simple_stackops from the simple library
set execute_simple_stackops = "simple_stackops stk=$LINE stk2=$ptcl_spi_file merge=yes"
tput bold; tput setaf 7; printf "Merging files...\n" ; tput sgr0;

$execute_simple_stackops

end

#Finally getting the info of the Merged stack using simple_iminfo 
tput bold;
tput setaf 2;
echo "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~"
tput setaf 7;
printf "Finally getting the stack info\n"          ; tput setaf 4;
printf "from the stack list using simple_iminfo: " ; tput setaf 1;
printf "$ptcl_spi_file\n"                          ; tput setaf 7;
printf "To initialise the start of the stack\n"    ; tput sgr0;
set execute_simple_iminfo = "simple_iminfo fname=$ptcl_spi_file"
tput bold;
tput setaf 2;
echo "~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~"
tput sgr0;

$execute_simple_iminfo

printf "Total $count Lines read"
echo `date`
