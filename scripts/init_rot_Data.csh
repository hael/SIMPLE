#!/bin/bash
set ARGV = `basename -a $1 $2 $3 $4 $5`
tput bold;
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!                   Make sure that you have sourced                     !"
echo "!               eman2.bashrc before running this srcipt                 !"
echo "!                                                                       !"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
tput sgr0;

set dir_data_in = $1
set dir_data_rot_in = $2           #Dir where the rotated Data
set file_framename_in = $3   #target_frames.txt
set ext_in = $4              #extension name input variable
set deg_in = $5              #rotation to undertake

tput bold;
tput setaf 1; printf "Input into initial_rotation_Data.csh:\n"; tput sgr0;
              printf "Data_MRC dir: "; tput bold;
tput setaf 2; printf "$dir_data_in\n"; tput sgr0;
              printf "Data Rotated by (deg) dir: "; tput bold;
tput setaf 2; printf "$dir_data_rot_in\n"; tput sgr0;
              printf "Data filename: "; tput bold;
tput setaf 3; printf "$file_framename_in \n"; tput sgr0;
              printf "Extension name for rotated files: "; tput bold;
tput setaf 3; printf "$ext_in \n"; tput sgr0;
              printf "Amount of rotation (anti-clockwise): "; tput bold;
tput setaf 4; printf "$deg_in (deg)\n"; tput sgr0;

set outdir = $dir_data_rot_in
set ext = $ext_in
set deg = $deg_in
set exeName = "e2proc2d.py"

#/media/frederic/LaCie/Data/Simple/Titan_Krios/run2-ribos_3262/temp/
set local = ./

mkdir $outdir

cd $local
pwd

set sptr = "/"
set ext_mrc = ".mrc"
set ext_txt = ".txt"

### this bit is already taken care in StripData.pl 
set file_framename = $file_framename_in
#ls $local/FoilHole_*frames.mrc > $local$file_framename
###
set FILENAME = $local$file_framename

@ count = 0

foreach LINE (`cat $FILENAME`)
@ count ++
tput setaf 2;
echo "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
tput bold;
echo $LINE
tput setaf 4; echo "Processing image(mrc) number: $count"; tput sgr0;
tput setaf 2;
echo "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
tput sgr0

set s = $LINE
set my_raw = `echo $LINE:q | sed 's/.mrc/ /g'`
foreach word($my_raw:q)
#echo $word:q
#"Data_MRC/FoilHole_1422566_Data_1393858_1393859_20150502_2204"
set raw_file = $word:q

set s = $raw_file
set my_new_raw_file = `echo $raw_file:r | sed 's/Data_MRC\// /g'`
foreach newword($my_new_raw_file:r)
#echo $newword:r
#"FoilHole_1422566_Data_1393858_1393859_20150502_2204"
set new_raw_file = $newword:r

set new_file = $new_raw_file$ext$ext_mrc

end
end

echo $outdir$sptr$new_file 
$exeName $LINE $new_file --rotate=$deg
mv $new_file $outdir

end

#ls $outdir

echo -e "\nTotal $count Lines read"
echo `date`
