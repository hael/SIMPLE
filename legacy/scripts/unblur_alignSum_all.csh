#!/bin/bash
set ARGV = `basename -a $1 $2 $3 $4`
tput bold;
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!                   Make sure that you have sourced                     !"
echo "!               eman2.bashrc before running this srcipt                 !"
echo "!                                                                       !"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
tput sgr0;

set rundir_in = $1           #/opt/Devel_tools/unblur_1.0/bin/
set file_framename_in = $2   #target_frames.txt
set nframes_in = $3          #7
set pixel_in = $4            #1.77

tput bold;
tput setaf 1; printf "Input coming into unblur_alignSum.csh:\n"; tput sgr0;
              printf "Unblur exec dir: "; tput bold;
tput setaf 2; printf "$rundir_in\n"; tput sgr0;
              printf "Frame filename: "; tput bold;
tput setaf 3; printf "$file_framename_in \n"; tput sgr0;
              printf "Number of frames: "; tput bold;
tput setaf 4; printf "$nframes_in \n"; tput sgr0;
              printf "Pixel x: "; tput bold;
tput setaf 6; printf "$pixel_in\n"; tput sgr0;

set rundir = $rundir_in
set exeName = $rundir"unblur"

#/media/frederic/LaCie/Data/Simple/Titan_Krios/run2-ribos_3262/temp/
set local = ./

set als_dir = "Als_rot"
set shf_dir = "Shf_txt"
set rotated_frames_als_dir = "Rotated_Frames_Als"

mkdir $local$als_dir
mkdir $local$shf_dir
mkdir $local$rotated_frames_als_dir

set align_rot_dir = $local$als_dir
set shift_txt_dir = $local$shf_dir
set rotated_align_dir = $local$rotated_frames_als_dir

cd $local
pwd

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
echo $word:q

#"FoilHole_1422566_Data_1393858_1393859_20150502_2204_frames"
set raw_file = $word:q

set ext_mrc = ".mrc"
set ext_txt = ".txt"
set ext_als = "_als"
set ext_shi = "_shf"
set ext_rot = "_rot"

set ext_als_mrc = $ext_als$ext_mrc

set file_in = $local$raw_file$ext_mrc
set nframes = $nframes_in #"7"
set file_out = $local$raw_file$ext_als_mrc
set file_out_shift = $local$raw_file$ext_shi$ext_txt
set pixel = $pixel_in #"1.77"
set dose = "NO"
set expert = "NO"

set file_out_rot = $local$raw_file$ext_als$ext_rot$ext_mrc

echo `date`
#mprun -Mf hostfile -p myr -np $nproc $exeName << ....END
nice +10 $exeName << ....END
$file_in
$nframes
$file_out
$file_out_shift
$pixel
$dose
$expert
....END

e2proc2d.py $file_out $file_out_rot --rotate=90
#rm $file_out
mv $file_out $rotated_align_dir

mv $file_out_rot $align_rot_dir
mv $file_out_shift $shift_txt_dir

end
end 

#ls $shift_txt_dir/

echo -e "\nTotal $count Lines read"
echo `date`
