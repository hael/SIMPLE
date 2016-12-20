#!/bin/bash
set ARGV = `basename -a $1 $2 $3 $4 $5 $6 $7`
tput bold;
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
echo "!                                                                       !"
echo "!                   Make sure that you have sourced                     !"
echo "!               eman2.bashrc before running this srcipt                 !"
echo "!                                                                       !"
echo "! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !"
tput sgr0;

set rundir_in = $1         #/opt/Devel_tools/CTFfind/
set file_framename_in = $2 #target_als_rot.txt
set pixel_in = $3          #1.77
set AccVlt_in = $4         #300 (KeV) acceleration volatage
set sph_abe_in = $5        #sphererical aberration 2.7 (default)
set amp_con_in = $6        #Amplitude constrast default 0.07
set sze_pwr_spc_in = $7    #Size of prower spectrum default 512

tput bold;
tput setaf 1; printf "Input coming into unblur_alignSum.csh:\n"; tput sgr0;
              printf "Unblur exec dir: "; tput bold;
tput setaf 2; printf "$rundir_in\n"; tput sgr0;
              printf "Frame filename: "; tput bold;
tput setaf 3; printf "$file_framename_in \n"; tput sgr0;
              printf "Pixel x: "; tput bold;
tput setaf 6; printf "$pixel_in\n"; tput sgr0;
              printf "Acceleration Voltage (KeV): "; tput bold;
tput setaf 6; printf "$AccVlt_in\n"; tput sgr0;

set rundir = $rundir_in
set exeName = $rundir"ctffind"

#/media/frederic/LaCie/Data/Simple/Titan_Krios/run2-ribos_3262/temp/
set local = ./

set ctf_dir = "CTFfind"
set ctf_txt_dir = "CTFfind_txt"

mkdir $local$ctf_dir
mkdir $local$ctf_txt_dir

set ctffind_dir = $local$ctf_dir
set ctf_txt_dir = $local$ctf_txt_dir

cd $local
pwd

### this is already taken care in the StripData.pl script
set file_framename = $file_framename_in #"target_als_rot.txt"
#ls $local/Als_rot/FoilHole_*frames_als_rot.mrc > $local$file_framename
###

set FILENAME = $local$file_framename #"target_als_rot.txt"

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
set ext_rot = "_rot"
set ext_avr = "_avrot"
set ext_ctf = "_ctffd" #"_pspec" previously called

set file_in = $local$raw_file$ext_mrc          ; echo $file_in
set file_out = $local$raw_file$ext_ctf$ext_mrc ; echo $file_out

set Pix_Sze = $pixel_in           #"1.77"
set Acc_Vlt = $AccVlt_in          #"300"
set Sph_Abe = $sph_abe_in         #sphererical aberration 2.7 (default) "2.7"
set Amp_Con = $amp_con_in         #Amplitude constrast default 0.07 "0.07"
set Sze_Pwr_Spc = $sze_pwr_spc_in #Size of prower spectrum default 512 "512"
set Min_Res = "50.0"
set Max_Res = "3.0"
set Min_Dfc = "5000.0"
set Max_Dfc = "50000.0"
set Dfc_Srh = "500.0"
set Exp_Atg = "100.0"
set Add_Phs = "no"

set file_out_ctf_avr = $local$raw_file$ext_ctf$ext_avr$ext_txt
set file_out_ctf_txt = $local$raw_file$ext_ctf$ext_txt
set file_out_ctf_mrc = $local$raw_file$ext_ctf$ext_mrc

echo $file_out_ctf_mrc

echo `date`
#mprun -Mf hostfile -p myr -np $nproc $exeName << ....END
nice +10 $exeName << ....END
$file_in
$file_out
$Pix_Sze
$Acc_Vlt
$Sph_Abe
$Amp_Con
$Sze_Pwr_Spc
$Min_Res
$Max_Res
$Min_Dfc
$Max_Dfc
$Dfc_Srh
$Exp_Atg
$Add_Phs
....END

mv $file_out_ctf_mrc $ctffind_dir
mv $file_out_ctf_avr $file_out_ctf_txt $ctf_txt_dir

end
end 

echo -e "\nTotal $count Lines read"
echo `date`
