#!/usr/bin/perl
use lib './';
use warnings;
use strict;
use Cwd qw(getcwd);
use Env;
use Config;
use Tie::File;
use Term::ANSIColor;
use XML::Parser;
use Data::Dumper;
################################################################################
# Perl script to strip the data and recast in proper directory structure       #
#   **** Do not modify this script, this sript automates the compilation ****  #
#   ****                         procedure                               ****  #
#                                                                              #
# Input variables should be modifed in the file:                               #
#                                                                              #
#                          simple_StripData_input.pm                           #
#                                                                              #
################################################################################
use simple_StripData_input;
################################################################################
# Declare global variables                                                     #
################################################################################
my $sptr    = "/";           # set the name for the separator.
my $b_sptr  = "\\";          # set the name for the separator.
my $all     = "*";           # unix commands for all files eg: file*
my $eq      = "=";           # equal sign for file writting
my $my_path = getcwd();      # getting the command path
my $true_path;               # true path
my $execdir;                 # the fixed execdir path
#files extension names
my $ext;                           #generic extension filename
my $ext_txt        = ".txt";       # txt file name extension: file*.txt
my $ext_xml        = ".xml";       # xml file name extension: file*.xml
my $ext_jpg        = ".jpg";       # jpg file name extension: file*.jpg
my $ext_mrc        = ".mrc";       # mrc file name extension: file*.mrc
my $ext_dmf        = ".dm";        # dm file name extension: file*.dm
my $ext_str        = ".star";      # extension name for star files in RELION
my $ext_asc        = ".asc";       #file extension for the *.asc files
my $ext_log        = ".log";       #file extension for the *.log files
my $ext_frames_mrc = "frames.mrc"; # file*_frames.mrc
my $ext_Data_mrc   = "Data";       # Data.mrc filename extension: file*Data*.mrc
my $ext_Trg        = "TargetLocation"; # TargetLocation file name extension
my $ext_Grd        = "GridSquare"; # GridSquare file name extension
my $ext_box        = ".box";       # .box extension file variable
my $ext_ctffd      = "_ctffd";     #"_pspec"; # TODO:change to "_ctffd"
my $ext_def        = "_def";       # def tabs extension name
my $ext_als        = "_als";       # extension name for als (aligned sum)
my $ext_shf        = "_shf";       # shf tabs extension name
my $ext_rot        = "_rot";       # extension for the rotation
#key words in the file
my $ext_frames      = "_frames";     #key word frames in the file
#Directories name
my $dir;                               #generic directory
my $dir_loc             = ".";         # local dir for all of the file*.*
my $dir_xml             = "XML";        # XML dir for all of the xml file*.xml
my $dir_jpg             = "JPG";        # JPG dir for all of the xml file*.jpg
my $dir_frames_mrc      = "Frames_MRC"; # *.frames.mrc dir for file*_frames.mrc
my $dir_Data_mrc        = "Data_MRC";   # *Data*.mrc dir for file*_frames.mrc
my $dir_Data_rot_mrc    = "DRot_MRC";   # *Data*.mrc dir for file*_frames.mrc
my $dir_Othr_mrc        = "Othr_MRC";   # Other *.mrc dir for file*other*.mrc
my $dir_TargetLocation  = "TargetLocation"; # TargetLocations dir
my $dir_GridSquare      = "GridSquare"; # GridSquare
my $dir_Boxfile         = "Boxfile";    # Boxfile directory
my $dir_Als_rot_mrc     = "Als_rot";    #*_als_rot.mrc dir for file*_als_rot.mrc
my $dir_CTFfind         = "CTFfind";    #CTFFind dir
my $dir_CTFfind_txt     = "CTFfind_txt";# CTFfind_txt directory
my $dir_Shf             = "Shf_txt";    # Shf_txt directory
my $dir_Deftabs         = "Deftabs";    # Deftabs directory
my $dir_Rotd_Frames_Als = "Rotated_Frames_Als"; #rotated frames als_rot dir
#XML parsed values
my $x_pixSze;                 #The parsed x PixelSize
my $y_pixSze;                 #The parsed y PixelSize
my $TemMag;                   #Nominal Magnification
my $AccVlt;                   #Acceleration Voltage (KeV)
my $message;
#The organiser value
my $has_been_organised=0; #stamp variable for file organisation default: "un"
my $orga = 0;          #file organisation variable 0:not 1:organised   
#The graphing program
my $graph   = "/usr/bin/graph -T X --bg-color black --frame-color white";
#Arrays to handles files and the like
my @targetlocation;  #array for all of the directory structure
my @nonemptyString;  #the non empty string containing file nane with **
my @nonemptyTarget;  #The non empty TargetLocation array
my @datastatTarget;  #the quality control array                 
my @targetFrames;    #array for the targetFrames from the target_frames.txt
my @targetDataMRC;   #array for the targetDataMRC from the target_Data_MRC.txt
my @xml_file_array;  #xml file array 
my @xml_dir_array;   #xml dir array
my @lines;           #Array for the XML parsing files
my @nparticles;      #array for the nparticles per boxed files
my @boxed_file;      #array of files names of the bnoxed files
my @file_array_ctf;  #array for ctf files 
my @file_array_shf;  #array for Shifted files.

my @file_array_ctf_recov;  #array for ctf files for recovery
my @file_array_shf_recov;  #array for Shifted files  for recovery

my @ctf_values;      #ctf values from the CTFFind_txt files
my @targetShiftVal;  #array for target shifted values
my @shf_values_x;    #shifted values from the Sht_txt files for x
my @shf_values_y;    #shifted values from the Sht_txt files for y
my @nx_array;        #n xi array
my @ny_array;        #n yi array
#Files names for the target folders
my $target_frames_file    = "target_frames";
my $target_data_mrc_file  = "target_Data_MRC";
my $target_organiser_file = "target_organiser";
my $target_als_rot_file   = "target_als_rot";
my $target_Boxfile_file   = "target_Boxfile";
my $target_ShiftVal_file  = "target_ShiftVal";
my $dftb_stk_file         = "Deftabs_stack.asc";
my $dftb_stk_full_file    = "Deftabs_full_stack.asc";
my $dftb_stk_full_nptcls_file    = "Deftabs_full_stack_nptcls.asc";
my $stat_box_file         = "dfxdfy_box_stat.asc";
my $stat_shf_file         = "Shftxy_shf_stat.asc";
my $stat_format_shf_file  = "Shftxy_formated_shf_stat.asc";
my $star_file             = "project.star";
my $simple_movies_file    = "simple_movies.txt";
my $simple_Boxfile_file   = "simple_Boxfile.txt";
#File hadnlers
my $xml_file;             #XML file to be parsed
my $f1;                   #file handler for file generic
my $f2;                   #file handler for file generic
my $fscn;
my $ftar;
my $ffoo;
my $fxml;
my $ffrm;                 #file handler for file target_frames.txt
my $forg;                 #file handler for file organaiser target_organiser
my $fstk;                 #file handler for file Deftabs_stack.asc
my $fstk_full;            #file handler for file Deftabs_stack.asc
my $fstk_full_nptcls;     #file handler for file Deftabs_full_stack_nptcls.asc
my $fbox;                 #file handler for file dfxdfy_box_stat.asc
my $fctf;                 #file handler for file file_ctf
my $fdef;                 #file handler for file file_def
my $fstr;                 #file handler for file file_star
my $fdef_full;            #file handler for file file_def_full
my $fshf;                 #file handler for file target_ShiftVal.txt
my $fshf_txt;             #file handler for file Shf_txt/*.txt
my $fshf_stat;            #file handler for file Shftxy_shf_stat.asc
my $fshf_format_stat;     #file handler for file Shftxy_formated_shf_stat.asc
my $fdta_mrc;             #file handler for file target_Data_MRC.txt
my $fsmv;                 #file handler for file simple_movies.txt
my $fsbx;                  #file handler for file simple_Boxfile.txt
#enumerators and counters
my $nrows = 0;            #number of lines in the file target.txt
my $irows;
my $cnt_str = 0;          #counter for the star file
#input variables for the c shell script unblur and ctffinder4
my $rundir_in;            #/opt/Devel_tools/unblur_1.0/bin/
my $ext_in;               #extension name input for c-shell scripts
my $file_framename_in;    #target_frames.txt
my $nframes_in;           #7
my $pixel_in;             #1.77
my $AccVlt_in;            #300 (KeV) acceleration volatage 
my $sph_abe_in;           #sphererical aberration 2.7 (default)
my $amp_con_in;           #Amplitude constrast default 0.07
my $sze_pwr_spc_in;       #Size of prower spectrum default 512
my $deg_in;               #Roation angle in degrees for the c-shell scripts
my $dir_data_in;          #directory where Data_MRC files lives in
#handlers varaibles for the parsing of the box files and the extraction
my $file_box;
my $file_ctf;
my $file_def;
#handlers variables for the parsing of the Shf_txt files
my $file_shf;
#recivery handlers for parsing recovery
my $file_shf_reco;
my $file_ctf_reco;
#utilities string 
my $die_open_string = "cannot open input file";
# The x_i and y_i for file output from 
my $xi;
my $yi;
my $pixel_in_x;                #1.77
my $pixel_in_y;                #1.77
my $shf_values_x_normalised;
my $shf_values_y_normalised;
################################################################################
# Start of the execution commands                                              #
################################################################################
################################################################################
# First check the integrity of the path and see if there are no spaces in it   #
################################################################################
#$execdir = $my_path;
get_check_fix_path($my_path);
$execdir = $true_path;
################################################################################
# Now we can start if it is all good                                           #
################################################################################
print color('bold blue');
print qq[-------------------------------------------------------------------\n];
print qq[   Perl script to strip and analyse mrc files from Krios-Titan     \n];
print qq[                  By Frederic D.R. Bonnet major dates              \n];
print qq[           11th May. 2015: file distribution and regrouping        \n];
print qq[           08th Aug. 2015: addition of the organisation            \n];
print qq[           09th Aug. 2015: addition of the XML parsing             \n];
print qq[                           intergration of the c-shell script for  \n];
print qq[                           unblur and CTFFind4 calculation         \n];
print qq[           14th Aug. 2015: creation of star files for RELION       \n];
print qq[-------------------------------------------------------------------\n];
print qq[Our working path is                         : $SIMPLE_DATA_PATH\n];
print qq[The execution dir                           : $execdir\n];
print qq[Number of frames file_frames.mrc file       : $NFRAMES\n];
print qq[mrc files to be stripped in the dataset     : $NCONFIG\n];
print qq[Sphererical Aberration                      : $SPH_ABE\n];
print qq[Amplitude constrast                         : $AMP_CON\n];
print qq[Size of power spectrum in CTFFind4          : $SZE_PWR_SPC\n];
print qq[Creating directories(distribute) 0:no 1:yes : $CREATE_DIR\n];
print qq[File Handler                                : $FILE_PRE_HANDLER\n];
print qq[File organisation                           : $FILE_ORGANISE_DIR\n];
print qq[Preprocessing (unblur) 0:no 1:yes           : $UNBLUR\n];
print qq[Extracting CTF params (CTFFind4) 0:no 1:yes : $CTFFIND\n];
print qq[Unblur system directory                     : $UNBLUR_DIR\n];
print qq[CTFFind 4 systrem directory                 : $CTFFIND_DIR\n];
print qq[File Post-Preprocessing handler variables   : $FILE_PST_HANDLER\n];
print qq[Boxing variable 0:no 1: yes (in Boxfile dir): $HAVE_BOXFILE\n];
print qq[Deleting dirs (when distributed) 0:no 1:yes : $DELETE_DIR\n];
print qq[-------------------------------------------------------------------\n];

print color('bold white');
if ( $SIMPLE_DATA_PATH ne $execdir ) {
    print qq[\n];
    print qq[-----Our working path is not the same as The execution dir-----\n];
    print qq[---------Check that the paths are the same in input module-----\n];
    print qq[-----------------The script is terminated----------------------\n];
    print qq[\n];
    die;
}
print color('reset');
################################################################################
# Setting up the target file                                                   #
################################################################################
$f1 = $TARGET_FILE;  # setting up the target file name
if ($FILE_PRE_HANDLER =~ /distribute/ || $FILE_PRE_HANDLER =~ /lookupdata/ ) {
    generate_TargetLocation_list($f1);
}
if ($FILE_PRE_HANDLER =~ /preprocess/ ) {
    generate_TargetLocation_list($f1);
}
################################################################################
# Extract the directory names from the target.txt file                         #
################################################################################
get_target_file($f1);
################################################################################
# Scan the directory against the target.txt file                               #
################################################################################
if ( $FILE_PRE_HANDLER ne "preprocess" ) {
    if ( $FILE_PST_HANDLER ne "extract_data" ) {
	scan_directory();
    }
}
################################################################################
# Creating the directory structure for the file system                         #
################################################################################
if ($FILE_PRE_HANDLER =~ /distribute/) {
    if ( $CREATE_DIR == 1 ) {
	create_TargetDirectories();
	#system("ls -sF");
    }
}
################################################################################
# Now organizing the data in their associated locations                        #
################################################################################
if ($FILE_PRE_HANDLER =~ /lookupdata/ ) {
    print qq[Nothing to be done\n];
    print qq[Just the graph stats of the data\n];
    $has_been_organised = 0;  #setting the organiser file 0:default
    set_organiser_file($has_been_organised);
} elsif ($FILE_PRE_HANDLER =~ /distribute/) {
    distribute_files();
    #system("ls -sF");
} elsif ($FILE_PRE_HANDLER =~ /regroup/ ) {
    regroup_files();
    #system("ls -sF");
} elsif ($FILE_PRE_HANDLER =~ /preprocess/) {
    #first arrange the data in a more suitable format with the directories
    if ($FILE_ORGANISE_DIR =~ /un/) {
	print_unorganised_message();
	unorganise();
    } elsif ($FILE_ORGANISE_DIR =~ /or/) {
	print_organised_message();
	#organising the files in manageable directory structure
	organise();
    } elsif ($FILE_ORGANISE_DIR =~ /re/) {
	print_reorganised_message();
	#reorganising the files into its orginal structure
	reorganise();
    } elsif ($FILE_ORGANISE_DIR =~ /no/) {
	$f1 = $target_organiser_file.$ext_txt;
	get_target_oragniser_file($f1);
	if ($orga == 0) {
	    $has_been_organised = 0; #stamp variables to keep track organisation
	} elsif ($orga == 1) {
	    $has_been_organised = 1; #stamp variables to keep track organisation
	}
    }
}
################################################################################
# Now prepropcesing the data with unblur and CTFfinder4                        #
################################################################################
if ($FILE_PRE_HANDLER =~ /preprocess/) {
    print_preprocessing_message();
    #first generate the files list from the Frames folder
    $f1 = $target_frames_file.$ext_txt;
    if ($has_been_organised == 0) {
	generate_Target_frames_list($f1,$dir_loc,$ext_frames_mrc);
    } elsif ($has_been_organised == 1) {
	generate_Target_frames_list($f1,$dir_frames_mrc,$ext_mrc);
    }
    #Next extract the file names from the target_frames.txt file
    get_target_Frames_file($f1);
    #Parsing the xml files for the target_files only need the first 1 or 2 
    parse_XML_file();
    #processing the data using unblur
    if ($UNBLUR == 1) { preprocessing_unblur(); }
    #after having unblur proceeding with CTFFind4 on the unblur images
    $f1 = $target_als_rot_file.$ext_txt;
    if ($has_been_organised == 0) {
	generate_Target_als_rot_list($f1,$dir_loc,$ext_frames_mrc);
    } elsif ($has_been_organised == 1) {
	generate_Target_als_rot_list($f1,$dir_Als_rot_mrc,$ext_mrc);
    }
    if ($CTFFIND == 1) { preprocessing_ctffind(); }

    if    ( $HAVE_BOXFILE == 0 ) {print_message_need_box();}
    elsif ( $HAVE_BOXFILE == 1 ) {print_message_have_box();}
}
################################################################################
# Now Parsing the files and generating input files from the prepropcesing      #
# coming from unblur and CTFfinder4                                            #
################################################################################
if ($FILE_PST_HANDLER =~ /extract_data/) {

    if ( $HAVE_BOXFILE == 1 ) {
	print_message_have_box();
	#creatign the target_Boxfile.txt
	create_boxfile_target();
	#getting the number of particles from the target_Boxfile.txt file
	$f1 = $target_Boxfile_file.$ext_txt;
	get_nboxed_2_defTab_boxfile($f1);
	#creatign the target_ShiftVal.txt
	create_shiftval_target();
	$f1 = $target_ShiftVal_file.$ext_txt;
	get_target_ShiftVal_file($f1);

	#parsing the entire data set and printing the stats from all files
	system("mkdir $dir_Deftabs");
	construct_Deftabs_dfxdfy_asc_file();
    } elsif ( $HAVE_BOXFILE == 0 ) {
	print_message_need_box();
	$f1 = $target_data_mrc_file.$ext_txt;
	if ($has_been_organised == 0) {
	    generate_Target_data_list($f1,$dir_loc,$ext_mrc);
	} elsif ($has_been_organised == 1) {
	    generate_Target_data_list($f1,$dir_Data_mrc,$ext_mrc);
	}
	$f1 = $target_data_mrc_file.$ext_txt;
	get_target_Data_MRC_file($f1);
	rotate_Data_MRC_files();
	print_message_hasbeen_rot_box();
    }

}
################################################################################
# Deleting the directory structure for the file system                         #
################################################################################
if ($FILE_PRE_HANDLER =~ /regroup/) {
    if ( $DELETE_DIR == 1 ) {delete_TargetDirectories();}
}
################################################################################
# Deleting the generated preprocessing files                                   #
################################################################################
if ($FILE_PRE_HANDLER =~ /preprocess/) {
    if ($FILE_ORGANISE_DIR =~ /re/) {
	if ( $DELETE_DIR == 1 ) {delete_PreprocessFiles();}
    }
}
################################################################################
#                           Subroutines                                        #
################################################################################
################################################################################
#subroutine to e2proc2d.py the rotated ./Data_MRC/*.mrc to match the correct   #
# _frames.mrc orientation. Invoking init_rot_Data.csh c-shell script           #
#Frederic D.R. Bonnet, Date: 17th of Aug. 2015.                                #
################################################################################
sub rotate_Data_MRC_files {

    $dir_data_in = $dir_Data_mrc;            #Data in where the Data_MRC lives
    $rundir_in = $dir_Data_rot_mrc;          #Dir where the rotated Data
    $file_framename_in = $target_data_mrc_file.$ext_txt;
    $ext_in = $ext_frames.$ext_als.$ext_rot;
    $deg_in = $DATA_ROT;
    print_input_var_initial_rotation();
    #Here rotation script to rotate data files for a given angle to match
    #the orientation of the *_frames.mrc files.
    #It is these files that need to be boxed.
    system("csh init_rot_Data.csh $dir_data_in $rundir_in $file_framename_in $ext_in $deg_in");

    print color('bold red');
    print qq[init_rot_Data.csh Done...\n];
    print color('reset');

    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to create the \$fstk=Deftabs_stack.asc                             #
#                         \$fbox=dfxdfy_box_stat.asc                           #
#                         \$fdef=./Deftabs/*.txt                               #
# files for the shifted values simple input                                    #
#Frederic D.R. Bonnet, Date: 10th of Aug. 2015.                                #
################################################################################
sub construct_Deftabs_dfxdfy_asc_file {

    my $file_movie_txt;
    my @Fld;
    my $ipart;
    my ($dfx, $dfy, $angast);

    #setting the nx and ny arrays for later use in case
    setup_nx_ny_arrays();

    # and creating the defteb files
    open($fstk,'>', $dftb_stk_file) or 
	die "$die_open_string $dftb_stk_file";
    open($fstk_full,'>', $dftb_stk_full_file) or 
	die "$die_open_string $dftb_stk_full_file";
    open($fstk_full_nptcls,'>', $dftb_stk_full_nptcls_file) or 
	die "$die_open_string $dftb_stk_full_nptcls_file";
    open($fbox,'>', $stat_box_file) or 
	die "$die_open_string $stat_box_file";
    open($fshf_stat,'>', $stat_shf_file) or 
	die "$die_open_string $stat_box_file";
    open($fshf_format_stat,'>', $stat_format_shf_file) or 
	die "$die_open_string $stat_box_file";
    open($fstr,'>', $star_file) or 
	die "$die_open_string $star_file";
    open($fsmv,'>', $simple_movies_file) or 
	die "$die_open_string $simple_movies_file";
    open($fsbx,'>', $simple_Boxfile_file) or 
	die "$die_open_string $simple_Boxfile_file";

    print_header_star_file();

    $cnt_str = 0;
    for ($irows = 0 ; $irows < ($#nparticles) ; $irows++ ) {
	$file_box = "$boxed_file[$irows]";
	#writing to the simple_Boxfile.txt
	print $fsbx qq[$file_box\n];
	#writing to the simple_movies.txt
	$file_movie_txt = $file_box;
	@Fld = split("$ext_als", $file_box);
	@Fld = split("$dir_Boxfile", $Fld[0]);
	$file_movie_txt = $dir_frames_mrc.$Fld[1].$ext_mrc;
	print $fsmv qq[$file_movie_txt\n];

	@file_array_ctf = split("$ext_box", $file_box);
	@file_array_ctf = split("$dir_Boxfile", $file_array_ctf[0]);
	@file_array_ctf = split("$sptr", $file_array_ctf[1]);

	$file_ctf = $dir_CTFfind_txt.$sptr.
	            $file_array_ctf[1].$ext_ctffd.$ext_txt;
	#system("cat $file_ctf");
	parsing_target_CTFfind_pspec_file($file_ctf);

        #Constructing the file name for parsing
	$file_shf = $file_array_ctf[1];
	@file_array_shf = split("$ext_als", $file_shf);

	$file_shf = $dir_Shf.$sptr.$file_array_shf[0].$ext_shf.$ext_txt;
	#system("cat $file_shf");
	parsing_target_Shf_txt_file($file_shf);

	$file_def = "$dir_Deftabs$sptr$file_array_ctf[1]"."$ext_def$ext_txt";
	write_nboxed_2_Deftabs_boxfile($nparticles[$irows],$file_def);

        $dfx = $ctf_values[1] * (0.0001);
        $dfy = $ctf_values[2] * (0.0001);
        $angast=$ctf_values[3];

        for ($ipart = 0 ; $ipart < ($nparticles[$irows]) ; $ipart++ ) {
            print $fstk_full_nptcls qq[dfx=$dfx dfy=$dfy angast=$angast\n];
        }

    }

    close $fbox;
    close $fstk;
    close $fstr;
    close $fstk_full;
    close $fstk_full_nptcls;
    close $fshf_stat;
    close $fshf_format_stat;
    close $fsmv;
    close $fsbx;
}
################################################################################
#subroutine to extract the number of particles in each box files               #
#Frederic D.R. Bonnet, Date: 3rd of Aug. 2015.                                 #
################################################################################
sub write_nboxed_2_Deftabs_boxfile
{
    (my $nptc,$f1)=@_;
    my @Fld;
    my $full_f1;    #file name for fukll output
  
    #constructing the file name for the full output
    @Fld = split("$ext_txt",$f1);
    @Fld = split("$dir_Deftabs",$Fld[0]);
    @Fld = split("$sptr",$Fld[1]);
    $full_f1 = $dir_Deftabs.$sptr.$Fld[1]."_full".$ext_txt;

    #Writting variables files for SIMPLE
    print qq[~~~~~~~~~~~~~~~~~~Writting in the Deftabs dir~~~~~~~~~~~~~~~~~~\n];
    print_writting_message($nptc,$f1);
    print color ('bold blue');
    print qq[The SIMPLE parameters files in repeated $nptc in ];
    print qq[$dir_Deftabs$sptr\n]; print color('reset');
    #write write_2_Deftabs_parsed_files nptc times for full and non-full
    write_2_Deftabs_parsed_files($nptc,$full_f1);
    print qq[~~~~~~~~~~Writting in the $dftb_stk_full_file~~~~~~~~~~~~~~~\n];
    print_writting_message($nptc,$dftb_stk_full_file);
    print color ('bold blue');
    print qq[The SIMPLE parameters files in $dftb_stk_full_file\n];
    print color('reset');
    write_stack_Deftabs_parsed_files();

    print qq[~~~~~~~~~~~~~~~~Writting in the Star file (RELION)~~~~~~~~~~~~~\n];
    print_writting_message($nptc,$star_file);
    print color('bold blue');
    print qq[Writting the RELION parameters for $star_file\n]; print color('reset');
    write_star_files($nptc);

    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to write out what is going on message                              #
#                      Date: 15th of Aug. 2015.                                #
################################################################################
sub print_writting_message {
    my ($nptc, $f1) = @_;
    print qq[Wrtting the parsed line according to the number of particles.  \n];
    print color('bold blue');
    print qq[Number of particles: ]; print color('bold green');
    print qq[$nptc]; print color('reset'); print qq[, ];
    print qq[generates $nptc lines in file: ]; print color('bold green');
    print qq[$f1  \n]; print color('reset');
}
################################################################################
#subroutine to write out the star files for RELION                             #
#                      Date: 14th of Aug. 2015.                                #
################################################################################
sub write_star_files {

    (my $nptc)=@_;

    my @Fld;
    my $index;
    my $amp_con_in = $AMP_CON;
    my $sph_abe_in = $SPH_ABE;
    my $my_file;

    $cnt_str++;
    $index = sprintf("%06d", $cnt_str);
    $my_file = $boxed_file[$cnt_str-1];

    @Fld = split("$dir_Boxfile",$my_file);
    @Fld = split("$sptr",$Fld[1]);
    @Fld = split("$ext_frames",$Fld[1]);

    $my_file = $dir_Data_mrc.$sptr.$Fld[0].$ext_mrc;

    print $fstr qq[$index\@$my_file ];
    print $fstr qq[$ctf_values[1] $ctf_values[2] $ctf_values[3] ];
    print $fstr qq[$AccVlt $amp_con_in $sph_abe_in\n];

    print color ('bold green');
    print qq[counter: $cnt_str\n]; print color('reset');
}
sub print_header_star_file {

    print qq[~~~~~~~~~~~~Writting Header in the Star file (RELION)~~~~~~~~~~\n];
    print qq[Wrtting the Header in the Star file for (RELION): ];
    print color('bold green');
    print qq[$star_file  \n]; print color('reset');

    print $fstr qq[data_images\n];
    print $fstr qq[loop_\n];
    print $fstr qq[_rlnImageName\n];
    print $fstr qq[_rlnDefocusU\n];
    print $fstr qq[_rlnDefocusV\n];
    print $fstr qq[_rlnDefocusAngle\n];
    print $fstr qq[_rlnVoltage\n];
    print $fstr qq[_rlnAmplitudeContrast\n];
    print $fstr qq[_rlnSphericalAberration\n];
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to write out the asc file input for SIMPLE                         #
#                      Date: 14th of Aug. 2015.                                #
################################################################################
sub write_stack_Deftabs_parsed_files{

    my ($dfx, $dfy, $angast, $ctfres, $pixel_in_x, $pixel_in_y);
    my $ival;

    $pixel_in_x = $x_pixSze * (1E+10);                #1.77
    $pixel_in_y = $y_pixSze * (1E+10);                #1.77

    #print qq[$nptc $ctf_values[1] $ctf_values[2] $ctf_values[3]\n];
    #converting from Angstrom to micron meters
    $dfx = $ctf_values[1] * (0.0001);
    $dfy = $ctf_values[2] * (0.0001);
    $angast=$ctf_values[3];
    $ctfres=$ctf_values[6];

    # normal printout 
    print $fstk qq[dfx=$dfx dfy=$dfy angast=$angast ctfres=$ctfres\n];
    # full printout 
    for ($ival = 0 ; $ival < ($#shf_values_x+1) ; $ival++ ) {
	get_xi($ival);
	$shf_values_x_normalised = $shf_values_x[$ival] / $pixel_in_x;
	print $fstk_full qq[$xi$eq$shf_values_x_normalised  ];
    }
    for ($ival = 0 ; $ival < ($#shf_values_y) ; $ival++ ) {
	get_yi($ival);
	$shf_values_y_normalised = $shf_values_y[$ival] / $pixel_in_y;
	print $fstk_full qq[$yi$eq$shf_values_y_normalised  ];
    }
    get_yi(($#shf_values_y));
    $shf_values_y_normalised = $shf_values_y[$#shf_values_y] / $pixel_in_y;
    print $fstk_full qq[$yi$eq$shf_values_y_normalised ];
    print $fstk_full qq[dfx=$dfx dfy=$dfy angast=$angast ctfres=$ctfres\n];

}
################################################################################
#subroutine to write out the the txt files in the Deftabs directory            #
#                      Date: 14th of Aug. 2015.                                #
################################################################################
sub write_2_Deftabs_parsed_files{

    my ($nptc,$full_f1) = @_;
    my($dfx, $dfy, $angast, $ctfres, $pixel_in_x, $pixel_in_y);
    my $iptc;
    my $ival;

    $pixel_in_x = $x_pixSze * (1E+10);                #1.77
    $pixel_in_y = $y_pixSze * (1E+10);                #1.77

    #print qq[$nptc $ctf_values[1] $ctf_values[2] $ctf_values[3]\n];
    #converting from Angstrom to micron meters
    $dfx = $ctf_values[1] * (0.0001);
    $dfy = $ctf_values[2] * (0.0001);
    $angast=$ctf_values[3];
    $ctfres=$ctf_values[6];

    open($fdef,'>', $f1) or die "$die_open_string $f1";
    open($fdef_full,'>', $full_f1) or die "$die_open_string $full_f1";
    for ($iptc = 0 ; $iptc < $nptc ; $iptc++ ) {
	print $fdef qq[dfx=$dfx dfy=$dfy angast=$angast ctfres=$ctfres\n];
	#full write out
	for ($ival = 0 ; $ival < ($#shf_values_x+1) ; $ival++ ) {
	    get_xi($ival);
	    $shf_values_x_normalised = $shf_values_x[$ival] / $pixel_in_x;
	    print $fdef_full qq[$xi$eq$shf_values_x_normalised  ];
	}
	for ($ival = 0 ; $ival < ($#shf_values_y) ; $ival++ ) {
	    get_yi($ival);
	    $shf_values_y_normalised = $shf_values_y[$ival] / $pixel_in_y;
	    print $fdef_full qq[$yi$eq$shf_values_y_normalised  ];
	}
	get_yi(($#shf_values_y));
	$shf_values_y_normalised = $shf_values_y[$#shf_values_y] / $pixel_in_y;
	print $fdef_full qq[$yi$eq$shf_values_y_normalised ];

	print $fdef_full qq[dfx=$dfx dfy=$dfy angast=$angast ctfres=$ctfres\n];
    }
    close $fdef;
    close $fdef_full;

}
################################################################################
#subroutine to print out the xi and yi for a given value                       #
#Frederic D.R. Bonnet, Date: 10th of Aug. 2015.                                #
################################################################################
sub get_xi_yi {
    my ($ival) = @_;
    my $xs = "x";
    my $ys = "y";
    my $count = $ival + 1;
    $xi = $xs.$count;
    $yi = $ys.$count;
    return ($xi,$yi);
}
sub get_xi {
    my ($ival) = @_;
    my $xs = "x";
    my $count = $ival + 1;
    $xi = $xs.$count;
    return ($xi);
}
sub get_yi {
    my ($ival) = @_;
    my $ys = "y";
    my $count = $ival + 1;
    $yi = $ys.$count;
    return ($yi);
}
################################################################################
#subroutine to setup the x,y={1,..,\$#shf_values_(x,y)} arrays according to    #
# the shifted values                                                           #
#Frederic D.R. Bonnet, Date: 10th of Aug. 2015.                                #
################################################################################
sub setup_nx_ny_arrays () {

    my $ival;
    my $xs = "x";
    my $ys = "y";

    my $count = 1;
    for ($ival = 0 ; $ival < ($#shf_values_x+1) ; $ival++ ) {
	$nx_array[$ival] = $xs.$count;
	#print qq[nx_array[$ival]= $nx_array[$ival]\n];
	$count++;
    }
    $count = 1;
    for ($ival = 0 ; $ival < ($#shf_values_y+1) ; $ival++ ) {
	$ny_array[$ival] = $ys.$count;
	#print qq[ny_array[$ival]= $ny_array[$ival]\n];
	$count++;
    }
}
################################################################################
#subroutine to extract the number of particles in each box files               #
#Frederic D.R. Bonnet, Date: 3rd of Aug. 2015.                                 #
################################################################################
sub get_nboxed_2_defTab_boxfile
{
    ($f1)=@_;
    my @Fld;
    print color('bold magenta');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    print qq[Getting the various entries from: ];
    print color('bold green');print qq[$f1\n]; print color('reset');
    open($ffoo,'<', $f1) or die "$die_open_string $f1";
    my $end = 0;
    my $irows = 0;
    while ( $end == 0 )
    {
	$_= <$ffoo>;
	chop $_[0];
	@Fld = split(' ', $_);
	$nparticles[$irows] = $Fld[0];
	$boxed_file[$irows] = $Fld[3];
#	print qq[$irows nptcl: $nparticles[$irows] file: $boxed_file[$irows]];
	if ( $Fld[3] eq 'total' ){$end = 1 ;}
	if ( $end == 0 ) {$irows = $irows + 1;}
    }

    $nrows = $irows;

    close $ffoo;
    print color('bold magenta');       print color('bold green');
    print qq[$f1 ];                    print color('bold magenta');
    print qq[reports nrows: ];         print color('bold green');
    print qq[$nrows ];                 print color('bold magenta');
    print qq[with ];                   print color('bold green');
    print qq[$nrows ];                 print color('bold magenta');
    print qq[files to be analysed.\n];
    print qq[And a ];                  print color('bold green');
    print qq[$boxed_file[$nrows] ];    print color('bold magenta');
    print qq[of ];                     print color('bold green');
    print qq[$nparticles[$nrows] ];    print color('bold magenta');
    print qq[boxed particles\n];
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    print color ('reset');
    return($nrows,@nparticles,@boxed_file);

}
################################################################################
#subroutine to parse the lines from the old pspec files.                       #
#Frederic D.R. Bonnet, Date: 31th of Jul. 2015.                                #
################################################################################
sub parsing_target_Shf_txt_file
{
    ($f1)=@_;
    my $ival;

    if ($VERBOSE == 1) {
	print color('bold yellow');
	print qq[* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n];
	print color('bold cyan');
	print qq[Processing file:\n];
	print qq[$f1\n];
    }

    open($fshf_txt,'<', $f1) or die "$die_open_string $f1";
    my @lines = <$fshf_txt>;
    close $fshf_txt;
    if ($VERBOSE == 1) {
	print color('bold green');
	print qq[$lines[$#lines-1]];
	print qq[$lines[$#lines]];
	print color ('reset');
    }
    @shf_values_x = split(' ',$lines[$#lines-1]);
    @shf_values_y = split(' ',$lines[$#lines]);

    for ($ival = 0 ; $ival < ($#shf_values_x+1) ; $ival++ ) {
	print $fshf_stat qq[$shf_values_x[$ival]  ];
	printf $fshf_format_stat "%12.8f  ",$shf_values_x[$ival];
    }

    for ($ival = 0 ; $ival < ($#shf_values_y) ; $ival++ ) {
	print $fshf_stat qq[$shf_values_y[$ival]  ];
	printf $fshf_format_stat "%12.8f  ",$shf_values_y[$ival];
    }
    print $fshf_stat qq[$shf_values_y[$#shf_values_y]\n];
    printf $fshf_format_stat "%12.8f  \n",$shf_values_y[$#shf_values_y];

    if ($VERBOSE == 1) {
	print color('bold yellow');
	print qq[* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n];
	print color ('reset');
    }
    return (@shf_values_x,@shf_values_y);
}
################################################################################
#subroutine to create the target file for the shifted values                   #
#Frederic D.R. Bonnet, Date: 10th of Aug. 2015.                                #
################################################################################
sub create_shiftval_target
{
    system("ls $dir_Shf$sptr$all$ext_txt > $target_ShiftVal_file$ext_txt");
}
################################################################################
#subroutine to parse the lines from the old pspec files. Includes the recovery #
# in case a file is missing due to unstable computation from unblur and/or     #
#ctffind4                                                                      #
#Frederic D.R. Bonnet, Date: 31th of Jul. 2015.                                #
################################################################################
sub parsing_target_CTFfind_pspec_file 
{
    ($f2)=@_;
    my $ival;
    my $irecover = 0;

    if ($VERBOSE == 1) {
	print color('bold yellow');
	print qq[* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n];
	print color('bold cyan');
	print qq[Processing file:\n];
	print qq[$f2\n];
    }

    $irecover = 0;

RECOVERY:
    #open($fctf,'<', $f2) or die "$die_open_string $f2";
    if (open($fctf,'<', $f2)) {

        #if  ($irecover > 0) {die;}

	my @lines = <$fctf>;
	close $fctf;
	if ($VERBOSE == 1) {
	    print color('bold green');
	    print qq[$lines[$#lines]];
	    print color ('reset');
	}
	@ctf_values = split(' ',$lines[$#lines]);

	for ($ival = 0 ; $ival < ($#ctf_values-1) ; $ival++ ) {
	    print $fbox qq[$ctf_values[$ival]  ];
	}

	print $fbox qq[$ctf_values[$#ctf_values]\n];

    } else {
	print color('bold red');
	print qq[= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n];
	print qq[the file: ];
	print color('bold cyan');print qq[$f2 ];print color('bold red');
	print qq[does not exist!\n];
	print qq[This is because unblur did not generate the file.\n];
	print qq[We will regenerate the file and parse again!!\n];
	print qq[= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n];

	$file_shf_reco = $f2;
	@file_array_shf_recov = split("$ext_als", $file_shf_reco);
	@file_array_shf_recov = split("$dir_CTFfind_txt$sptr",
				      $file_array_shf_recov[0]);

	$file_shf_reco = $dir_frames_mrc.$sptr.
	                 $file_array_shf_recov[1].$ext_mrc;
	system("ls -al $file_shf_reco");
	preprocessing_unblur_one($file_shf_reco);

	$file_ctf_reco = $f2;
	@file_array_ctf_recov = split("$ext_als", $file_ctf_reco);
	@file_array_ctf_recov = split("$dir_CTFfind_txt$sptr",
				      $file_array_ctf_recov[0]);

	$file_ctf_reco = $dir_Als_rot_mrc.$sptr.
	                 $file_array_ctf_recov[1].
	                 $ext_als.$ext_rot.$ext_mrc;
	system("ls -al $file_ctf_reco");
	preprocessing_ctffind_one($file_ctf_reco);

        $irecover = $irecover + 1;

goto RECOVERY;
    }

    if ($VERBOSE == 1) {
	print color('bold yellow');
	print qq[* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n];
	print color ('reset');
    }
    return (@ctf_values);
}
################################################################################
#subroutine to extract the number of particles in each box files               #
#Frederic D.R. Bonnet, Date: 3rd of Aug. 2015.                                 #
################################################################################
sub create_boxfile_target
{
    system("wc $dir_Boxfile$sptr$all$ext_box > $target_Boxfile_file$ext_txt");
}
################################################################################
#subroutine to print message for need to have the Boxfile                      #
#Frederic D.R. Bonnet, Date: 09th of Aug. 2015.                                #
################################################################################
sub print_message_hasbeen_rot_box {
    my @Fld;
    my $my_file;
    @Fld = split("$ext_frames",$targetFrames[0]); 
    $my_file = $Fld[0].$ext_frames.$ext_als.$ext_rot.$ext_box;
    print color('bold red');
    print qq[The Data files in: ];
    print color('bold green');print qq[$dir_Data_mrc \n];
    print color('bold red'); print qq[Have been rotated by an angle: ];
    print color('bold blue');print qq[$DATA_ROT (deg)\n];
    print color('bold red');
    print qq[The files have been moved to the directory: ];
    print color('bold green');print qq[$dir_Data_rot_mrc\n];
    print color('bold red');
    print qq[It is these files that need to be boxed using EMAN2, save all\n];
    print qq[the $all$ext_box files in a directory called Boxfile.\n];
    print qq[\n];
    print qq[Keep the same file structure as the raw files\n];
    print qq[\n];
    print qq[i.e.: ];print color('green');print qq[$my_file\n];
    print color('bold red');
    print qq[\n];
    print qq[You can then rerun this script with the variable:\n];
    print qq[\$FILE_PRE_HANDLER = "extract_data"\n];
    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to print message for need to have the Boxfile                      #
#Frederic D.R. Bonnet, Date: 09th of Aug. 2015.                                #
################################################################################
sub print_message_need_box {
    my @Fld;
    my $my_file;
    @Fld = split("$ext_frames",$targetFrames[0]); 
    $my_file = $Fld[0].$ext_box;
    print color('bold red');
    print qq[The next step is to box particles using EMAN2 and save all of\n];
    print qq[the $all$ext_box files in a directory called Boxfile.\n];
    print qq[\n];
    print qq[Keep the same file structure as the raw files\n];
    print qq[\n];
    print qq[i.e.: $my_file\n];
    print qq[\n];
    print qq[You can then rerun this script with the variable:\n];
    print qq[\$FILE_PRE_HANDLER = "extract_data"\n];
    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to print message for need to have the Boxfile                      #
#Frederic D.R. Bonnet, Date: 09th of Aug. 2015.                                #
################################################################################
sub print_message_have_box {
    print color('bold red');
    print qq[You have selected you have the box files ready to be parsed.\n];
    print qq[The next step is to parse these to extract variables\n];
    print qq[the $all$ext_box files in a directory called Boxfile.\n];
    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to preprocess with unblur                                          #
#Frederic D.R. Bonnet, Date: 9th of Aug. 2015.                                 #
################################################################################
sub preprocessing_unblur {

    $rundir_in = $UNBLUR_DIR;           #/opt/Devel_tools/unblur_1.0/bin/
    $file_framename_in = $target_frames_file.$ext_txt; #target_frames.txt
    $nframes_in = $NFRAMES;             #7
    $pixel_in = $x_pixSze * (1E+10);    #1.77

    print_input_var_unblur();
    system("csh unblur_alignSum_all.csh $rundir_in $file_framename_in $nframes_in $pixel_in");
    print color('bold red');
    print qq[unblur_alignSum_all.csh is Done...\n];
    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];

}
################################################################################
#subroutine to preprocess with unblur                                          #
#Frederic D.R. Bonnet, Date: 9th of Aug. 2015.                                 #
################################################################################
sub preprocessing_ctffind {

    $rundir_in = $CTFFIND_DIR;          #/opt/Devel_tools/CTFfind/
    $file_framename_in = $target_als_rot_file.$ext_txt;
    $pixel_in = $x_pixSze * (1E+10);    #1.77
    $AccVlt_in = $AccVlt;               #300 (KeV) acceleration volatage
    $sph_abe_in = $SPH_ABE;             #sphererical aberration 2.7 (default)
    $amp_con_in = $AMP_CON;             #Amplitude constrast default 0.07
    $sze_pwr_spc_in = $SZE_PWR_SPC;     #Size of prower spectrum default 512

    print_input_var_ctffind();
    system("csh ctffind_alignSum_all.csh $rundir_in $file_framename_in $pixel_in $AccVlt_in $sph_abe_in $amp_con_in $sze_pwr_spc_in");
    print color('bold red');
    print qq[ctffind_alignSum_all.csh Done...\n];
    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];

}
################################################################################
#subroutine to preprocess with unblur                                          #
#Frederic D.R. Bonnet, Date: 9th of Aug. 2015.                                 #
################################################################################
sub preprocessing_unblur_one {
    ($f1) = @_;
    $rundir_in = $UNBLUR_DIR;           #/opt/Devel_tools/unblur_1.0/bin/
    $file_framename_in = $f1;           #input file
    $nframes_in = $NFRAMES;             #7
    $pixel_in = $x_pixSze * (1E+10);    #1.77

    print_input_var_unblur();
    system("csh unblur_alignSum_one.csh $rundir_in $file_framename_in $nframes_in $pixel_in");
    print color('bold red');
    print qq[unblur_alignSum_all.csh is Done...\n];
    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];

}
################################################################################
#subroutine to preprocess with unblur                                          #
#Frederic D.R. Bonnet, Date: 9th of Aug. 2015.                                 #
################################################################################
sub preprocessing_ctffind_one {
    ($f1) = @_;
    $rundir_in = $CTFFIND_DIR;          #/opt/Devel_tools/CTFfind/
    $file_framename_in = $f1;           # input file
    $pixel_in = $x_pixSze * (1E+10);    #1.77
    $AccVlt_in = $AccVlt;               #300 (KeV) acceleration volatage
    $sph_abe_in = $SPH_ABE;             #sphererical aberration 2.7 (default)
    $amp_con_in = $AMP_CON;             #Amplitude constrast default 0.07
    $sze_pwr_spc_in = $SZE_PWR_SPC;     #Size of prower spectrum default 512

    print_input_var_ctffind();
    system("csh ctffind_alignSum_one.csh $rundir_in $file_framename_in $pixel_in $AccVlt_in $sph_abe_in $amp_con_in $sze_pwr_spc_in");
    print color('bold red');
    print qq[ctffind_alignSum_all.csh Done...\n];
    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];

}
################################################################################
#subroutine to set the values for the stanp variables organiser values         #
#Frederic D.R. Bonnet, Date: 9th of Aug. 2015.                                 #
################################################################################
sub set_organiser_file {
    my ($organiser) = @_;

    $f1 = $target_organiser_file.$ext_txt;

    open($forg,"> $f1") or die ("cannot open $f1");

    print qq[Generating the target file: ];
    print color('bold green');print qq[$f1\n];print color('reset');
    print qq[Setting the organiser variable 0:not orgnaised 1:organised\n];
    print qq[organiser: $organiser\n];

    $orga  = $organiser;
    print $forg qq[$orga $all\n];
    close $forg;

    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    return ($orga);
}

################################################################################
#subroutine to set the input values for the init_rot_Data.csh                  #
#Frederic D.R. Bonnet, Date: 18th of Aug. 2015.                                #
################################################################################
sub print_input_var_initial_rotation {

    print color('bold red');
    print qq[Input variables from StripData.pl for init_rot_Data.csh :\n];

    print color('reset'); print qq[Rotated Data dir: ];
    print color('bold green'); print qq[$rundir_in\n];

    print color('reset');print qq[Data filename: ];print color('bold yellow');
    print qq[ $file_framename_in \n];

    print color('reset'); print qq[Extension name: ];
    print color('bold yellow'); print qq[$ext_in\n];

    print color('reset'); print qq[Amount of rotation (anti-clockwise): ];
    print color('bold blue'); print qq[$deg_in (deg)\n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to set the input values for the unblur program                     #
#Frederic D.R. Bonnet, Date: 8th of Aug. 2015.                                 #
################################################################################
sub print_input_var_ctffind {

    print color('bold red');
    print qq[Input variables from StripData.pl for ];
    print qq[ctffind_alignSum_all.csh :\n]; print color('reset');

    print qq[Unblur exec dir: ]; print color('bold green');
    print qq[$rundir_in\n]; print color('reset');

    print qq[Frame filename: ];print color('bold yellow');
    print qq[ $file_framename_in \n]; print color('reset'); 

    print qq[Pixel x: ]; print color('bold cyan'); print qq[$pixel_in\n];
    print color('reset'); 

    print qq[Acceleration Voltage (KeV): ];
    print color('bold cyan'); print qq[$AccVlt_in \n];print color('reset');

    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to set the input values for the unblur program                     #
#Frederic D.R. Bonnet, Date: 8th of Aug. 2015.                                 #
################################################################################
sub print_input_var_unblur {

    print color('bold red');
    print qq[Input variables from StripData.pl for ];
    print qq[unblur_alignSum_all.csh :\n]; print color('reset');

    print qq[Unblur exec dir: ]; print color('bold green');
    print qq[$rundir_in\n]; print color('reset');

    print qq[Frame filename: ];print color('bold yellow');
    print qq[ $file_framename_in \n]; print color('reset');

    print qq[Number of frames: ]; print color('bold blue');
    print qq[$nframes_in \n]; print color('reset');

    print qq[Pixel x: ]; print color('bold cyan'); print qq[$pixel_in\n];
    print color('reset');

    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to parse the lines from the old pspec files.                       #
#Frederic D.R. Bonnet, Date: 8th of Aug. 2015.                                 #
################################################################################
sub parse_XML_file {
    # only need the kfirst one or two files in the list as the parsed param are
    # are fixed for the entire measurment

    for ($irows = 0; $irows < ($#targetFrames-($#targetFrames-2)); $irows++ ) {
	#constructing files with the directory structure
	@xml_file_array = split("$ext_mrc",$targetFrames[$irows]);
	@xml_file_array = split("$ext_frames",$xml_file_array[0]);

	#parsing the xml files.
	if ($has_been_organised == 0) {
	    $xml_file = "$xml_file_array[0]$ext_xml";	    
	} elsif ($has_been_organised == 1) {
	    $xml_file = "$dir_xml$sptr$xml_file_array[0]$ext_xml";
	}

	#print qq[$xml_file\n];
	parsing_target_xml_file($xml_file);

    }
    print color('bold cyan');
    print qq[The parsed PixelSize x: $x_pixSze, y: $y_pixSze, ];
    print qq[Nominal Magnification: $TemMag and ];
    print qq[Acceleration Voltage (KeV): $AccVlt\n];
    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];

}
#parsing method
sub parsing_target_xml_file
{
    ($f2)=@_;

    print color('bold yellow');
    print qq[* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n];
    print color('bold cyan');
    print qq[Processing file:\n];
    print qq[$f2\n];

    open($fxml,'<', $f2) or die "$die_open_string $f2";
    @lines = <$fxml>;
    close $fxml;
    print color('magenta');
    #print qq[$lines[$#lines]\n];
    print color ('reset');

    print qq[~~~~~~~~~~~~~~~~~~~~~~Parsing XML file~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    parse_PixelSize(@lines);
    parse_TemMagnification(@lines);
    parse_AccelerationVoltage(@lines);
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];

    print color('bold yellow');
    print qq[* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n];
    print color ('reset');
}
#The parsers
sub parse_AccelerationVoltage {
    (@lines)=@_;
    #doing the parsing, creating our parser object
    my $psr = new XML::Parser(
	Handlers => {
	    Start => \&hdl_start_AccVlt, End => \&hdl_end_AccVlt,
	    Char  => \&hdl_char,     Default => \&hdl_def, });
    $psr->parse($lines[$#lines]);
}
sub parse_TemMagnification {
    (@lines)=@_;
    #doing the parsing, creating our parser object
    my $psr = new XML::Parser(
	Handlers => {
	    Start => \&hdl_start_TemMag, End => \&hdl_end_TemMag,
	    Char  => \&hdl_char,     Default => \&hdl_def, });
    $psr->parse($lines[$#lines]);
}
sub parse_PixelSize {
    (@lines)=@_;
    #doing the parsing, creating our parser object
    my $psr = new XML::Parser(
	Handlers => {
	    Start => \&hdl_start_pixSze, End => \&hdl_end_pixSze,
	    Char  => \&hdl_char,     Default => \&hdl_def, });
    $psr->parse($lines[$#lines]);
}
#The Handlers
sub hdl_start_AccVlt {
    my ($p, $elt, %atts) = @_;
    return unless $elt eq 'AccelerationVoltage';
    $atts{'_str'} = '';
    $message = \%atts;
}
sub hdl_end_AccVlt {
    my ($p, $elt) = @_;
    format_message_AccVlt($message) if $elt eq 'AccelerationVoltage' 
	&& $message && $message->{'_str'} =~ /\S/;
}
sub format_message_AccVlt {
    my $atts = shift;
    $atts->{'_str'} =~ s/\n//g;
    $AccVlt = $atts->{'_str'} / 1000;
    print "Acceleration Voltage (KeV): ";
    print color('bold cyan'); print "$AccVlt\n";
    print color ('reset');
}
sub hdl_start_TemMag {
    my ($p, $elt, %atts) = @_;
    return unless $elt eq 'TemMagnification';
    $atts{'_str'} = '';
    $message = \%atts;
}
sub hdl_end_TemMag {
    my ($p, $elt) = @_;
    format_message_TemMag($message) if $elt eq 'TemMagnification'
	&& $message && $message->{'_str'} =~ /\S/;
}
sub format_message_TemMag {
    my $atts = shift;
    $atts->{'_str'} =~ s/\n//g;
    $TemMag = $atts->{'_str'};
    print "Nominal Magnification: ";
    print color('bold cyan'); print "$TemMag\n";
    print color('reset');
}
sub hdl_start_pixSze {
    my ($p, $elt, %atts) = @_;
    return unless $elt eq 'pixelSize';
    $atts{'_str'} = '';
    $message = \%atts;
}
sub hdl_end_pixSze {
    my ($p, $elt) = @_;
    format_message_pixSze($message) if $elt eq 'pixelSize' && $message
	&& $message->{'_str'} =~ /\S/;
}
sub format_message_pixSze {
    my $atts = shift;
    $atts->{'_str'} =~ s/\n//g;

    my @Fld = split('1m',$atts->{'_str'});
    $x_pixSze = $Fld[0];
    $y_pixSze = $Fld[1];
    print "PixelSize(Angstroms): x: ";
    print color('bold cyan'); print "$x_pixSze ";
    print color('reset'); print "y: ";
    print color('bold cyan'); print "$y_pixSze\n";
    print color('reset');
}
#common handlers
sub hdl_char {
    my ($p, $str) = @_;
    $message->{'_str'} .= $str;
}
sub hdl_def{}
################################################################################
#subroutine to generate the target file for the Frames folder                  #
#Frederic D.R. Bonnet, Date: 08th of Aug. 2015.                                #
################################################################################
sub generate_Target_data_list
{
    ($f1,$dir,$ext)=@_;
    print qq[Generating the target Data_MRC file: ];
    print color('bold green');print qq[$f1\n];print color('reset');

    system("ls $dir$sptr$all$ext > $f1");
}
################################################################################
#subroutine to generate the target file for the Frames folder                  #
#Frederic D.R. Bonnet, Date: 08th of Aug. 2015.                                #
################################################################################
sub generate_Target_frames_list
{
    ($f1,$dir,$ext)=@_;
    print qq[Generating the target frames file: ];
    print color('bold green');print qq[$f1\n];print color('reset');

    system("ls $dir$sptr$all$ext > $f1");
}
################################################################################
#subroutine to generate the target file for the Als_rot folder                 #
#Frederic D.R. Bonnet, Date: 11th of Aug. 2015.                                #
################################################################################
sub generate_Target_als_rot_list
{
    ($f1,$dir,$ext)=@_;
    print qq[Generating the target als rot file: ];
    print color('bold green');print qq[$f1\n];print color('reset');
    system("ls $dir$sptr$all$ext > $f1");
}
################################################################################
#subroutine to organise the data in a manageable directory structure           #
#Frederic D.R. Bonnet, Date: 9th of Aug. 2015.                                 #
################################################################################
sub unorganise {
    $has_been_organised = 0;
    set_organiser_file($has_been_organised);
}
################################################################################
#subroutine to organise the data in a manageable directory structure           #
#Frederic D.R. Bonnet, Date: 9th of Aug. 2015.                                 #
################################################################################
sub organise {
    organise_xml(); organise_jpg();
    #Here the order matters
    organise_frames_mrc(); # 1st the *frames.mrc files
    organise_Data_mrc();   # 2nd the Data*.mrc files
    organise_Othr_mrc();   # 3rd the Other *.mrc files
    #Finaly the TargetLocations and GridSquare
    organise_TargetLocation_dm(); organise_GridSquare_dm();
    $has_been_organised = 1;
    set_organiser_file($has_been_organised);
}
################################################################################
#subroutine to reorganise the Data into its original form                      #
#Frederic D.R. Bonnet, Date: 9th of Aug. 2015.                                 #
################################################################################
sub reorganise {
    reorganise_xml(); reorganise_jpg();
    #Here the order matters
    reorganise_frames_mrc(); # 1st the *frames.mrc files
    reorganise_Data_mrc();   # 2nd the Data*.mrc files
    reorganise_Othr_mrc();   # 3rd the Other *.mrc files
    #Finaly the TargetLocations and GridSquare
    reorganise_TargetLocation_dm(); reorganise_GridSquare_dm();
    if ( $DELETE_DIR == 1 ) {
	delete_PreprocessDirectories();
    }
    $has_been_organised = 0; #stamp variables to keep track organisation
    set_organiser_file($has_been_organised);
}
################################################################################
#subroutine to delete the preprocess directories structure for the file system #
#Frederic D.R. Bonnet, Date: 08th of Aug. 2015.                                #
################################################################################
sub delete_PreprocessDirectories {
    my $ans;
    my $final_ans;

    print color('bold red');
    print qq[!!!!!!!!!!! Warning you have selected to delete !!!!!!!!!!!!!!!\n];
    print qq[!!!!!!!!!!!     the prepocess directories!?     !!!!!!!!!!!!!!!\n];
    print qq[\n];
    print qq[Have you made sure to have emptied the following directories?  \n];
    print qq[                      "];print color('reset');
    print qq[yes];print color ('bold red');
    print qq[" or "]; print color('reset');print qq[no];print color('bold red');
    print qq["                            \n];
    print color('bold cyan');  print qq[$dir_xml\n]; 
    print color('bold blue');  print qq[$dir_jpg\n]; 
    print color('bold yellow');print qq[$dir_frames_mrc\n$dir_Data_mrc\n]; 
    print color('bold yellow');print qq[$dir_Othr_mrc\n];
    print color('bold magenta');print qq[$dir_TargetLocation\n$dir_GridSquare\n];

    print color('reset');
    $ans= <>;
    print color('bold red');
    if ($ans =~ /es/ ) {
	print qq[Directories and all of its content will be deleted!!!\n];
	print qq[Proceed: "yes" or no"?\n];
	print color('reset');
	$final_ans = <>;
	print color('bold red');
	if ($final_ans =~ /es/ ) {
	    print qq[Deleting the directories...\n];
	    print color('bold cyan');  print qq[$dir_xml ]; 
	    print color('bold blue');  print qq[$dir_jpg ]; 
	    print color('bold yellow');print qq[$dir_frames_mrc $dir_Data_mrc]; 
	    print color('bold yellow');print qq[$dir_Othr_mrc ];
	    print color('bold magenta');print qq[$dir_TargetLocation $dir_GridSquare\n];

	    system("rm -rf $dir_xml $dir_jpg $dir_frames_mrc $dir_Data_mrc $dir_Othr_mrc $dir_TargetLocation $dir_GridSquare");
	    print color('bold red');
	    print qq[Done...\n];
	} elsif ($final_ans =~ /no/ ) {
	    print color('reset');
	    print qq[Make sure that you have run the script with the option\n];
	    print qq[\$FILE_ORGANISE_DIR = "re" in simple_StripData_input.pm\n];	    
	}
	print qq[Do you want to delete the folowing directories as well?  \n];
	print qq[                      "];print color('reset');
	print qq[yes];print color ('bold red');
	print qq[" or "]; print color('reset');print qq[no];print color('bold red');
	print qq["                            \n];
	print color('bold yellow'); print qq[$dir_Als_rot_mrc\n]; 
	print color('bold yellow'); print qq[$dir_CTFfind\n]; 
	print color('bold white');  print qq[$dir_CTFfind_txt\n]; 
	print color('bold white');  print qq[$dir_Shf\n]; 
	print color('bold yellow'); print qq[$dir_Deftabs\n];
	print color('bold yellow');  print qq[$dir_Rotd_Frames_Als\n];
	print color('bold yellow');  print qq[$dir_Data_rot_mrc\n];

	print color('reset');
	$ans= <>;
	print color('bold red');
	if ($final_ans =~ /es/ ) {
	    print qq[Deleting the directories...\n];
	    print color('bold yellow'); print qq[$dir_Als_rot_mrc ];
	    print color('bold yellow'); print qq[$dir_CTFfind ];
	    print color('bold white');  print qq[$dir_CTFfind_txt ];
	    print color('bold white');  print qq[$dir_Shf ];
	    print color('bold yellow'); print qq[$dir_Deftabs ];
	    print color('bold yellow');  print qq[$dir_Rotd_Frames_Als\n];
	    system("rm -rf $dir_Als_rot_mrc $dir_CTFfind $dir_CTFfind_txt $dir_Shf $dir_Deftabs $dir_Rotd_Frames_Als $dir_Data_rot_mrc");
	    print color('bold red');
	    print qq[Done...\n];
	} elsif ($final_ans =~ /no/ ) {
	    print color('reset');
	    print qq[Make sure that you have run the script with the option\n];
	    print qq[\$FILE_ORGANISE_DIR = "re" in simple_StripData_input.pm\n];	    
	}

    } elsif ($ans =~ /no/ ) {
	print color('reset');
	print qq[Make sure that you have run the script with the option\n];
	print qq[\$FILE_ORGANISE_DIR = "re" in simple_StripData_input.pm\n];
    }

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to delete the preprocess directories structure for the file system #
#Frederic D.R. Bonnet, Date: 15th of Aug. 2015.                                #
################################################################################
sub delete_PreprocessFiles {
    my $ans;
    my $final_ans;

    print color('bold red');
    print qq[!!!!!!!!!!! Warning you have selected to delete !!!!!!!!!!!!!!!\n];
    print qq[!!!!!!!!!!!         the prepocess files!?       !!!!!!!!!!!!!!!\n];
    print qq[\n];
    print qq[Have you made sure that you do not need these files?           \n];
    print qq[                      "];print color('reset');
    print qq[yes];print color ('bold red');
    print qq[" or "]; print color('reset');print qq[no];print color('bold red');
    print qq["                            \n];
    print color('bold white'); print qq[all the $all$ext_txt\n];
    print color('bold red');   print qq[all the $all$ext_asc\n];
    print color('bold red');   print qq[all the $all$ext_log\n];
    print color('bold red');   print qq[all the $all$ext_str\n];

    print color('reset');
    $ans= <>;
    print color('bold red');
    if ($ans =~ /es/ ) {
	print qq[Preporcessing generated files will be deleted!!!\n];
	print qq[Proceed: "yes" or no"?\n];
	print color('reset');
	$final_ans = <>;
	print color('bold red');
	if ($final_ans =~ /es/ ) {
	    print qq[Deleting the directories...\n];
	    print color('bold white'); print qq[all the $all$ext_txt\n];
	    print color('bold red');   print qq[all the $all$ext_asc\n];
	    print color('bold red');   print qq[all the $all$ext_log\n];
	    print color('bold red');   print qq[all the $all$ext_str\n];
	    system("rm $all$ext_txt $all$ext_asc $all$ext_log $all$ext_str");
	    print color('bold red');
	    print qq[Done...\n];
	} elsif ($final_ans =~ /no/ ) {
	    print color('reset');
	    print qq[Make sure that you have run the script with the option\n];
	    print qq[\$FILE_ORGANISE_DIR = "re" in simple_StripData_input.pm\n];	    
	}

    } elsif ($ans =~ /no/ ) {
	print color('reset');
	print qq[Make sure that you have run the script with the option\n];
	print qq[\$FILE_ORGANISE_DIR = "re" in simple_StripData_input.pm\n];
    }

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the GridSquare*.dm back to execdir                         #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub reorganise_GridSquare_dm {

    print color('bold magenta');
    print qq[Moving all of the GridSquare*.dm files from $dir_GridSquare back: $execdir ... \n];
    system("mv $execdir$sptr$dir_GridSquare$sptr$all$ext_dmf $execdir");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the GridSquare*.dm to TargetLocation folder                #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub organise_GridSquare_dm {

    print color('bold magenta');
    print qq[Creating $dir_GridSquare directory for $ext_Grd$all$ext_dmf files\n];
    system("mkdir $execdir$sptr$dir_GridSquare");

    print qq[Moving all $ext_Grd$all$ext_dmf files to $dir_GridSquare ...\n];
    system("mv $execdir$sptr$all$ext_Grd$all$ext_dmf $execdir$sptr$dir_GridSquare");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the other files*.dm back to execdir                        #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub reorganise_TargetLocation_dm {

    print color('bold magenta');
    print qq[Moving all of the Other *.dm files from $dir_TargetLocation back: $execdir ... \n];
    system("mv $execdir$sptr$dir_TargetLocation$sptr$all$ext_dmf $execdir");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the Other files*.dm to TargetLocation folder               #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub organise_TargetLocation_dm {

    print color('bold magenta');
    print qq[Creating TargetLocation directory for TargetLocation*.dm files\n];
    system("mkdir $execdir$sptr$dir_TargetLocation");

    print qq[Moving all TargetLocation*.dm files to $dir_TargetLocation ...\n];
    system("mv $execdir$sptr$all$ext_Trg$all$ext_dmf $execdir$sptr$dir_TargetLocation");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the other files*.dm back to execdir                        #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub reorganise_Othr_mrc {

    print color('bold yellow');
    print qq[Moving all of the Other *.mrc files from $dir_Othr_mrc back: $execdir ... \n];
    system("mv $execdir$sptr$dir_Othr_mrc$sptr$all$ext_mrc $execdir");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the Other files*.dm to TargetLocation folder               #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub organise_Othr_mrc {

    print color('bold yellow');
    print qq[Creating $dir_Othr_mrc directory for other $all.$ext_mrc files\n];
    system("mkdir $execdir$sptr$dir_Othr_mrc");

    print qq[Moving all $all.$ext_mrc files to $dir_Othr_mrc ...\n];
    system("mv $execdir$sptr$all$ext_mrc $execdir$sptr$dir_Othr_mrc");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the files*_Data.mrc back to execdir                        #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub reorganise_Data_mrc {

    print color('bold yellow');
    print qq[Moving all of the *_Data.mrc files from $dir_Data_mrc back: $execdir ... \n];
    system("mv $execdir$sptr$dir_Data_mrc$sptr$all$ext_mrc $execdir");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the files*_Data.mrc to Frames_MRC folder                   #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub organise_Data_mrc {

    print color('bold yellow');
    print qq[Creating the MRC directory for all of the *_Data.mrc files\n];
    system("mkdir $execdir$sptr$dir_Data_mrc");

    print qq[Moving all of the *_Data*.mrc files to $dir_Data_mrc ...\n];
    system("mv $execdir$sptr$all$ext_Data_mrc$all$ext_mrc $execdir$sptr$dir_Data_mrc");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the files*_frames.mrc back to execdir                      #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub reorganise_frames_mrc {

    print color('bold yellow');
    print qq[Moving all of the *_frames.mrc files from $dir_frames_mrc back: $execdir ... \n];
    system("mv $execdir$sptr$dir_frames_mrc$sptr$all$ext_mrc $execdir");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the files*_frames.mrc to Frames_MRC folder                 #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub organise_frames_mrc {

    print color('bold yellow');
    print qq[Creating the MRC directory for all of the *_frames.mrc files\n];
    system("mkdir $execdir$sptr$dir_frames_mrc");

    print qq[Moving all of the *.mrc files to $dir_frames_mrc ...\n];
    system("mv $execdir$sptr$all$ext_frames_mrc $execdir$sptr$dir_frames_mrc");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the files.jpg to JPG folder                                #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub reorganise_jpg {

    print color('bold blue');
    print qq[Moving all of the *.jpg files from $dir_jpg back: $execdir ... \n];
    system("mv $execdir$sptr$dir_jpg$sptr$all$ext_jpg $execdir");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the files.jpg to JPG folder                                #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub organise_jpg {

    print color('bold blue');
    print qq[Creating the JPG directory for all of the *.jpg files\n];
    system("mkdir $execdir$sptr$dir_jpg");

    print qq[Moving all of the *.jpg files to $dir_jpg ...\n];
    system("mv $execdir$sptr$all$ext_jpg $execdir$sptr$dir_jpg");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the files.xml to XML folder                                #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub reorganise_xml {

    print color('bold cyan');
    print qq[Moving all of the *.xml files from $dir_xml back: $execdir ... \n];
    system("mv $execdir$sptr$dir_xml$sptr$all$ext_xml $execdir");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to move the files.xml to XML folder                                #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub organise_xml {

    print color('bold cyan');
    print qq[Creating the XML directory for all of the *.xml files\n];
    system("mkdir $execdir$sptr$dir_xml");

    print qq[Moving all of the *.xml files to $dir_xml ...\n];
    system("mv $execdir$sptr$all$ext_xml $execdir$sptr$dir_xml");
    print qq[Done... \n];

    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to print organised message                                         #
#Frederic D.R. Bonnet, Date: 08th of Aug. 2015.                                #
################################################################################
sub print_preprocessing_message {
    print color('bold red');
    print qq[Now processing the Data using unblur and CTFfinder4\n];
    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to print organised message                                         #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub print_reorganised_message {
    print color('bold red');
    print qq[Reorganising the files back to their original dir structure\n];
    print qq[in a suitable way\n];
    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to print organised message                                         #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub print_organised_message {
    print color('bold red');
    print qq[Organising the directory structure in a suitable way\n];
    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to print unorganised message                                       #
#Frederic D.R. Bonnet, Date: 07th of Aug. 2015.                                #
################################################################################
sub print_unorganised_message {
    print color('bold red');
    print qq[Leaving the directory structure untouched\n];
    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to check the path if the path has spaces in it add \ in the spaces #
#so that it can be processed.                                                  #
#Frederic D.R. Bonnet, Date: 06th of Aug. 2015.                                #
################################################################################
sub get_check_fix_path {
    (my $path) = @_;
    my $t_p;
    my $ifld;

    my @Fld = split(' ', $path);

    $true_path = "";
    for ($ifld = 0 ; $ifld < ($#Fld+1) ; $ifld++ ) {
	$true_path .= $Fld[$ifld].$b_sptr." ";
    }

    my @ta = split (//,$true_path);
    for ($ifld = 0 ; $ifld < ($#ta-1) ; $ifld++ ) {
	$t_p .= $ta[$ifld];
    }

    $true_path = $t_p;

    return($true_path);
}
################################################################################
#subroutine to distribute the file into their current locations after the      #
#directories have been created                                                 #
#Frederic D.R. Bonnet, Date: 11th of May. 2015.                                #
################################################################################
sub distribute_files { 
    my $iemp;

    #first setting the organiser to default value 0
    $has_been_organised = 0;  #setting the organiser file 0:default
    set_organiser_file($has_been_organised);

    if ( ($#nonemptyString+1) != ($#nonemptyTarget+1) ) {die};

    for ($iemp = 0 ; $iemp < ($#nonemptyTarget+1) ; $iemp++ ) {
	system("mv $nonemptyString[$iemp][0] $nonemptyTarget[$iemp]");
	system("mv $nonemptyString[$iemp][1] $nonemptyTarget[$iemp]");
	system("mv $nonemptyString[$iemp][2] $nonemptyTarget[$iemp]");
#	system("ls -R $nonemptyString[$iemp][0] $nonemptyTarget[$iemp]");
    }
    system("mkdir TargetLocation GridSquare");
    system("mv TargetLocation_* TargetLocation");
    system("mv GridSquare_* GridSquare");
}
################################################################################
#subroutine to regroup the file into their current locations after the         #
#directories have been created                                                 #
#Frederic D.R. Bonnet, Date: 11th of May. 2015.                                #
################################################################################
sub regroup_files { 
    my $irows;
    
    #first setting the organiser to default value 0
    $has_been_organised = 0;  #setting the organiser file 0:default
    set_organiser_file($has_been_organised);

    for ($irows = 0 ; $irows < $nrows ; $irows++ ) {
	system("mv $targetlocation[$irows]$sptr$all $execdir");
#	system("ls -R *");
    }
    system("mv TargetLocation/TargetLocation_* $execdir");
    system("mv GridSquare/GridSquare_* $execdir");
}
################################################################################
#subroutine to create the directory structure for the file system              #
#Frederic D.R. Bonnet, Date: 11th of May. 2015.                                #
################################################################################
sub create_TargetDirectories {
    my $iemp;

    if ( ($#nonemptyString+1) != ($#nonemptyTarget+1) ) {die};

    for ($iemp = 0 ; $iemp < ($#nonemptyTarget+1) ; $iemp++ ) {
	#print qq[$iemp $nonemptyTarget[$iemp]\n];
	system("mkdir $nonemptyTarget[$iemp]");
    }    
}
################################################################################
#subroutine to delete the directory structure for the file system              #
#Frederic D.R. Bonnet, Date: 11th of May. 2015.                                #
################################################################################
sub delete_TargetDirectories {
    my $iemp;

    if ( ($#nonemptyString+1) != ($#nonemptyTarget+1) ) {die};

    for ($iemp = 0 ; $iemp < ($#nonemptyTarget+1) ; $iemp++ ) {
	#print qq[$iemp $nonemptyTarget[$iemp]\n];
	system("rm -rf $nonemptyTarget[$iemp]");
    }
    #deleting the unwanted directories TargetLocation and  
    system("rm -rf TargetLocation GridSquare");
}
################################################################################
#subroutine to scan existing files in the current dierctory                    #
#Frederic D.R. Bonnet, Date: 11th of May. 2015.                                #
################################################################################
sub scan_directory { 
    my $irows;
    my $string_mrc;
    my $string_xml;
    my $string_jpg;
    my $rc;            #return code for system call

    my $iemp;
    
    my $dataQualityfile = "quality.asc";
    open($fscn,"> $dataQualityfile") or die ("cannot open $dataQualityfile ");

    $iemp = 0;
    for ($irows = 0 ; $irows < $nrows ; $irows++ ) {
	$string_mrc = "$execdir$sptr$all$targetlocation[$irows]$all".".mrc";
	$string_xml = "$execdir$sptr$all$targetlocation[$irows]$all".".xml";
	$string_jpg = "$execdir$sptr$all$targetlocation[$irows]$all".".jpg";

	$datastatTarget[$irows] = 0;

	$rc = system("ls -sF $string_mrc");
	if ($rc == 0 ) {
	    $nonemptyString[$iemp][0] = $string_mrc;
	    $nonemptyString[$iemp][1] = $string_xml;
	    $nonemptyString[$iemp][2] = $string_jpg;
	    $nonemptyTarget[$iemp] = $targetlocation[$irows];
	    $datastatTarget[$irows] = 
		$targetlocation[$irows]/$targetlocation[0];
	    $iemp ++;
	}

	printf $fscn "%6i  %15i\n",$irows,$datastatTarget[$irows];
    }
    #Before making a mess check for coherence
    if ( ($#nonemptyString+1) != ($#nonemptyTarget+1) ) {die};
    #saving a log to disk for the operations
    for ($iemp = 0 ; $iemp < ($#nonemptyString+1) ; $iemp++) {
	system("ls -sF $nonemptyString[$iemp][0] >> nonemptyString.log");
	system("echo TargetLocation_$nonemptyTarget[$iemp].dm >> nonemptyTarget.log ");
    }

    close $fscn;
    my $stat = int( ($iemp / $nrows) *100.0); 
    my $title="Data Quality control: $stat\% usable ($iemp".'/'."$nrows files)"; 
    my $xaxis = 'number files';
    my $yaxis = 'Sample tag files';
    my $yaxis_lowHigh = '--y-limits 0 2 1';
    my $axes_lab = "--x-label \"$xaxis\" --y-label \"$yaxis\"";
    my $lines = "--line-mode 0 --symbol 4";
    system("$graph -L \"$title\" $axes_lab $yaxis_lowHigh $lines < $dataQualityfile&");

    return (@nonemptyTarget, @nonemptyString);
}
################################################################################
#subroutine to generate the target file                                        #
#Frederic D.R. Bonnet, Date: 11th of May. 2015.                                #
################################################################################
sub generate_TargetLocation_list
{
    ($f1)=@_;
    print qq[Generating the target file: ];
    print color('bold green');print qq[$f1\n];print color('reset');
    system("ls TargetLocation_* > $f1");
}
################################################################################
#subroutine to extract the prices from the old TickData files.                 #
#Frederic D.R. Bonnet, Date: 31st of Jul. 2015.                                #
################################################################################
sub get_target_Data_MRC_file
{
    ($f1)=@_;
    my @Fld;
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    print qq[Getting the various entries from the file name: ];
    print color('bold green');print qq[$f1\n];print color('reset');
    print qq[We are now storing the details from file: $f1                  \n];
    open($fdta_mrc,'<', $f1) or die "$die_open_string $f1";
    my $end = 0;
    my $irows = 0;
    while ( $end == 0 )
    {
	$_= <$fdta_mrc>;
	chop $_[0];
	@Fld = split("$sptr", $_);
	$targetDataMRC[$irows] = $Fld[1];
	#print qq[$irows $targetDataMRC[$irows]];
	if ( $Fld[0] eq '' ){$end = 1 ;}
	if ( $end == 0 ) {$irows = $irows + 1;}
    }

    $nrows = $irows;

    close $fdta_mrc;
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    return($nrows,@targetDataMRC);
}
################################################################################
#subroutine to extract the prices from the old TickData files.                 #
#Frederic D.R. Bonnet, Date: 31st of Jul. 2015.                                #
################################################################################
sub get_target_oragniser_file
{
    ($f1)=@_;
    my @Fld;
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    print qq[Getting the various entries for file name: ];
    print color('bold green');print qq[$f1\n];print color('reset');
    print qq[We are now storing the details from file: $f1                  \n];

    open($ffoo,'<', $f1) or die "$die_open_string $f1";

    $_= <$ffoo>;
    chop $_[0];
    @Fld = split(' ', $_);
    $orga = $Fld[0];

    close $ffoo;

    print qq[In get_target_oragniser_file orga: ];print color('bold green');
    print qq[$orga\n];print color('reset');

    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    return($orga);
}
################################################################################
#subroutine to extract the prices from the old TickData files.                 #
#Frederic D.R. Bonnet, Date: 31st of Jul. 2015.                                #
################################################################################
sub get_target_ShiftVal_file
{
    ($f1)=@_;
    my @Fld;
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    print qq[Getting the various entries from the file name: ];
    print color('bold green');print qq[$f1\n];print color('reset');
    print qq[We are now storing the details from file: $f1                  \n];
    open($fshf,'<', $f1) or die "$die_open_string $f1";
    my $end = 0;
    my $irows = 0;
    while ( $end == 0 )
    {
	$_= <$fshf>;
	chop $_[0];
	@Fld = split('/', $_);
	$targetShiftVal[$irows] = $Fld[1];
	#print qq[$irows $targetShiftVal[$irows]];
	if ( $Fld[0] eq '' ){$end = 1 ;}
	if ( $end == 0 ) {$irows = $irows + 1;}
    }

    $nrows = $irows;

    close $fshf;
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    return($nrows,@targetShiftVal);
}
################################################################################
#subroutine to extract the prices from the old TickData files.                 #
#Frederic D.R. Bonnet, Date: 31st of Jul. 2015.                                #
################################################################################
sub get_target_Frames_file
{
    ($f1)=@_;
    my @Fld;
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    print qq[Getting the various entries from the file name: ];
    print color('bold green');print qq[$f1\n];print color('reset');
    print qq[We are now storing the details from file: $f1                  \n];
    open($ffrm,'<', $f1) or die "$die_open_string $f1";
    my $end = 0;
    my $irows = 0;
    while ( $end == 0 )
    {
	$_= <$ffrm>;
	chop $_[0];
	@Fld = split('/', $_);
	$targetFrames[$irows] = $Fld[1];
	#print qq[$irows $targetFrames[$irows]];
	if ( $Fld[0] eq '' ){$end = 1 ;}
	if ( $end == 0 ) {$irows = $irows + 1;}
    }

    $nrows = $irows;

    close $ffrm;
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    return($nrows,@targetFrames);
}
################################################################################
#subroutine to extract the prices from the old TickData files.                 #
#Frederic D.R. Bonnet, Date: 11th of May. 2015.                                #
################################################################################
sub get_target_file            
{
    ($f1)=@_;
    my @Fld;
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    print qq[Getting the various entries from the file name: ];
    print color('bold green');print qq[$f1\n];print color('reset');
    print qq[We are now storing the details from file: $f1                  \n];
    open($ftar,'<', $f1) or die "$die_open_string $f1";
    my $end = 0;
    my $irows = 0;
    while ( $end == 0 )
    {
	$_= <$ftar>;
	chop $_[0];
	@Fld = split('_', $_);
	@Fld[1] = split('.dm', @Fld[1]);
	$targetlocation[$irows] = $Fld[1];
	#print qq[$irows $targetlocation[$irows]\n];
	if ( $Fld[0] eq '' ){$end = 1 ;}
	if ( $end == 0 ) {$irows = $irows + 1;}
    }
    $nrows = $irows;
    close $ftar;
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    return($nrows,@targetlocation);
}
