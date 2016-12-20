#!/usr/bin/perl
################################################################################
# Input variables for the compilation driving perl script                      #
#                                                                              #
#       Package simple_user_input                                              #
#                                                                              #
# perl script that lets the user input the global variables that defines the   #
# global variables for the entire SIMPLE Application. This is the only place   #
# where the user needs to intervene in the compilation step of the library.    #
################################################################################
#
package simple_StripData_input;  # saved as simple_StripData_input.pm
use lib './';
use Exporter;
use warnings;
use Config;

@ISA = ('Exporter');
@EXPORT = qw($NFRAMES $NCONFIG $SIMPLE_DATA_PATH $TARGET_FILE $CREATE_DIR $FILE_PRE_HANDLER $FILE_ORGANISE_DIR $DELETE_DIR $UNBLUR $CTFFIND $UNBLUR_DIR $CTFFIND_DIR $FILE_PST_HANDLER $HAVE_BOXFILE $SPH_ABE $AMP_CON $SZE_PWR_SPC $VERBOSE $DATA_ROT);

#enter the path where the data is hiding and would like to strip
$SIMPLE_DATA_PATH='/Volumes/donia/processing/RSC/GridSquare_26623187/Data';

#enter the path where the data is hiding and would like to strip
$NFRAMES = 7;
$NCONFIG = 452;             #not used for checking only, real nconfig extracted
$SPH_ABE = 2.7;             #sphererical aberration
$AMP_CON = 0.07;            #Amplitude constrast
$SZE_PWR_SPC = 512;         #Size of prower spectrum 
$DATA_ROT = -90;            #Rotation angle for the Size of prower spectrum 
$TARGET_FILE = 'target.txt';#Name of output file list

#if FILE_PRE_HANDLER set to distribute then create directories
#Create the directory structure: 0 :no 1: yes 
$CREATE_DIR = 0;  #set to 1 if they not already created otherwise set to 0 
# File Preprocessing handler variables: 
# "lookupdata" : set to get a first estimate of usable data (assums no fls orga)
# "distribute" : set when wanting to organise the file in directory sequence
# "regroup",   : set when files have been distributed and regroup is needed
# "preprocess" : set to preprocess data with unblur and CTFfind4
$FILE_PRE_HANDLER = "preprocess";
#if preprocess directory can be:
# un="unorganise" or="organise" re="reorganise" no="donothing"
$FILE_ORGANISE_DIR = "or";
# The preprocessing parameter 0:no 1: yes
$UNBLUR = 1;
$CTFFIND = 1;
$UNBLUR_DIR = "/opt/Devel_tools/unblur_1.0.2/bin";
$CTFFIND_DIR = "/opt/Devel_tools/ctffind-4.0.16/build";
# File Post-Preprocessing handler variables: 
# "no"          :   set when preprocessing is not done
# "extract_data":   set when extracting the parameters from file*.box
$FILE_PST_HANDLER = "extract_data";
#the boxing variable 0:no 1: yes
#If yes the box files MUST be in a directory called "Boxfile" with the 
# nameming convention: "raw_filename"+"_als_rot".box
# It will not work otherwise
$HAVE_BOXFILE = "0";
#print parsing to screen 0:no 1: yes
$VERBOSE = "1";
################################################################################
#!!!!Warning if set to 1 it will delete all the directories and its content!!!!#
#     if you do not want to delete files but just the directories              #
#                   set $FILE_PRE_HANDLER = "regroup";                         #
#   this will copy all the files from directories back to working directory    #
################################################################################
#Clean up the directory structure: 0 :no 1: yes 
$DELETE_DIR = 0;
#regroup can be used with DELETE_DIR
1;