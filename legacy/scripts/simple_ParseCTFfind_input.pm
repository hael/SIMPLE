#!/usr/bin/perl
################################################################################
# Input variables for the compilation driving perl script                      #
#                                                                              #
#       Package simple_ParseCTFfind_input                                      #
#                                                                              #
# perl script that lets the user input the global variables that defines the   #
# global variables for the entire SIMPLE Application. This is the only place   #
# where the user needs to intervene in the compilation step of the library.    #
################################################################################
#
package simple_ParseCTFfind_input;  # saved as simple_ParseCTFfind_input.pm
use lib './';
use Exporter;
use warnings;
use Config;

@ISA = ('Exporter');
@EXPORT = qw($NCONFIG $SIMPLE_DATA_PATH $TARGET_FILE $CREATE_DIR $FILE_HANDLER $DELETE_DIR);

#enter the path where the data is hiding and would like to strip
$SIMPLE_DATA_PATH="/media/frederic/LaCie/Data/Simple/Titan_Krios/run2-ribos_3262";

#enter the path where the data is hiding and would like to strip
$NCONFIG=452;             #not used for checking only, real nconfig extracted
$TARGET_FILE = 'target.txt';#Name of output file list

#Clean up  the directory structure: 0 :no 1: yes 
$CREATE_DIR = 0;
#The file handler variables: "distribute" or "regroup" "lookupdata"
$FILE_HANDLER = "distribute";
################################################################################
#!!!!Warning if set to 1 it will delete all the directories and its content!!!!#
#     if you do not want to delete files but just the directories              #
#                   set $FILE_HANDLER = "regroup";                             #
#   this will copy all the files from directories back to working directory    #
################################################################################
#Clean up  the directory structure: 0 :no 1: yes 
$DELETE_DIR = 0;
#regroup can be used with DELETE_DIR

1;

