#!/usr/bin/perl
use lib './';
use warnings;
use strict;
use warnings;
use Term::ANSIColor;
our (@ISA, @EXPORT);
use Exporter;
@ISA = qw(Exporter);
use Cwd qw(getcwd);
use File::Path qw(make_path);
use Tie::File;
################################################################################
# Calling the input module file simple_user_input.pm                           #
#   **** Do not modify this script, this sript automates the compilation ****  #
#   ****                         procedure                               ****  #
#                                                                              #
# Input variables should be modifed in the file:                               #
#                                                                              #
#                          simple_user_input.pm                                #
#                                                                              #
################################################################################
use simple_user_input;
################################################################################
# Declare global variables                                                     #
################################################################################
my @modnames_all;
my @prgnames_all;
my @prgnames_all_short;
my @modfiles_glob;
my %ffiles;
my $ending  = "f90";
my $execdir = getcwd();
my $mainprog;
my @prod_dirs;
my $SIMPLE_SCRIPTS_PATH   = $SIMPLE_PATH.'/scripts';
my $SIMPLE_SRC_PATH       = $SIMPLE_PATH.'/src/simple_main';
my $SIMPLE_PROD_PATH      = $SIMPLE_PATH."/production";
my $SIMPLE_TEST_PROD_PATH = $SIMPLE_PATH."/production/simple_tests";
################################################################################
# Declare local variables                                                      #
################################################################################
my $sptr    = "/";           # set the name for the separator.
my $b_sptr  = "\\";          # set the name for the separator.
my $all     = "*";           # unix commands for all files eg: file*
my $eq      = "=";           # equal sign for file writting
################################################################################
# Set compiler options                                                         #
################################################################################
my $option;
my $mkma_gcc_flags;
my $mkma_gpp_flags;
my $mkma_f90_flags;
my $mkma_f77_flags;
my $mkma_mpif90_flags;
my $mkma_mpif77_flags;

my $option_in;
#compiling options for the Makefile_macros
my $option_mkma_gcc_flags;
my $option_mkma_gpp_flags;
my $option_mkma_f90_flags;
my $option_mkma_f77_flags;
my $option_mkma_mpif90_flags;
my $option_mkma_mpif77_flags;

#optimization level
my $DPLAT; #made global variable from setCompiling_options
my $opti;
my $dbg_lvl_f;
#compiler optimisation varaibales (This goes into the Makefile_macros file) this options are 
#for [linux]
#linking options (this goes into the f90g95_local file)

# test if all required simple directories exist
my$dirtest_pass = 1; check_lib_paths();
# make sure that the object folder is created
make_path("$SIMPLE_PATH/obj/GFORTgpu/");

if( $FCOMPILER =~ /gfortran/ ) {
    if ( $PLATFORM == 0 ) {
    	# debugging options
    	if( $DEBUG eq 'yes' ) {
    	    $dbg_lvl_f = "-g -Og -Wall -fbounds-check";
    	    if( $DEBUG_LEVEL =~ /high/ ) {
    		    # $dbg_lvl_f = $dbg_lvl_f." -Wextra -pedantic -fbacktrace -fmax-errors=25";
                $dbg_lvl_f = $dbg_lvl_f." -Wextra -pedantic -fbacktrace -fdump-parse-tree -fdump-core -frecord-marker=4 -ffpe-summary=all -ffpe-trap=zero,overflow,underflow";
    	    }
    	} else { 
            $dbg_lvl_f = "";
        }
    	# for [MacOSX]
    	$option_in                = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore';
    	# compiling options for the Makefile_macros
    	$option_mkma_gcc_flags    = '-DCUBLAS_GFORTRAN -DADD_';
    	$option_mkma_gpp_flags    = '-DADD_';
    	$option_mkma_f90_flags    = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore';
    	$option_mkma_f77_flags    = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore';
    	$option_mkma_mpif90_flags = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fno-second-underscore';
    	$option_mkma_mpif77_flags = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore';
    } elsif ( $PLATFORM == 1 ) {
    	# debugging options
    	if( $DEBUG eq 'yes' ) {
    	    $dbg_lvl_f = "-g -Og -Wall -fbounds-check";
    	    if( $DEBUG_LEVEL =~ /high/ ) {
    		    # $dbg_lvl_f = $dbg_lvl_f." -Wextra -pedantic -fdump-parse-tree -fbacktrace -fdump-core";
                $dbg_lvl_f = $dbg_lvl_f." -Wextra -pedantic -fbacktrace -fdump-parse-tree -fdump-core -frecord-marker=4 -ffpe-summary=all -ffpe-trap=zero,overflow,underflow";
    		}
	    } else { 
            $dbg_lvl_f = "";
        }
        #for [Linux]
        $option_in                = '-ffree-form -cpp -fPIC -fno-second-underscore';
    	#compiling options for the Makefile_macros
    	$option_mkma_gcc_flags    = '-DCUBLAS_GFORTRAN -DADD_';
    	$option_mkma_gpp_flags    = '-DADD_';
    	$option_mkma_f90_flags    = '-ffree-form -cpp -fno-second-underscore';
    	$option_mkma_f77_flags    = '-ffixed-form -cpp -fPIC -fno-second-underscore';
    	$option_mkma_mpif90_flags = '-ffree-form -cpp -fno-second-underscore';
    	$option_mkma_mpif77_flags = '-ffixed-form -cpp -fPIC -fno-second-underscore';
    }
} elsif( $FCOMPILER =~ /ifort/ ) {
    if( $PLATFORM == 0 ){
    	# debugging options
    	if( $DEBUG eq 'yes' ) {
    	    $dbg_lvl_f = "-g -O0 -check bounds";
    	    if( $DEBUG_LEVEL =~ /high/ ) {
    		    $dbg_lvl_f = $dbg_lvl_f." -traceback -check all -check pointers -debug all";
    	    }
	    } else {
            $dbg_lvl_f = "";
        }
    	# for [MacOSX]
    	$option_in                = '-fimplicit-none -fall-intrinsics -free-form -cpp -fpic -assume no2underscores';
    	# compiling options for the Makefile_macros
    	$option_mkma_gcc_flags    = '-DCUBLAS_GFORTRAN -DADD_';
    	$option_mkma_gpp_flags    = '-DADD_';
    	$option_mkma_f90_flags    = '-fimplicit-none -fall-intrinsics -free-form -cpp -fpic -assume no2underscores';
    	$option_mkma_f77_flags    = '-fimplicit-none -fall-intrinsics -free-form -cpp -fpic -assume no2underscores';
    	$option_mkma_mpif90_flags = '-fimplicit-none -fall-intrinsics -free-form -cpp -fpic -assume no2underscores';
    	$option_mkma_mpif77_flags = '-fimplicit-none -fall-intrinsics -free-form -cpp -fpic -assume no2underscores';
    } elsif( $PLATFORM == 1 ) {
    	#debuging options
    	if( $DEBUG eq 'yes' ) {
    	    $dbg_lvl_f = "-debug -O0";
    	    if( $DEBUG_LEVEL =~ /high/ ) {
    		    $dbg_lvl_f = $dbg_lvl_f."";
    	    }
    	} else { 
            $dbg_lvl_f = "";
        }
    	#for [Linux]
    	$option_in                = '-implicitnone -132 -ffree-form -cpp -fPIC -nus';
    	#compiling options for the Makefile_macros
    	$option_mkma_gcc_flags    = '-DCUBLAS_GFORTRAN -DADD_';
    	$option_mkma_gpp_flags    = '-DADD_';
    	$option_mkma_f90_flags    = '-implicitnone -132 -ffree-form -cpp -fPIC -nus';
    	$option_mkma_f77_flags    = '-implicitnone -132 -ffree-form -cpp -fPIC -nus';
    	$option_mkma_mpif90_flags = '-implicitnone -132 -ffree-form -cpp -fPIC -nus';
    	$option_mkma_mpif77_flags = '-implicitnone -132 -ffree-form -cpp -fPIC -nus';
    }
}
#setting up the options for both Makefile_macros and compile_and_link_local.csh
setCompiling_options();
################################################################################
# Execute SIMPLE code generators                                               #
################################################################################
# goto simple_src directory of the repo
print color('bold blue');
print"Moving to dir: "; print color('bold green'); 
print"$SIMPLE_SRC_PATH\n"; print color('reset');
chdir($SIMPLE_SRC_PATH);
print color('bold blue');
print"Executing simple_args_generator.pl in dir: ";print color('bold green');
print"$SIMPLE_SRC_PATH\n"; print color('reset');
system("$SIMPLE_SCRIPTS_PATH/simple_args_generator.pl");
################################################################################
# Execute SIMPLE code modifiers                                                #
################################################################################
# tie to the simple_args.f90 file
my @simple_args;
tie @simple_args, 'Tie::File', 'simple_args.f90' or die "Could not open file simple_args.f90 $!";
# substitute in the right filename for the list of command-line variables
foreach my $line (@simple_args){
    if( $line =~ /spath\s=\s'(.+)'/ ){
	    $line =~ s/$1/$SIMPLE_PATH/;
    }
}
untie @simple_args;
# goto script directory of the repo
print color('bold blue');
print"Moving to dir: "; print color('bold green'); 
print"$SIMPLE_SCRIPTS_PATH\n"; print color('reset');
chdir($SIMPLE_SCRIPTS_PATH);
print color('bold blue');
print"Tie::distr_simple.pl in dir: ";print color('bold green');
print"$SIMPLE_SCRIPTS_PATH\n"; print color('reset');
my @distr_simple;
tie @distr_simple, 'Tie::File', 'distr_simple.pl' or die "Could not open file distr_simple.pl $!";
$distr_simple[1] =~ s/use lib '.+'/use lib '$SIMPLE_SCRIPTS_PATH'/;
$distr_simple[2] =~ s/use lib '.+'/use lib '$SIMPLE_PATH'/;
untie @distr_simple;
################################################################################
# Generate the Makefile_macros compilation script using the variables          #
# inputed by the user (in simple_user_input.pm)                                #
################################################################################
# goto root directory of the repo
print color('bold blue');
print"Moving to dir: "; print color('bold green'); 
print"$SIMPLE_PATH\n"; print color('reset');
chdir($SIMPLE_PATH);
print color('bold blue');
print"Generating Makefile_macros: ";print color('bold green');
print"$SIMPLE_PATH\n"; print color('reset');
my $filename_mkma = 'Makefile_macros';
open(my $mkma, '>', $filename_mkma) or die "Could not open file '$filename_mkma' $!";
make_Makefile_macros();
close $mkma;
################################################################################
# Generate the Makefile_inject Makefile on script using the variables          #
# inputed by the user (in simple_user_input.pm) CPU test codes                 #
################################################################################
# goto root directory of the repo
my $cpu_test_code_path ="checks/simple_LibTests/cpu_test_code";
$cpu_test_code_path = "$SIMPLE_PATH$sptr$cpu_test_code_path";
print color('bold blue');
print"Moving to dir: "; print color('bold green');
print"$cpu_test_code_path\n"; print color('reset');
chdir($cpu_test_code_path);
print color('bold blue');
print"Generating Makefile_inject: ";print color('bold green');
print"$cpu_test_code_path\n"; print color('reset');
my $f_mkma_injct_cpu = 'Makefile_inject';
open(my $mkma_injct_cpu, '>', $f_mkma_injct_cpu) or die "Could not open file '$f_mkma_injct_cpu' $!";
make_Makefile_inject_cpu($f_mkma_injct_cpu);
close $mkma_injct_cpu;
print color('bold blue');
print"Moving back to dir: "; print color('bold green');
print"$SIMPLE_PATH\n"; print color('reset');
chdir($SIMPLE_PATH);
################################################################################
# Generate the Makefile_inject Makefile on script using the variables          #
# inputed by the user (in simple_user_input.pm) GPU test codes                 #
################################################################################
# goto root directory of the repo
my $gpu_test_code_path ="checks/simple_LibTests/gpu_test_code";
$gpu_test_code_path = "$SIMPLE_PATH$sptr$gpu_test_code_path";
print color('bold blue');
print"Moving to dir: "; print color('bold green');
print"$gpu_test_code_path\n"; print color('reset');
chdir($gpu_test_code_path);
print color('bold blue');
print"Generating Makefile_inject: ";print color('bold green');
print"$gpu_test_code_path\n"; print color('reset');
my $f_mkma_injct_gpu = 'Makefile_inject';
open(my $mkma_injct_gpu, '>', $f_mkma_injct_gpu) or die "Could not open file '$f_mkma_injct_gpu' $!";
make_Makefile_inject_gpu($f_mkma_injct_gpu);
close $mkma_injct_gpu;
print color('bold blue');
print"Moving back to dir: "; print color('bold green');
print"$SIMPLE_PATH\n"; print color('reset');
chdir($SIMPLE_PATH);
################################################################################
# Now compile  the library codes using the Makefile_macros script              #
################################################################################
# first compile library and put all object files in ./obj/GFORTgpu folder
if( $ICOMPILE == 0 ){
    system("make cleanall");
    system("make");
} elsif ( $ICOMPILE == 1 ) {
    system("make clean");
    system("make");
} elsif ( $ICOMPILE == 2 ) {
    # add directive here, at the moment do nothing, just continue.
    system("make");
} elsif ( $ICOMPILE == 3 ) {
    # add directive here, at the moment do nothing, just continue.
}
################################################################################
# Setting up the environment for compile_and_link script generation            #
# for the production codes                                                     #
################################################################################
print "\n";
print "$op_sys, Platform = $PLATFORM\n";
print "Architecture: $architecture\n";
print "\n";
# goto the production folder to compile the production code
print color('bold blue');
print"Moving to dir: "; print color('bold green'); 
print"$SIMPLE_PROD_PATH\n"; print color('reset');
chdir($SIMPLE_PROD_PATH);
print color('bold blue');
print"Generating compile_and_link: ";print color('bold green');
print"$SIMPLE_PROD_PATH\n"; print color('reset');
# recursively walk the directory tree to produce the hash
# $ffiles{simple_prog} = /absolute/path/simple_prog.f90
dir_walk('.',\&ffiles2hash);
# extract all module names
@prgnames_all = keys %ffiles;
# substitute the initial "./" with the absolute path
my $subst = $SIMPLE_PROD_PATH.'/';
print color('bold blue');
print "*********************************************************\n";
print "* These are all programs in                             *\n";
print "* Simple-v2.1                                           *\n";
print "*********************************************************\n";
print color('reset');
foreach my $i (0 .. $#prgnames_all){
    $ffiles{$prgnames_all[$i]} =~ s/^\.\//$subst/;
	my $j = $i+1;
    print "$j: $prgnames_all[$i]\n";
    my @parts = split('/', $ffiles{$prgnames_all[$i]});
    pop @parts;
    $prod_dirs[$i] = join('/', @parts);
    $prgnames_all_short[$i] = $prgnames_all[$i];
	$prgnames_all[$i] = $ffiles{$prgnames_all[$i]};
}
################################################################################
# Compiling  production codes                                                  #
################################################################################
# Generate the database of all modules in the library
print color('bold blue');
print"Moving to dir: "; print color('bold green'); 
print"$SIMPLE_PATH\n"; print color('reset');
chdir($SIMPLE_PATH);
print color('bold blue');
print"Compiling production codes: ";print color('bold green');
print"$SIMPLE_PATH\n"; print color('reset');
# First, forget all about the old %ffiles
undef %ffiles;
# Then, recursively walk the directory tree to produce the hash
# $ffiles{simple_prog} = /absolute/path/simple_prog.f90
dir_walk('.',\&ffiles2hash);
# extract all file names in the library
@modnames_all = keys %ffiles;              
# substitute the initial "./" with the absolute path
$subst = $SIMPLE_PATH.'/';
foreach my $i (0 .. $#modnames_all){
	$ffiles{$modnames_all[$i]} =~ s/^\.\//$subst/;
	my $j = $i+1;
}
# make sure than bin and bin dir are created
make_path("$SIMPLE_PATH/bin/bin_tests/");
# Generate the compile_and_link scripts at the right locations
my $fh;          # filehandle
my $toBcompiled; # file to be compiled
my $filename;    # filename of compile_and_link script
foreach my $j(0 .. $#prod_dirs){
	# goto the production directory
	chdir($SIMPLE_PROD_PATH);
	# this is the program file:
    $toBcompiled = $prgnames_all[$j];
    my$bindir;
    if( $toBcompiled =~ /$SIMPLE_TEST_PROD_PATH/ ){
        $bindir = 'bin/bin_tests';
    }else{
        $bindir = 'bin';
    }
	# Forget all about the old @modfiles_glob
	undef @modfiles_glob;
    # recursively process the source files to include all dependencies
	process_fsource($toBcompiled);
	# goto the program directory
	chdir($prod_dirs[$j]);
    # produce the script
	$filename = './compile_and_link_local.csh';
	open($fh, '>', $filename) or die "Could not open file '$filename' $!";
    my $basename = $toBcompiled;
    $basename =~ s/\.f90//;
	if( $PLATFORM == 0 ){                               # MacOSX
	    gen_compile_and_link_OSX10_9_5($basename, $bindir);
	} elsif( $PLATFORM == 1 ){                          # LINUX
	    gen_compile_and_link_LINUX($basename, $bindir); # for Linux: Ubuntu
	}
	close $fh;
	system("chmod a+x $filename");
	print ">>> COMPILING & LINKING: $prgnames_all_short[$j]\n";
	system("$filename");
}
# produce addonfiles for the bash (.bashrc) and tcsh (.tcshrc) shells
chdir($SIMPLE_PATH);
open(BASH, ">add2.bashrc") or die "Cannot open file add2.bashrc for writing: $!\n";
print BASH "export SIMPLEPATH=$SIMPLE_PATH\n";
print BASH 'export PATH=${SIMPLEPATH}/scripts:${SIMPLEPATH}/bin:$PATH', "\n";
close(BASH);
open(TCSH, ">add2.tcshrc") or die "Cannot open file add2.tcshrc for writing: $!\n";
print TCSH "setenv SIMPLEPATH $SIMPLE_PATH\n";
print TCSH 'set path=(${SIMPLEPATH}/scripts ${SIMPLEPATH}/bin $path)', "\n";
close(TCSH);
my $finishdir = getcwd();
finishing_up_compilation($finishdir);
################################################################################
# Subroutines                                                                  #
################################################################################
#subroutine to finish up script.
sub finishing_up_compilation{
    my ($finishdir_in) = @_;
    print color('bold blue');
    print"SIMPLE library has finished compilation in dir: ";print color('bold green');
    print"$finishdir_in\n"; print color('bold blue');
    print "*******************************************************************************\n";
    print "* Compilation and linkage is complete for Simple-v2.1                         *\n";
    print "* You may run all simple checks  --- List check options                       *\n";
    print "* ";print color('bold yellow');print"> make [OPTIONS] [VAR=VALUE]";
    print color('bold blue');print"   --- ";
    print color('bold yellow');print"> make check_help";
    print color('bold blue');print"                        *\n";
    print "* ";print color('bold yellow');print"            ";
    print color('bold blue');print"                   --- ";
    print color('bold yellow');print"> make check_news";
    print color('bold blue');print"                        *\n";
    print "* ";print color('bold yellow');print"[OPTIONS]:";
    print color('bold blue');print" --- ";
    print color('bold yellow');print"> {check, check_cpu, check_gpu, "; 
    print color('bold blue');print"                             *\n";
    print color('bold yellow'); print"                    ";
    print"bench, check_news, wc, tar}";
    print color('bold blue');print"                               *\n";
    print "* ";print color('bold yellow');print"[VAR]:";
    print color('bold blue');print"     --- ";
    print color('bold yellow');print"> {use_gpu,bench_gpu,fix_gpu,set_gpu,help}";
    print color('bold blue');print"                   *\n";
    print "* ";print color('bold yellow');print"[VALUE]:";
    print color('bold blue');print"   --- ";
    print color('bold yellow');print"> {yes,no},{0,..,MAX_N_GPU}";
    print color('bold blue');print"                                  *\n";
    print "* ";print color('bold yellow');print"[Example]:";
    print color('bold blue');print" --- ";
    print color('bold yellow');print"> make check use_gpu=yes help=no";
    print color('bold blue');print"                             *\n";
    print "* ";print color('bold yellow');print"          ";
    print color('bold blue');print" --- ";
    print color('bold yellow');print"> make bench use_gpu=yes fix_gpu=yes set_gpu=0 bench_gpu=no";
    print color('bold blue');print"  *\n";
    print "* ";print color('bold blue');print"         ";
    print color('bold blue');print"     ";
    print color('bold yellow');print"                                  ";
    print color('bold blue');print"                            *\n";
    print "* ";print color('bold blue');print"Cleaners:";
    print color('bold blue');print"  --- ";
    print color('bold yellow');print"> make {clean,cleanall,clean_check_cpu,"; 
    print color('bold blue');print"                      *\n";
    print color('bold yellow'); print"                         ";
    print"clean_check_gpu}";
    print color('bold blue');print"                                     *\n";
    print "* ";print color('bold blue');print"New Rel.:";
    print color('bold blue');print"  --- ";
    print color('bold yellow');print"> make check_news             ";
    print color('bold blue');print"                               *\n";
    print "* ";print color('bold blue');print"Lne Cntr:";
    print color('bold blue');print"  --- ";
    print color('bold yellow');print"> make wc             ";
    print color('bold blue');print"                                       *\n";
    print "*******************************************************************************\n";
    print color('reset');
}
# Subroutine to check library paths
sub check_lib_paths{
    print color('bold blue');
    print "*********************************************************\n";
    print "* Checking and printing the input directories...        *\n";
    print "*********************************************************\n";
    print color('reset');
    if( -d $SIMPLE_PATH ){
        print color('bold blue');
        print"SIMPLE_PATH          : "; print color('bold red'); 
        print"$SIMPLE_PATH\n"; print color('reset');
    }else{
        print "The SIMPLE root directory: $SIMPLE_PATH does not exist, have you configured simple_user_input.pm correctly?\n";
        $dirtest_pass = 0;
    }
    if( -d $SIMPLE_SRC_PATH ){
        print color('bold blue');
        print"SIMPLE_SRC_PATH      : "; print color('bold red'); 
        print"$SIMPLE_SRC_PATH\n"; print color('reset');
    }else{
        print "The SIMPLE/SRC directory: $SIMPLE_SRC_PATH does not exist, have you configured simple_user_input.pm correctly?\n";
        $dirtest_pass = 0;
    }
    if( -d $SIMPLE_PROD_PATH ){
        print color('bold blue');
        print"SIMPLE_PROD_PATH     : "; print color('bold red'); 
        print"$SIMPLE_PROD_PATH\n"; print color('reset');
    }else{
        print "The directory: $SIMPLE_PROD_PATH does not exist, have you configured simple_user_input.pm correctly?\n";
        $dirtest_pass = 0;
    }

    if( -d $SIMPLE_TEST_PROD_PATH ){
        print color('bold blue');
        print"SIMPLE_TEST_PROD_PATH: "; print color('bold red'); 
        print"$SIMPLE_TEST_PROD_PATH\n"; print color('reset');
    }else{
        print "The directory: $SIMPLE_TEST_PROD_PATH does not exist, have you configured simple_user_input.pm correctly?\n";
        $dirtest_pass = 0;
    }
    if( -d $SIMPLE_SCRIPTS_PATH ){
        print color('bold blue');
        print"SIMPLE_SCRIPTS_PATH  : "; print color('bold red'); 
        print"$SIMPLE_SCRIPTS_PATH\n"; print color('reset');
    }else{
        print "The SIMPLE/SCRIPTS directory: $SIMPLE_SCRIPTS_PATH does not exist, have you configured simple_user_input.pm correctly?\n";
        $dirtest_pass = 0;
    }
    if( $dirtest_pass == 0 ){
        die "The required directory structure is corrupted, exiting...\n";
    }
    print color('bold blue');
    print "*********************************************************\n";
    print color('reset');
}
# Subroutine to recursively process a Fortran source file
sub process_fsource{  
    my $mainprogfile = shift;
    open(MAINPROG, $mainprogfile) or die "Cannot find main program file $mainprogfile: $!\n";
    # Read through Fortran source looking for use statements
    # There should be nothing but whitespace before the use.
    # Sloppily, we allow tabs, although the standard does not.
    my @modulelist;
    while( my $line=<MAINPROG> ){ 
        if( $line =~ /^[\s|\t]*use[\s|\t]+(\w+)/i ){
            push @modulelist, $1;
        }
    }
    close(MAINPROG);
    MODLOOP: foreach my $module (@modulelist){ # looping over all modules in local list
        foreach my $name (@modnames_all){      # looping over all potential modules names in the library
            if( $name =~ /^$module$/i ){       # if matching module found
                # if(not defined($ffiles{$module})){
                #     die "key: $module not defined in %ffiles: \n";
                # }
                if(defined($ffiles{$module})){
                    if( -e $ffiles{$module} ){
                        open(SOURCEFILE, "$ffiles{$module}") or die "Cannot find source file $ffiles{$module}: $!\n";
                        while( my $line=<SOURCEFILE> ){
        	                if( $line =~ /^[\s|\t]*module[\s|\t]+(\w+)/i ){
        	                    if($1 =~ /^$module$/i){ # uses $module (abs path in $ffiles{$module})
        	                        if( grep (/$ffiles{$module}/,@modfiles_glob) ){
                                        # do nothing, the file is already there
        	                        } else {
        	                            push @modfiles_glob, $ffiles{$module};
        	                            process_fsource($ffiles{$module}); # recursion
        	                        }
        	                        close (SOURCEFILE);
        	                        next MODLOOP;	    
        	                    }
        	                }
                        }
                        close(SOURCEFILE);
                    }
                }
            }
        }
        # exhausted source files
        # print STDERR "Couldn't find source file for module $module\n";
    }
}
# Subroutine to set the compilation options
sub setCompiling_options {
    #my $DPLAT; # moved asd global variable
    #setting up the optimization levels for the compilation.
    if ($SET_OPTIMIZATION == 0 ) {
	    $opti = "-O0";
    } elsif ($SET_OPTIMIZATION == 1) {
	    $opti = "-O1";
    } elsif ($SET_OPTIMIZATION == 2) {
	    $opti = "-O2";
    } elsif ($SET_OPTIMIZATION == 3) {
	    if( $FCOMPILER =~ /gfortran/ ) {
	        $opti = "-O3";
	    } elsif( $FCOMPILER =~ /ifort/ ) {
	        $opti = "-fast";
	    }
    } elsif( $DEBUG eq 'yes' ) {
	    $opti = "";
    }

    if ( $PLATFORM == 0 ) {
	    $DPLAT = "-DMACOSX";
    } elsif ( $PLATFORM == 1 ) {
	    $DPLAT = "-DLINUX";
    }

    if( $FCOMPILER =~ /gfortran/ || $FCOMPILER =~ /ifort/ ) {
    	if( $DEBUG eq 'yes' ) {
    	    $option            = $option_in." $dbg_lvl_f $DPLAT $DOPENMP $DBENCH $DCUDA $DOPENCL $DMAGMA $DMKL $DSIMPLE_MPI";
    	    $mkma_gcc_flags    = $option_mkma_gcc_flags." -g $DPLAT $DOPENMP $DBENCH $DCUDA $DOPENCL $DMAGMA";
    	    $mkma_gpp_flags    = $option_mkma_gpp_flags." -g $DPLAT $DOPENMP $DBENCH $DCUDA $DOPENCL $DMAGMA";
    	    $mkma_f90_flags    = $option_mkma_f90_flags." $dbg_lvl_f $DPLAT $DOPENMP $DBENCH $DCUDA $DOPENCL $DMAGMA $DMKL $DSIMPLE_MPI";
    	    $mkma_f77_flags    = $option_mkma_f77_flags." $dbg_lvl_f $DPLAT $DOPENMP $DBENCH $DCUDA $DOPENCL $DMAGMA $DMKL $DSIMPLE_MPI";
    	    $mkma_mpif90_flags = $option_mkma_mpif90_flags." $DPLAT $DBENCH";
    	    $mkma_mpif77_flags = $option_mkma_mpif77_flags." $DPLAT $DBENCH";
    	} else {
    	    $option            = $option_in." $opti $DPLAT $DOPENMP $DBENCH $DCUDA $DOPENCL $DMAGMA $DMKL $DSIMPLE_MPI";
    	    $mkma_gcc_flags    = $option_mkma_gcc_flags." $opti -O $DPLAT $DOPENMP $DBENCH $DCUDA $DOPENCL $DMAGMA";
    	    $mkma_gpp_flags    = $option_mkma_gpp_flags." $opti $DPLAT $DOPENMP $DBENCH $DCUDA $DOPENCL $DMAGMA";
    	    $mkma_f90_flags    = $option_mkma_f90_flags." $opti $DPLAT $DOPENMP $DBENCH $DCUDA $DOPENCL $DMAGMA $DMKL $DSIMPLE_MPI";
    	    $mkma_f77_flags    = $option_mkma_f77_flags." $opti $DPLAT $DOPENMP $DBENCH $DCUDA $DOPENCL $DMAGMA $DMKL $DSIMPLE_MPI";
    	    $mkma_mpif90_flags = $option_mkma_mpif90_flags." $opti $DPLAT $DBENCH";
    	    $mkma_mpif77_flags = $option_mkma_mpif77_flags." $opti $DPLAT $DBENCH";
    	}
    } else {
        die "Unsupported compiler!";
    }
}
# Subroutine to walk through the directories
sub dir_walk{
    my ($top, $filefunc, $dirfunc) = @_;
    my $DIR;
    if( -d $top){
        my $file;
        unless(opendir $DIR, $top){
            warn "Could not open dir $top: $!; skipping.\n";
            return;   
        }
	    my @results;
	    while($file = readdir $DIR){
	       next if $file eq '.' || $file eq '..' || $file =~ /legacy/;
	       push @results, dir_walk("$top/$file", $filefunc, $dirfunc); # recursion
	    }
	    return $dirfunc ? $dirfunc->($top, @results) : () ;
	}else{
	    return $filefunc ? $filefunc->($top) : () ; 
	}
}
# Subroutine to store fortran files in the global %ffiles hash
sub ffiles2hash{
    chomp $_[0];
    if( $_[0] =~ /\.f90$/ ){
        my @arr = split('/', $_[0]);
        chomp(@arr);
        $arr[-1] =~ s/\.f90//;
        $ffiles{$arr[-1]} = $_[0]; 
    }
    return;
}
# Subroutine to generate the MacOSX compile_and_link compiler script
sub gen_compile_and_link_OSX10_9_5{
    my $fname  = shift;
    my $bindir = shift;
    my $white = '                                              ';
    print $fh qq[FLNAME=$fname\n];
    print $fh qq[\n];
    print $fh qq[GFORTRAN=$FCOMPILER\n];
    print $fh qq[\n];
    print $fh qq[SOURCE_DIR=$SIMPLE_PATH\n];
    print $fh qq[\n];
    print $fh qq[MODDIR=\$SOURCE_DIR/$MODDIR\n];
    print $fh qq[OBJDIR=\$SOURCE_DIR/$OBJDIR\n];
    print $fh qq[\n];
    if( $USE_GPU == 1 ){
	print $fh qq[MAGMADIR=$MAGMADIR\n];
	print $fh qq[MAGMADIR_LIB=\$MAGMADIR/lib\n];
	print $fh qq[\n];
    }
    print $fh qq[CUDADIR=$CUDADIR\n];
    print $fh qq[\n];
    print $fh qq[SHARE_LIB=$OPT_DEV_TOOLS_NVIDIA/$OPENCL_SHARE_LIB\n];
    print $fh qq[\n];
    print $fh qq[OPTION="$option"\n];
    print $fh qq[\n];
    print $fh qq[FFTW_LIB=$FFTW_LIB\n];
    print $fh qq[\n];
    print $fh qq[\$GFORTRAN \$OPTION \$LDFLAGS -I \$OBJDIR -I \$MODDIR -o \$FLNAME \\\n];
    print $fh qq[$white\$FLNAME.f90 \\\n];
    foreach my $name (@modfiles_glob){
        my @arr = split('/', $name);
        $arr[-1] =~ s/$ending/o/;
        print $fh qq[$white\$OBJDIR/$arr[-1] \\\n];
    }
    print $fh qq[$white\$OBJDIR/get_systemQuery_cpu.o \\\n];
    print $fh qq[$white\$OBJDIR/get_systemQuery_cpu_c.o \\\n];
    print $fh qq[$white\$OBJDIR/fortran.o \\\n];
    print $fh qq[$white\$OBJDIR/simple_cudnn_fortran.o \\\n];
    print $fh qq[$white\$OBJDIR/get_fft123D_cpu_c.o \\\n];
    print $fh qq[$white\$OBJDIR/get_fft123D_cpu.o \\\n];
    print $fh qq[$white\$OBJDIR/strlen.o \\\n];
    if( $DCUDA =~ /DCUDA/ ) {
	print $fh qq[$white\$OBJDIR/get_deviceQuery_gpu.o \\\n];
	print $fh qq[$white\$OBJDIR/get_deviceQuery_gpu_c.o \\\n];
        print $fh qq[$white\$OBJDIR/get_fft123D_gpu_c.o \\\n];
        print $fh qq[$white\$OBJDIR/get_fft123D_gpu.o \\\n];
        #CUDA kernel insertion here        
        print $fh qq[$white\$OBJDIR/simple_polarft_corr_gpu_c.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_polarft_corr_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_corr_Helpr_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_corr_Hadmr_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_corr_Hadmr_N_N.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_corr_Hadmr_F_N.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_corr_Hadmr_P_N.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_corr_Hadmr_X_N.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_gencorrAll_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_gencorrAll_Z_N.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_multi-GPUs_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_multi-GPUs_Z_M.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_krnl-Opti_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_krnl-Opti_O_N.o \\\n];
        print $fh qq[$white\$OBJDIR/carte2D_ftExt-corr_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/carte2D_ftExt-corr_C_N.o \\\n];
        print $fh qq[$white\$OBJDIR/carte2D_ftExt-corr_C_F.o \\\n];
        print $fh qq[$white\$OBJDIR/Global_polarft.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_pfts_Sizes.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_mesh3D.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_mesh1D.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_mesh1DV.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_img_2D_cart_Sizes.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_ErrorHandler.o \\\n];
        print $fh qq[$white\$OBJDIR/matmul_cuda_transfers.o \\\n];
        print $fh qq[$white\$OBJDIR/matmul_cuda_device.o \\\n];
        print $fh qq[$white\$OBJDIR/matvec_cuda_transfers.o \\\n];
        print $fh qq[$white\$OBJDIR/matvec_cuda_device.o \\\n];
        print $fh qq[$white\$OBJDIR/gen_polar_coords_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_math_gpu_c.o \\\n];
        if( $DMAGMA =~ /DMAGMA/ ) {
            print $fh qq[$white\$OBJDIR/magma_invert_gpu-v2.o \\\n];
            print $fh qq[$white\$OBJDIR/magma_zgetri_blocked_gpu-v2.o \\\n];
            print $fh qq[$white\$OBJDIR/magma_dgetri_blocked_gpu-v2.o \\\n];
            print $fh qq[$white\$OBJDIR/magma_get_getri_nb_gpu.o \\\n];
        }
    }
    print $fh qq[$white\$OBJDIR/get_cpu_time_c.o \\\n];
    print $fh qq[$white\$OBJDIR/timestamp.o \\\n];
    print $fh qq[$white\$OBJDIR/timming.o \\\n];
    print $fh qq[$white\$OBJDIR/timming_c.o \\\n];
    print $fh qq[$white\$OBJDIR/polarft_corr_Helpr_gpu_c.o \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3 \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3f \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3f_threads \\\n];
    print $fh qq[$white-lblas \\\n];
    print $fh qq[$white-llapack \\\n];
    print $fh qq[$white-lstdc++ \\\n];
    # If we use GPU
    if( $USE_GPU == 1 ){
    	if( $GPU_MODEL == 0 ){
    	    # for CUDA GPU driver invoke
    	    if( $DCUDA =~ /DCUDA/ ) {
                print $fh qq[$white-lcublas \\\n];
                print $fh qq[$white-lcudart \\\n];
                print $fh qq[$white-lcufft \\\n];
                print $fh qq[$white-lcudnn \\\n];
                print $fh qq[$white-L \$CUDADIR/lib \\\n];
    	    }
	} elsif( $GPU_MODEL == 1 ){
            # for OpenCL GPU driver invoke
	    if( $DOPENCL =~ /DOPENCL/ ) {
    		print $fh qq[$white-lOpenCL \\\n];
    		print $fh qq[$white-L \$SHARE_LIB -lGLEW \\\n];
	    }
	}
    	if( $DMAGMA =~ /DMAGMA/ ) {
    	    print $fh qq[$white-L -lmagma \\\n];
    	    print $fh qq[$white-L \$MAGMADIR_LIB -lmagma\\\n];
    	}
    }
    print $fh qq[$white-lpthread \\\n];
    print $fh qq[$white-lm\n];
    print $fh qq[mv \$FLNAME \$SOURCE_DIR/$bindir\n];
}
# Subroutine to generate the LINUX compile_and_link compiler script
sub gen_compile_and_link_LINUX{
    my $fname  = shift;
    my $bindir = shift;
    my $white = '                                           ';
    print $fh qq[FLNAME=$fname\n];
    print $fh qq[\n];
    print $fh qq[GFORTRAN=$FCOMPILER\n];
    print $fh qq[\n];
    print $fh qq[SOURCE_DIR=$SIMPLE_PATH\n];
    print $fh qq[\n];
    print $fh qq[MODDIR=\$SOURCE_DIR/$MODDIR\n];
    print $fh qq[OBJDIR=\$SOURCE_DIR/$OBJDIR\n];
    print $fh qq[\n];
    if( $USE_GPU == 1 ){
    	print $fh qq[MAGMADIR=$MAGMADIR\n];
    	print $fh qq[MAGMADIR_LIB=\$MAGMADIR/lib\n];
    	print $fh qq[\n];
    }
    print $fh qq[CUDADIR=$CUDADIR\n];
    print $fh qq[\n];
    print $fh qq[SHARE_LIB=$OPT_DEV_TOOLS_NVIDIA/$OPENCL_SHARE_LIB\n];
    print $fh qq[\n];
    print $fh qq[OPTION="$option"\n];
    print $fh qq[\n];
    print $fh qq[FFTW_LIB=$FFTW_LIB\n];
    print $fh qq[\n];
    print $fh qq[\$GFORTRAN \$OPTION \$LDFLAGS -I \$OBJDIR -I \$MODDIR -o \$FLNAME \\\n];
    print $fh qq[$white\$FLNAME.f90 \\\n];
    foreach my $name (@modfiles_glob){
        my @arr = split('/', $name);
        $arr[-1] =~ s/$ending/o/;
        print $fh qq[$white\$OBJDIR/$arr[-1] \\\n];
    }
    print $fh qq[$white\$OBJDIR/get_systemQuery_cpu.o \\\n];
    print $fh qq[$white\$OBJDIR/get_systemQuery_cpu_c.o \\\n];
    print $fh qq[$white\$OBJDIR/fortran.o \\\n];
    print $fh qq[$white\$OBJDIR/simple_cudnn_fortran.o \\\n];
    print $fh qq[$white\$OBJDIR/get_fft123D_cpu_c.o \\\n];
    print $fh qq[$white\$OBJDIR/get_fft123D_cpu.o \\\n];
    print $fh qq[$white\$OBJDIR/strlen.o \\\n];
    if( $DCUDA =~ /DCUDA/ ) {
        print $fh qq[$white\$OBJDIR/get_deviceQuery_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/get_deviceQuery_gpu_c.o \\\n];
        print $fh qq[$white\$OBJDIR/get_fft123D_gpu_c.o \\\n];
        print $fh qq[$white\$OBJDIR/get_fft123D_gpu.o \\\n];
        #CUDA kernel insertion here        
        print $fh qq[$white\$OBJDIR/simple_polarft_corr_gpu_c.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_polarft_corr_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_corr_Helpr_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_corr_Hadmr_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_corr_Hadmr_N_N.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_corr_Hadmr_F_N.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_corr_Hadmr_P_N.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_corr_Hadmr_X_N.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_gencorrAll_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_gencorrAll_Z_N.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_multi-GPUs_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_multi-GPUs_Z_M.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_krnl-Opti_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/polarft_krnl-Opti_O_N.o \\\n];
        print $fh qq[$white\$OBJDIR/carte2D_ftExt-corr_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/carte2D_ftExt-corr_C_N.o \\\n];
        print $fh qq[$white\$OBJDIR/carte2D_ftExt-corr_C_F.o \\\n];
        print $fh qq[$white\$OBJDIR/Global_polarft.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_pfts_Sizes.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_mesh3D.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_mesh1D.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_mesh1DV.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_img_2D_cart_Sizes.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_ErrorHandler.o \\\n];
        print $fh qq[$white\$OBJDIR/matmul_cuda_transfers.o \\\n];
        print $fh qq[$white\$OBJDIR/matmul_cuda_device.o \\\n];
        print $fh qq[$white\$OBJDIR/matvec_cuda_transfers.o \\\n];
        print $fh qq[$white\$OBJDIR/matvec_cuda_device.o \\\n];
        print $fh qq[$white\$OBJDIR/gen_polar_coords_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/simple_math_gpu_c.o \\\n];
        if( $DMAGMA =~ /DMAGMA/ ) {
            print $fh qq[$white\$OBJDIR/magma_invert_gpu-v2.o \\\n];
            print $fh qq[$white\$OBJDIR/magma_zgetri_blocked_gpu-v2.o \\\n];
            print $fh qq[$white\$OBJDIR/magma_dgetri_blocked_gpu-v2.o \\\n];
            print $fh qq[$white\$OBJDIR/magma_get_getri_nb_gpu.o \\\n];
        }
   }
    print $fh qq[$white\$OBJDIR/get_cpu_time_c.o \\\n];
    print $fh qq[$white\$OBJDIR/timestamp.o \\\n];
    print $fh qq[$white\$OBJDIR/timming.o \\\n];
    print $fh qq[$white\$OBJDIR/timming_c.o \\\n];
    print $fh qq[$white\$OBJDIR/polarft_corr_Helpr_gpu_c.o \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3 \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3f \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3f_threads \\\n];
    print $fh qq[$white-lblas \\\n];
    print $fh qq[$white-llapack \\\n];
    print $fh qq[$white-lstdc++ \\\n];
    print $fh qq[$white-lrt \\\n];
    print $fh qq[$white-lpthread \\\n];
    #If we use GPU
    if( $USE_GPU == 1 ){
        if( $GPU_MODEL == 0 ){
    	    #for CUDA GPU driver invoke
    	    if( $DCUDA =~ /DCUDA/ ) {
                print $fh qq[$white-lcuda \\\n];
                print $fh qq[$white-lcublas \\\n];
                print $fh qq[$white-lcudart \\\n];
                print $fh qq[$white-lcufft \\\n];
                print $fh qq[$white-lcudnn \\\n];
                print $fh qq[$white-L \$CUDADIR/lib64 \\\n];
    	    }
        } elsif( $GPU_MODEL == 1 ){
            #for OpenCL GPU driver invoke
            if( $DOPENCL =~ /DOPENCL/ ) {
                print $fh qq[$white-lOpenCL \\\n];
                print $fh qq[$white-L \$SHARE_LIB -lGLEW_x86_64 \\\n];
            }
        }
        if( $DMAGMA =~ /DMAGMA/ ) {
            print $fh qq[$white-L -lmagma \\\n];
            print $fh qq[$white-L \$MAGMADIR_LIB -lmagma \\\n];
        }
    }
    print $fh qq[$white-lm\n];
    print $fh qq[\n];
    print $fh qq[mv \$FLNAME \$SOURCE_DIR/$bindir\n];
}
#subroutine to make the Makefile_inject in cpu_test_code
# file to be cat into the Make CPU test
sub make_Makefile_inject_cpu{
    my ($f_mkma_injct_cpu_in) = @_;
    print $mkma_injct_cpu qq[##############################################\n];
    print $mkma_injct_cpu qq[# settings\n];
    print $mkma_injct_cpu qq[#\n];
    print $mkma_injct_cpu qq[# COMPILER = gfortran-4.9\n];
    print $mkma_injct_cpu qq[# CFLAGS   = -ffree-form -cpp -O3 -fno-second-underscore -fopenmp -DBENCH\n];
    print $mkma_injct_cpu qq[# LINK     = -lblas -llapack -lstdc++ -lpthread -lm\n];
    print $mkma_injct_cpu qq[#\n];
    print $mkma_injct_cpu qq[# .SUFFIXES: .f90 .o\n];
    print $mkma_injct_cpu qq[# .f90.o:\n];
    print $mkma_injct_cpu qq[#        \$(COMPILER) \$(CFLAGS) -c \$<\n];
    print $mkma_injct_cpu qq[#\n];
    print $mkma_injct_cpu qq[##############################################\n];
    print $mkma_injct_cpu qq[# Massive settings\n];
    print $mkma_injct_cpu qq[# module load fftw/3.3.4-gcc ;\n];
    print $mkma_injct_cpu qq[# module load gcc/4.9.1 ;\n];
    print $mkma_injct_cpu qq[# module load lapack/3.4.2 ;\n];
    print $mkma_injct_cpu qq[# module load openmpi ;\n];
    print $mkma_injct_cpu qq[#\n];
    print $mkma_injct_cpu qq[COMPILER_GFORTRAN = $FCOMPILER\n];
    print $mkma_injct_cpu qq[FFTW_LIB = $FFTW_LIB\n];
    print $mkma_injct_cpu qq[CFLAGS   = -ffree-form -cpp $opti -fno-second-underscore -fopenmp $DPLAT $DBENCH\n];
    print $mkma_injct_cpu qq[LINK     = -L \$(FFTW_LIB) -lfftw3                   \\\n];
    print $mkma_injct_cpu qq[           -L \$(FFTW_LIB) -lfftw3f                  \\\n];
    print $mkma_injct_cpu qq[           -L \$(FFTW_LIB) -lfftw3_threads           \\\n];
    print $mkma_injct_cpu qq[           -lblas -llapack -lstdc++ -lpthread -lm\n];
    system("cat $f_mkma_injct_cpu_in Makefile_input_cpu > Makefile");
    system("rm $f_mkma_injct_cpu_in");
}
#subroutine to make the Makefile_inject in cpu_test_code
# file to be cat into the Make CPU test
sub make_Makefile_inject_gpu{
    #TODO: create the GPU test_code Makefile
    my ($f_mkma_injct_gpu_in) = @_;
    print $mkma_injct_gpu qq[##############################################\n];
    print $mkma_injct_gpu qq[# settings\n];
    print $mkma_injct_gpu qq[#\n];
    print $mkma_injct_gpu qq[# COMPILER = gfortran-4.9\n];
    print $mkma_injct_gpu qq[# CFLAGS   = -ffree-form -cpp $opti -fno-second-underscore -fopenmp $DPLAT -DBENCH\n];
    print $mkma_injct_gpu qq[# LINK     = -lblas -llapack -lstdc++ -lpthread -lm\n];
    print $mkma_injct_gpu qq[#\n];
    print $mkma_injct_gpu qq[# .SUFFIXES: .f90 .o\n];
    print $mkma_injct_gpu qq[# .f90.o:\n];
    print $mkma_injct_gpu qq[#        \$(COMPILER) \$(CFLAGS) -c \$<\n];
    print $mkma_injct_gpu qq[#\n];
    print $mkma_injct_gpu qq[##############################################\n];
    print $mkma_injct_gpu qq[# Massive settings\n];
    print $mkma_injct_gpu qq[# module load fftw/3.3.4-gcc ;\n];
    print $mkma_injct_gpu qq[# module load gcc/4.9.1 ;\n];
    print $mkma_injct_gpu qq[# module load lapack/3.4.2 ;\n];
    print $mkma_injct_gpu qq[# module load openmpi ;\n];
    print $mkma_injct_gpu qq[#\n];
    print $mkma_injct_gpu qq[COMPILER_GFORTRAN = $FCOMPILER\n];
    print $mkma_injct_gpu qq[\n];
    print $mkma_injct_gpu qq[MAGMADIR     = $MAGMADIR\n];
    print $mkma_injct_gpu qq[MAGMADIR_LIB = $MAGMADIR/lib\n];
    print $mkma_injct_gpu qq[CUDADIR      = $CUDADIR\n];
    print $mkma_injct_gpu qq[SHARE_LIB    = $OPT_DEV_TOOLS_NVIDIA$OPENCL_SHARE_LIB\n];
    print $mkma_injct_gpu qq[FFTW_LIB     = $FFTW_LIB\n];
    print $mkma_injct_gpu qq[\n];
    print $mkma_injct_gpu qq[CFLAGS = -ffree-form -cpp -O3 -fno-second-underscore -fopenmp \\\n];
    print $mkma_injct_gpu qq[         $DPLAT $DBENCH $DCUDA $DMAGMA $DOPENCL\n];
    print $mkma_injct_gpu qq[LINK   = -L \$(FFTW_LIB) -lfftw3                   \\\n];
    print $mkma_injct_gpu qq[         -L \$(FFTW_LIB) -lfftw3f                  \\\n];
    print $mkma_injct_gpu qq[         -L \$(FFTW_LIB) -lfftw3_threads           \\\n];
    print $mkma_injct_gpu qq[         -lblas -llapack -lstdc++ -lpthread       \\\n];
    print $mkma_injct_gpu qq[         -lcuda -lcublas -lcudart -lcufft -lcudnn \\\n];
    print $mkma_injct_gpu qq[         -L \$(CUDADIR)/lib64                      \\\n];
    print $mkma_injct_gpu qq[         -lmagma                                  \\\n];
    print $mkma_injct_gpu qq[         -L \$(MAGMADIR_LIB)                       \\\n];
    print $mkma_injct_gpu qq[         -lm\n];
    
    system("cat $f_mkma_injct_gpu_in Makefile_input_gpu > Makefile");
    system("rm $f_mkma_injct_gpu_in");
}
# Subroutine to generate the Makefile_macros
sub make_Makefile_macros {
    print $mkma qq[# Makefile_macros for the gfortran/ifort compilers.\n];
    print $mkma qq[\n];
    print $mkma qq[####### The project path \#####\n];
    print $mkma qq[\n];
    print $mkma qq[Simple_source=$SIMPLE_PATH\n];
    print $mkma qq[\n];
    print $mkma qq[####### The switches \#####\n];
    print $mkma qq[\n];
    print $mkma qq[DOPENMP = $DOPENMP\n];
    print $mkma qq[DCUDA = $DCUDA\n];
    print $mkma qq[DOPENCL = $DOPENCL\n];
    print $mkma qq[DMAGMA = $DMAGMA\n];
    print $mkma qq[DMKL = $DMKL\n];
    print $mkma qq[DSIMPLE_MPI = $DSIMPLE_MPI\n];
    print $mkma qq[DBENCH = $DBENCH\n];
    print $mkma qq[\n];
    print $mkma qq[####### The MPI paths \#####\n];
    print $mkma qq[\n];
    print $mkma qq[MPIDIR=$MPI_DIR\n];
    print $mkma qq[MPIDIR_INCLUDE=$MPI_DIR_INCLUDE\n];
    print $mkma qq[\n];
    print $mkma qq[####### The compilers paths \#####\n];
    print $mkma qq[\n];
    print $mkma qq[CCCOMP=$CC_COMPILER\n];
    print $mkma qq[GCC=$GCC_COMPILER\n];
    print $mkma qq[MPIFORTRAN=$MPI_F_COMPILER\n];
    print $mkma qq[GFORTRAN=$FCOMPILER\n];
    print $mkma qq[\n];
    print $mkma qq[####### The CUDA and MAGMA paths \#####\n];
    print $mkma qq[\n];
    print $mkma qq[CUDADIR=$CUDADIR\n];
    print $mkma qq[\n];
    print $mkma qq[MAGMADIR=$MAGMADIR\n];
    print $mkma qq[MAGMADIR_CONTROL=\$(MAGMADIR)/control\n];
    print $mkma qq[MAGMADIR_INCLUDE=\$(MAGMADIR)/include\n];
    print $mkma qq[\n];
    print $mkma qq[####### The OpenCL paths \#####\n];
    print $mkma qq[\n];
    print $mkma qq[OPT_DEV_TOOLS_NVIDIA = $OPT_DEV_TOOLS_NVIDIA\n];
    print $mkma qq[\n];
    print $mkma qq[####### Set MAGMA-ADDONS={Magma version} to compile Magma addons. \#####\n];
    print $mkma qq[\n];
    print $mkma qq[GPU_MODEL=0.0     \#0: tesla arch, 1: fermi arch \n];
    print $mkma qq[\n];
    print $mkma qq[####### Modules and objects directories. \#####\n];
    print $mkma qq[\n];
    print $mkma qq[OBJDIR=$OBJDIR\n];
    print $mkma qq[MODDIR=$MODDIR\n];
    print $mkma qq[\n];
    print $mkma qq[##############\n];
    print $mkma qq[# C compiler.#\n];
    print $mkma qq[##############\n];
    print $mkma qq[\n];
    print $mkma qq[CC=\$(CCCOMP)\n];
    print $mkma qq[CFLAGS=$mkma_gcc_flags\n];
    print $mkma qq[CCLIB=\$(CC) -c \$(CFLAGS) -I \$(MODDIR)                                            \\\n];
    if( $DOPENCL =~ /DOPENCL/ ) {
    print $mkma qq[                         -I \$(OPT_DEV_TOOLS_NVIDIA)/OpenCL/common/inc/           \\\n];
    print $mkma qq[                         -I \$(OPT_DEV_TOOLS_NVIDIA)/shared/inc/                  \\\n];
    }
    if( $DCUDA =~ /DCUDA/ ) {
    print $mkma qq[                         -I \$(CUDADIR)/include/                                  \\\n];
    print $mkma qq[                         -I \$(CUDADIR)/src                                       \\\n];
    }
    if( $DMAGMA =~ /DMAGMA/ ) {
    print $mkma qq[                         -I \$(MAGMADIR_INCLUDE)                                  \\\n];
    print $mkma qq[                         -I \$(MAGMADIR_CONTROL)                                  \\\n];
    }
    if( $DSIMPLE_MPI =~ /DSIMPLE_MPI/ ){
	print $mkma qq[                         -I \$(MPIDIR_INCLUDE)                                    \\\n];
    }
    print $mkma qq[                         -I \$(Simple_source)/include                             \\\n];
    print $mkma qq[                         -I \$(Simple_source)/include/simple                      \\\n];
    print $mkma qq[                         -I \$(Simple_source)/include/OpenCL                      \\\n];
    print $mkma qq[                         -I \$(Simple_source)/include/mpe                         \\\n];
    print $mkma qq[                         -I $FFTW_INC\n];
    print $mkma qq[\n];
    print $mkma qq[################\n];
    print $mkma qq[# C++ compiler.#\n];
    print $mkma qq[################\n];
    print $mkma qq[\n];
    print $mkma qq[CPP=\$(GCC)\n];
    print $mkma qq[CPPFLAGS=$mkma_gpp_flags\n];
    print $mkma qq[CPPCLIB=\$(CPP) -c \$(CPPFLAGS) -I \$(MODDIR)                                            \\\n];
    if( $DOPENCL =~ /DOPENCL/ ) {
    print $mkma qq[                              -I \$(OPT_DEV_TOOLS_NVIDIA)/OpenCL/common/inc/           \\\n];
    print $mkma qq[                              -I \$(OPT_DEV_TOOLS_NVIDIA)/shared/inc/                  \\\n];
    }
    if( $DCUDA =~ /DCUDA/ ) {
    print $mkma qq[                              -I \$(CUDADIR)/include/                                  \\\n];
    print $mkma qq[                              -I \$(CUDADIR)/src                                       \\\n];
    }
    if( $DMAGMA =~ /DMAGMA/ ) {
    print $mkma qq[                              -I \$(MAGMADIR_INCLUDE)                                  \\\n];
    print $mkma qq[                              -I \$(MAGMADIR_CONTROL)                                  \\\n];
    }
    if( $DSIMPLE_MPI =~ /DSIMPLE_MPI/ ){
	print $mkma qq[                              -I \$(MPIDIR_INCLUDE)                                    \\\n];
    }
    print $mkma qq[                              -I \$(Simple_source)/include                             \\\n];
    print $mkma qq[                              -I \$(Simple_source)/include/simple                      \\\n];
    print $mkma qq[                              -I \$(Simple_source)/include/OpenCL                      \\\n];
    print $mkma qq[                              -I \$(Simple_source)/include/mpe                         \\\n];
    print $mkma qq[                              -I \$(Simple_source)/test_code                           \\\n];
    print $mkma qq[                              -I $FFTW_INC\n];
    if( $DCUDA =~ /DCUDA/ ) {
    print $mkma qq[\n];
    print $mkma qq[#################\n];
    print $mkma qq[# CUDA compiler.#\n];
    print $mkma qq[#################\n];
    print $mkma qq[\n];
    print $mkma qq[#ifeq (\$(GPU_MODEL),0.0)\n];
    print $mkma qq[#	PUSHMEM_GPU= -arch sm_13 -DGPUSHMEM=130 -gencode arch=compute_13,code=compute_13 -gencode arch=compute_10,code=compute_10\n];
    print $mkma qq[#else\n];
    print $mkma qq[#PUSHMEM_GPU= -arch sm_20 -DGPUSHMEM=200 -gencode arch=compute_20,code=compute_20\n];
    print $mkma qq[PUSHMEM_GPU= -arch sm_35 -DGPUSHMEM=350 -gencode arch=compute_35,code=compute_35\n];
    print $mkma qq[#endif\n];
    print $mkma qq[\n];
    print $mkma qq[NVCC=\$(CUDADIR)/bin/nvcc\n];
    if ($DOPENMP =~ /DOPENMP/) {
        print $mkma qq[NVCCFLAGS= --ptxas-options=-v \$(PUSHMEM_GPU) -O3 -DADD_ $DBENCH $DCUDA $DMAGMA -Xcompiler $DOPENMP\n];   
        print $mkma qq[NVCCFLAGS += -D_FORCE_INLINES -ccbin=\$(CXX) -Xcompiler -fPIC \$(COMMON_FLAGS)\n];
    } else {
        print $mkma qq[NVCCFLAGS= --ptxas-options=-v \$(PUSHMEM_GPU) -O3 -DADD_ $DBENCH $DCUDA $DMAGMA\n];
        print $mkma qq[NVCCFLAGS += -D_FORCE_INLINES -ccbin=\$(CXX) -Xcompiler -fPIC \$(COMMON_FLAGS)\n];
    }
    print $mkma qq[NVCCCLIB=\$(NVCC) -c \$(NVCCFLAGS) -I \$(MODDIR)                                \\\n];
    print $mkma qq[                                 -I \$(CUDADIR)/include                       \\\n];
    if( $DMAGMA =~ /DMAGMA/ ) {
    print $mkma qq[                                 -I \$(MAGMADIR_INCLUDE)                      \\\n];
    print $mkma qq[                                 -I \$(MAGMADIR_CONTROL)                      \\\n];
    }
    print $mkma qq[                                 -I \$(Simple_source)/src/simple_gpu/cuda/cub \\\n];
    print $mkma qq[                                 -I \$(Simple_source)/include/cuda            \\\n];
    print $mkma qq[                                 -I \$(Simple_source)/include                 \\\n];
    print $mkma qq[                                 -I \$(Simple_source)/include/simple\n];
    print $mkma qq[\n];
    }
    print $mkma qq[\n];
    print $mkma qq[#######################\n];
    print $mkma qq[# Fortran 77 compiler.#\n];
    print $mkma qq[#######################\n];
    print $mkma qq[\n];
    print $mkma qq[F77C=\$(GFORTRAN)\n];
    print $mkma qq[F77FLAGS=$opti -cpp\n];
    print $mkma qq[F77CLIB=\$(F77C) -c \$(F77FLAGS)\n];
    print $mkma qq[\n];
    print $mkma qq[#######################\n];
    print $mkma qq[# Fortran 90 compiler.#\n];
    print $mkma qq[#######################\n];
    print $mkma qq[\n];
    print $mkma qq[# \$(DLIB)\n];
    print $mkma qq[F90C=\$(GFORTRAN)\n];
    print $mkma qq[F90FLAGS=$mkma_f90_flags\n];
    print $mkma qq[F90FLAGS77=$mkma_f77_flags\n];
    print $mkma qq[F90CLIB=\$(F90C) -c \$(F90FLAGS) -I .                                \\\n];
    print $mkma qq[                               -I \$(MODDIR)                        \\\n];
    print $mkma qq[                               -J \$(MODDIR)                        \\\n];
    if( $DSIMPLE_MPI =~ /DSIMPLE_MPI/ ) {
	print $mkma qq[                               -I \$(MPIDIR_INCLUDE)                \\\n];
    }
    if( $DMAGMA =~ /DMAGMA/ ) {
    print $mkma qq[                               -I \$(MAGMADIR_INCLUDE)              \\\n];
    print $mkma qq[                               -I \$(MAGMADIR_CONTROL)              \\\n];
    }
    print $mkma qq[                               -I \$(Simple_source)/include/simple\n];
    print $mkma qq[\n];
    print $mkma qq[\n];
    print $mkma qq[F90CLIB77=\$(F90C) -c \$(F90FLAGS77) -I .                                \\\n];
    print $mkma qq[                                   -I \$(MODDIR)                        \\\n];
    print $mkma qq[                                   -J \$(MODDIR)                        \\\n];
    if( $DSIMPLE_MPI =~ /DSIMPLE_MPI/ ) {
    print $mkma qq[                                   -I \$(MPIDIR_INCLUDE)                \\\n];
    }
    if( $DMAGMA =~ /DMAGMA/ ) {
    print $mkma qq[                                   -I \$(MAGMADIR_INCLUDE)              \\\n];
    print $mkma qq[                                   -I \$(MAGMADIR_CONTROL)              \\\n];
    }
    print $mkma qq[                                   -I \$(Simple_source)/include/simple\n];
    print $mkma qq[\n];
    print $mkma qq[F90POST=\n];
    print $mkma qq[\n];
    print $mkma qq[#######################\n];
    print $mkma qq[# MPIF90 compiler.    \#\n];
    print $mkma qq[#######################\n];
    print $mkma qq[\n];
    print $mkma qq[MPIF90C=\$(MPIFORTRAN)\n];
    print $mkma qq[MPIF90FLAGS=$mkma_mpif90_flags\n];
    print $mkma qq[MPIF90FLAGS77=$mkma_mpif77_flags\n];
    print $mkma qq[MPIF90CLIB=\$(MPIF90C) -c \$(MPIF90FLAGS) -I .                      \\\n];
    print $mkma qq[                               -I \$(MODDIR)                       \\\n];
    print $mkma qq[                               -J \$(MODDIR)                       \\\n];
    print $mkma qq[                               -I \$(MPIDIR_INCLUDE)               \\\n];
    print $mkma qq[                               -I \$(Simple_source)/include/simple\n];
    print $mkma qq[\n];
    print $mkma qq[MPIF90CLIB77=\$(MPIF90C) -c \$(MPIF90FLAGS77) -I .                     \\\n];
    print $mkma qq[                                  -I \$(MODDIR)                       \\\n];
    print $mkma qq[                                  -J \$(MODDIR)                       \\\n];
    print $mkma qq[                                  -I \$(MPIDIR_INCLUDE)               \\\n];
    print $mkma qq[                                  -I \$(Simple_source)/include/simple\n];
    print $mkma qq[\n];
    print $mkma qq[MPIF90POST=\n];
    print $mkma qq[\n];
}
