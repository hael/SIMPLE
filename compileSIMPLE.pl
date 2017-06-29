#!/usr/bin/perl

################################################################################
# Setup the perl env                                                           #
################################################################################
use lib './';
use warnings;
use strict;
use Term::ANSIColor;
our (@ISA, @EXPORT);
use Exporter;
@ISA = qw(Exporter);
use Cwd qw(getcwd);
use File::Path qw(make_path);
use Tie::File;
use simple_user_input;
use File::Copy;

################################################################################
# Declare variables                                                            #
################################################################################
my @modnames_all;
my @prgnames_all;
my @prgnames_all_short;
my @modfiles_glob;
my %ffiles;
my $ending = "(f90|f95|f03|f08)";
my $execdir = getcwd();
my $mainprog;
my @prod_dirs;
my $SIMPLE_PATH           = $execdir;
my $SIMPLE_SCRIPTS_PATH   = $execdir.'/scripts';
my $SIMPLE_SRC_PATH       = $execdir.'/src/simple_main';
my $SIMPLE_PROD_PATH      = $execdir."/production";
my $SIMPLE_TEST_PROD_PATH = $execdir."/production/simple_tests";
my $option;
my $mkma_gcc_flags;
my $mkma_gpp_flags;
my $mkma_f90_flags;
my $option_in;
my $option_mkma_gcc_flags;
my $option_mkma_gpp_flags;
my $option_mkma_f90_flags;

################################################################################
# Set compiler options                                                         #
################################################################################
# optimization level
my $DPLAT;     # platform, from setCompiling_options
my $opti;      # compiler optimisation variable
my $dbg_lvl_f; # debugging directives
# test if all required simple directories exist
check_lib_paths();
# make sure that the object folder is created
make_path("$SIMPLE_PATH/obj/SIMPLEOFILES/");
if( $FCOMPILER =~ /pgfortran/ ){
    # PGI is supposedly platform in-dependent, so no need for diferent options for MacOSX/Linux
    # debugging options
    if( $DEBUG eq 'yes' ) {
        $dbg_lvl_f = "-C -g -Mbounds -Mchkptr -Mnodwarf -Mpgicoff -traceback";
        if( $DEBUG_LEVEL =~ /high/ ) {
            $dbg_lvl_f = $dbg_lvl_f."";
        }
    } else { 
        $dbg_lvl_f = "";
    }
    $option_in = '-fPIC -Mfree -Mpreprocess -DPGI -module obj/SIMPLEOFILES/';
    # compiling options for the Makefile_macros
    $option_mkma_gcc_flags = '';
    $option_mkma_gpp_flags = '-DADD_';
    $option_mkma_f90_flags = '-fPIC -Mfree -Mpreprocess -DPGI -module obj/SIMPLEOFILES/';
}elsif( $FCOMPILER =~ /gfortran/ ){
    if ( $PLATFORM == 0 ) {
    	# debugging options
    	if( $DEBUG eq 'yes' ) {
    	    $dbg_lvl_f = "-g -Og -Wall -fbounds-check -fbacktrace";
    	    if( $DEBUG_LEVEL =~ /high/ ) {
                $dbg_lvl_f = $dbg_lvl_f." -Wextra -pedantic -fdump-parse-tree -fdump-core -frecord-marker=4 -ffpe-summary=all -ffpe-trap=zero,overflow,underflow";
    	    }
    	} else { 
            $dbg_lvl_f = "";
        }
    	# for [MacOSX]
    	$option_in = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore -ffree-line-length-none';
    	# compiling options for the Makefile_macros
    	$option_mkma_gcc_flags = '';
    	$option_mkma_gpp_flags = '-DADD_';
    	$option_mkma_f90_flags = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore -ffree-line-length-none';
    } elsif ( $PLATFORM == 1 ){
    	# debugging options
    	if( $DEBUG eq 'yes' ) {
    	    $dbg_lvl_f = "-g -Og -Wall -fbounds-check";
    	    if( $DEBUG_LEVEL =~ /high/ ) {
                $dbg_lvl_f = $dbg_lvl_f." -Wextra -pedantic -fbacktrace -fdump-parse-tree -fdump-core -frecord-marker=4 -ffpe-summary=all -ffpe-trap=zero,overflow,underflow";
    		}
	    } else { 
            $dbg_lvl_f = "";
        }
        #for [Linux]
      $option_in = '-ffree-form -cpp -fPIC -fno-second-underscore -ffree-line-length-none';
    	#compiling options for the Makefile_macros
    	$option_mkma_gcc_flags = '';
    	$option_mkma_gpp_flags = '-DADD_';
    	$option_mkma_f90_flags = '-ffree-form -cpp -fPIC -fno-second-underscore -ffree-line-length-none';
    }
} elsif( $FCOMPILER =~ /ifort/ ) {
    if( $PLATFORM == 0 ){

        die "ifort compilation on Mac not yet supported";

    } elsif( $PLATFORM == 1 ) {
    	# debuging options
    	if( $DEBUG eq 'yes' ) {
    	    $dbg_lvl_f = "-g -debug -O0 -check bounds -traceback -check uninit -ftrapuv -debug all";
    	    if( $DEBUG_LEVEL =~ /high/ ) {
    		    $dbg_lvl_f = $dbg_lvl_f."";
    	    }
    	} else { 
            $dbg_lvl_f = "";
        }
    	#for [Linux]
    	$option_in = '-implicitnone -80 -cpp -fPIC -DINTEL -module obj/SIMPLEOFILES/';
    	#compiling options for the Makefile_macros
        $option_mkma_gcc_flags = '';
    	$option_mkma_gpp_flags = '-DADD_';
    	$option_mkma_f90_flags = '-implicitnone -80 -cpp -fPIC -DINTEL -module obj/SIMPLEOFILES/';
    }
}
# setting up the options for both Makefile_macros and compile_and_link_local.csh
setCompiling_options();

################################################################################
# Execute SIMPLE code generators/modifiers                                     #
################################################################################
# goto src directory of the repo
print color('bold blue');
print"Moving to dir: "; print color('bold green'); 
print"$SIMPLE_SRC_PATH\n"; print color('reset');
chdir($SIMPLE_SRC_PATH);
print color('bold blue');
print"Executing simple_args_generator.pl in dir: ";print color('bold green');
print"$SIMPLE_SRC_PATH\n"; print color('reset');
system("$SIMPLE_SCRIPTS_PATH/simple_args_generator.pl");
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
# Compile the library codes using the Makefile_macros script                   #
################################################################################
if( $ICOMPILE == 0 ){
    system("make -i cleanall") ;
    system("make") ;
} elsif ( $ICOMPILE == 1 ) {
    system("make clean");
    system("make");
} elsif ( $ICOMPILE == 2 ) {
    system("make");
} elsif ( $ICOMPILE == 3 ) {
    # just continue.
}

################################################################################
# Set up the environment for compile_and_link script generation                #
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
dir_walk('.',\&ffiles2hash);
# extract all module names
@prgnames_all = keys %ffiles;
# substitute the initial "./" with the absolute path
my $subst = $SIMPLE_PROD_PATH.'/';
foreach my $i (0 .. $#prgnames_all){
    $ffiles{$prgnames_all[$i]} =~ s/^\.\//$subst/;
	my $j = $i+1;
    my @parts = split('/', $ffiles{$prgnames_all[$i]});
    pop @parts;
    $prod_dirs[$i] = join('/', @parts);
    $prgnames_all_short[$i] = $prgnames_all[$i];
	$prgnames_all[$i] = $ffiles{$prgnames_all[$i]};
}

################################################################################
# Compile the production codes                                                 #
################################################################################
# Generate the database of all modules in the library
print color('bold blue');
print"Moving to dir: "; print color('bold green'); 
print"$SIMPLE_PATH\n"; print color('reset');
chdir($SIMPLE_PATH);
print color('bold blue');
print"Compiling production codes:\n";print color('reset');
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
# make sure than bin and bin/bin_tests directories are created
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
    $basename =~ s/\.f[09][038]//;
    gen_compile_and_link($basename, $bindir);
	close $fh;
	system("chmod a+x $filename")==0 or die 'You do not have premision to change attributes.  Make sure you have read/write access to the build directory';
	print ">>> COMPILING & LINKING: $prgnames_all_short[$j]\n";
	system("$filename") == 0 or die $filename.' returned an error';
}
# produce addonfiles for the bash (.bashrc) and tcsh (.tcshrc) shells
chdir($SIMPLE_PATH);
open(BASH, ">add2.bashrc") or die "Cannot open file add2.bashrc for writing: $!\n";
print BASH "export SIMPLE_EMAIL=\"my.name\@uni.edu\"\n";
print BASH "export SIMPLE_QSYS=\"local\"\n";
print BASH "export SIMPLE_PATH=$SIMPLE_PATH\n";
print BASH 'export PATH=${SIMPLE_PATH}/scripts:${SIMPLE_PATH}/bin:$PATH', "\n";
close(BASH);
open(TCSH, ">add2.tcshrc") or die "Cannot open file add2.tcshrc for writing: $!\n";
print TCSH "setenv SIMPLE_EMAIL=\"my.name\@uni.edu\"\n";
print TCSH "setenv SIMPLE_QSYS=\"local\"\n";
print TCSH "setenv SIMPLE_PATH $SIMPLE_PATH\n";
print TCSH 'set path=(${SIMPLE_PATH}/scripts ${SIMPLE_PATH}/bin $path)', "\n";
close(TCSH);

# Compile and install gui
compile_gui();

my $finishdir = getcwd();
finish_message($finishdir);

################################################################################
# Subroutines                                                                  #
################################################################################

sub compile_gui{
	my $guidir = $SIMPLE_PATH . "/gui";
	print color('bold blue');
	print "*********************************************************\n";
    print "* Compiling the SIMPLE GUI...                           *\n";
    print "*********************************************************\n";
	print color('bold green');
    make_path($guidir.'/bin');
	if($PLATFORM == 0){
		copy($guidir . "/src/ext/websocketd-mac", $guidir . "/bin/websocketd") or die "Failed: $!";
	}elsif($PLATFORM == 1){
		copy($guidir . "/src/ext/websocketd-linux", $guidir . "/bin/websocketd") or die "Failed: $!";
	}
	chmod 0755, $guidir . "/bin/websocketd";
	print color('bold blue');
	print color('reset');
	system("g++ -DSIMPLE_DIR=" . $SIMPLE_PATH . " " . 
				$guidir . "/src/lodepng.cpp " . 
				$guidir . "/src/base64.cpp " . 
				$guidir . "/src/simple_ws.cpp " . 
				"-ansi -pedantic -Wall -Wextra -O3 -pthread " . 
				"-o " . $guidir . "/bin/simple_ws");	
	system("g++ -DSIMPLE_DIR=" . $SIMPLE_PATH . " " . 
				$guidir . "/src/simple.cpp " . 
				"-ansi -pedantic -Wall -Wextra -O3 -pthread " . 
				"-o " . $SIMPLE_PATH . "/bin/simple");
	print color('bold green');			
	print color('reset');
}

sub finish_message{
    my $finishdir_in = shift;
    print color('bold blue');
    print "Compilation of SIMPLE completed in dir: ";print color('bold green');
    print "$finishdir_in\n"; print color('bold blue');
    print color('reset');
}

sub check_lib_paths{
    my $pass_or_fail = 1;
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
        $pass_or_fail = 0;
    }
    if( -d $SIMPLE_SRC_PATH ){
        print color('bold blue');
        print"SIMPLE_SRC_PATH      : "; print color('bold red'); 
        print"$SIMPLE_SRC_PATH\n"; print color('reset');
    }else{
        print "The SIMPLE/SRC directory: $SIMPLE_SRC_PATH does not exist, have you configured simple_user_input.pm correctly?\n";
        $pass_or_fail = 0;
    }
    if( -d $SIMPLE_PROD_PATH ){
        print color('bold blue');
        print"SIMPLE_PROD_PATH     : "; print color('bold red'); 
        print"$SIMPLE_PROD_PATH\n"; print color('reset');
    }else{
        print "The directory: $SIMPLE_PROD_PATH does not exist, have you configured simple_user_input.pm correctly?\n";
        $pass_or_fail = 0;
    }
    if( -d $SIMPLE_TEST_PROD_PATH ){
        print color('bold blue');
        print"SIMPLE_TEST_PROD_PATH: "; print color('bold red'); 
        print"$SIMPLE_TEST_PROD_PATH\n"; print color('reset');
    }else{
        print "The directory: $SIMPLE_TEST_PROD_PATH does not exist, have you configured simple_user_input.pm correctly?\n";
        $pass_or_fail = 0;
    }
    if( -d $SIMPLE_SCRIPTS_PATH ){
        print color('bold blue');
        print"SIMPLE_SCRIPTS_PATH  : "; print color('bold red'); 
        print"$SIMPLE_SCRIPTS_PATH\n"; print color('reset');
    }else{
        print "The SIMPLE/SCRIPTS directory: $SIMPLE_SCRIPTS_PATH does not exist, have you configured simple_user_input.pm correctly?\n";
        $pass_or_fail = 0;
    }
    if( $pass_or_fail == 0 ){
        die "The required directory structure is corrupted, exiting...\n";
    }
    print color('bold blue');
    print "*********************************************************\n";
    print color('reset');
}

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
    }
}

sub setCompiling_options {

    if ($SET_OPTIMIZATION == 0 ) {
	    $opti = "-O0";
    } elsif ($SET_OPTIMIZATION == 1) {
	    $opti = "-O1";
    } elsif ($SET_OPTIMIZATION == 2) {
	    $opti = "-O2";
    } elsif ($SET_OPTIMIZATION == 3) {
        if( $FCOMPILER =~ /pgfortran/ ){
            # $opti = "-fast -Mipa=fast,inline";
            $opti = "-fast"; 
	    }elsif( $FCOMPILER =~ /gfortran/ ){
	        $opti = "-O3";
	    }elsif( $FCOMPILER =~ /ifort/ ){
            $opti = "-O3 -no-prec-div -static -fp-model fast=2 -xHost";
	    }
    } elsif( $DEBUG eq 'yes' ) {
	    $opti = "";
    }
    if ( $PLATFORM == 0 ) {
	    $DPLAT = "-DMACOSX";
    } elsif ( $PLATFORM == 1 ) {
	    $DPLAT = "-DLINUX";
    }
    if ($DOPENMP ne "" && !($DOPENMP =~ /DOPENMP/) ) {$DOPENMP .= " -DOPENMP ";}
    if( $FCOMPILER =~ /gfortran/ || $FCOMPILER =~ /ifort/ || $FCOMPILER =~ /pgfortran/ ) {
    	if( $DEBUG eq 'yes' ) {
    	    $option            = $option_in." -D_DEBUG $dbg_lvl_f $DPLAT $DOPENMP $DBENCH $DCUDA";
    	    $mkma_gcc_flags    = $option_mkma_gcc_flags." -g $DPLAT $DOPENMP $DBENCH $DCUDA";
    	    $mkma_gpp_flags    = $option_mkma_gpp_flags." -g $DPLAT $DOPENMP $DBENCH $DCUDA";
    	    $mkma_f90_flags    = $option_mkma_f90_flags." -D_DEBUG $dbg_lvl_f $DPLAT $DOPENMP $DBENCH $DCUDA";
    	} else {
    	    $option            = $option_in." $opti $DPLAT $DOPENMP $DBENCH $DCUDA";
    	    $mkma_gcc_flags    = $option_mkma_gcc_flags." $opti -O $DPLAT $DOPENMP $DBENCH $DCUDA";
    	    $mkma_gpp_flags    = $option_mkma_gpp_flags." $opti $DPLAT $DOPENMP $DBENCH $DCUDA";
    	    $mkma_f90_flags    = $option_mkma_f90_flags." $opti $DPLAT $DOPENMP $DBENCH $DCUDA";
    	}
    } else {
        die "Unsupported compiler!";
    }
}

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

sub ffiles2hash{
    chomp $_[0];
    if( $_[0] =~ /\.$ending$/ ){
        my @arr = split('/', $_[0]);
        chomp(@arr);
        $arr[-1] =~ s/\.$ending//;
        $ffiles{$arr[-1]} = $_[0]; 
    }
    return;
}

sub gen_compile_and_link{
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
    print $fh qq[CUDADIR=$CUDADIR\n];
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
    # FFTW
    print $fh qq[$white-L \$FFTW_LIB -lfftw3 \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3f \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3f_threads \\\n];
    # CUDA
    if( $DCUDA =~ /DCUDA/ ) {

        # nothing for now
    }
    print $fh qq[$white-lm\n];
    print $fh qq[\n];
    print $fh qq[mv \$FLNAME \$SOURCE_DIR/$bindir\n];
}

sub make_Makefile_macros {
    print $mkma qq[# Makefile_macros\n];
    print $mkma qq[\n];
    print $mkma qq[####### The project path \#####\n];
    print $mkma qq[\n];
    print $mkma qq[Simple_source=$SIMPLE_PATH\n];
    print $mkma qq[\n];
    print $mkma qq[####### The switches \#####\n];
    print $mkma qq[\n];
    print $mkma qq[DOPENMP = $DOPENMP\n];
    print $mkma qq[DCUDA = $DCUDA\n];
    print $mkma qq[DBENCH = $DBENCH\n];
    print $mkma qq[\n];
    print $mkma qq[####### The compilers paths \#####\n];
    print $mkma qq[\n];
    print $mkma qq[CCCOMP=$CC_COMPILER\n];
    print $mkma qq[GCC=$GCC_COMPILER\n];
    print $mkma qq[GFORTRAN=$FCOMPILER\n];
    print $mkma qq[\n];
    print $mkma qq[####### The CUDA paths \#####\n];
    print $mkma qq[\n];
    print $mkma qq[CUDADIR=$CUDADIR\n];
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
    print $mkma qq[CCLIB=\$(CC) -c \$(CFLAGS) -I \$(MODDIR)  -I\$(Simple_source)/simple_utils        \\\n];
    if( $DCUDA =~ /DCUDA/ ) {
    print $mkma qq[                         -I \$(CUDADIR)/include/                                  \\\n];
    print $mkma qq[                         -I \$(CUDADIR)/src                                       \\\n];
    }
    print $mkma qq[                         -I $FFTW_INC\n];
    print $mkma qq[\n];
    print $mkma qq[##############\n];
    print $mkma qq[# Preprocessor .#\n];
    print $mkma qq[##############\n];
    print $mkma qq[\n];
    print $mkma qq[FPP=cpp-5 \n];
    print $mkma qq[FPPFLAGS= -P \$(DOPENMP) -J \$(MODDIR) -I \$(OBJDIR) -I \$(MODDIR) -I\$(Simple_source)/simple_utils  \\\n];
    if( $DCUDA =~ /DCUDA/ ) {
        print $mkma qq[     \$(DCUDA)           -I \$(CUDADIR)/include/                                  \\\n];
        print $mkma qq[                         -I \$(CUDADIR)/src                                       \\\n];
    }
    print $mkma qq[                         -I $FFTW_INC\n];
    print $mkma qq[\n];
    print $mkma qq[################\n];
    print $mkma qq[# C++ compiler.#\n];
    print $mkma qq[################\n];
    print $mkma qq[\n];
    print $mkma qq[CPP=\$(GCC)\n];
    print $mkma qq[CPPFLAGS=$mkma_gpp_flags\n];
    print $mkma qq[CPPCLIB=\$(CPP) -c \$(CPPFLAGS) -I \$(MODDIR)  -I\$(Simple_source)/simple_utils        \\\n];
    if( $DCUDA =~ /DCUDA/ ) {
    print $mkma qq[                              -I \$(CUDADIR)/include/                                  \\\n];
    print $mkma qq[                              -I \$(CUDADIR)/src                                       \\\n];
    }
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
        print $mkma qq[NVCCFLAGS= --ptxas-options=-v \$(PUSHMEM_GPU) -O3 -DADD_ $DBENCH $DCUDA -Xcompiler $DOPENMP\n];   
        print $mkma qq[NVCCFLAGS += -D_FORCE_INLINES -ccbin=\$(CXX) -Xcompiler -fPIC \$(COMMON_FLAGS)\n];
    } else {
        print $mkma qq[NVCCFLAGS= --ptxas-options=-v \$(PUSHMEM_GPU) -O3 -DADD_ $DBENCH $DCUDA\n];
        print $mkma qq[NVCCFLAGS += -D_FORCE_INLINES -ccbin=\$(CXX) -Xcompiler -fPIC \$(COMMON_FLAGS)\n];
    }
    print $mkma qq[NVCCCLIB=\$(NVCC) -c \$(NVCCFLAGS) -I \$(MODDIR)    -I\$(Simple_source)/simple_utils  \\\n];
    print $mkma qq[                                 -I \$(CUDADIR)/include                       \\\n];
    print $mkma qq[                                 -I \$(Simple_source)/src/simple_gpu/cuda/cub \\\n];
    print $mkma qq[\n];
    }
    print $mkma qq[\n];
    print $mkma qq[#######################\n];
    print $mkma qq[# Fortran 90 compiler.#\n];
    print $mkma qq[#######################\n];
    print $mkma qq[\n];
    print $mkma qq[# \$(DLIB)\n];
    print $mkma qq[F90C=\$(GFORTRAN)\n];
    print $mkma qq[F90FLAGS=$mkma_f90_flags \$(FFLAGS)\n];
    print $mkma qq[F90CLIB=\$(F90C) -c \$(F90FLAGS) -I \$(MODDIR)  -I\$(Simple_source)/simple_utils  \\\n];
    if( $FCOMPILER !~ /pgfortran/ and $FCOMPILER !~ /ifort/ ){
    print $mkma qq[                             -J \$(MODDIR)                      \\\n];
    }
    print $mkma qq[\n];
    print $mkma qq[\n];
    print $mkma qq[F90CLIB77=\$(F90C) -c \$(F90FLAGS77) -I .   -I\$(Simple_source)/simple_utils  \\\n];
    print $mkma qq[                                   -I \$(MODDIR)                        \\\n];
    if( $FCOMPILER !~ /pgfortran/ and $FCOMPILER !~ /ifort/ ){
    print $mkma qq[                                   -J \$(MODDIR)                        \\\n];
    }
    print $mkma qq[\n];
    print $mkma qq[F90POST=\n];
    print $mkma qq[\n];
}
