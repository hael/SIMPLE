#!/usr/bin/perl
use lib './';
use warnings;
use strict;
use Cwd qw(getcwd);
use Env;
use Config;
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
my @modfiles_glob;
my %ffiles;
my $ending  = "f90";
my $execdir = getcwd();
my $mainprog;
my @prod_dirs;
my $SIMPLE_SCRIPTS_PATH       = $SIMPLE_PATH.'/scripts';
my $SIMPLE_SRC_PATH           = $SIMPLE_PATH.'/src/simple';
# my $SIMPLE_VARLIST_FNAME      = $SIMPLE_SRC_PATH.'/simple_varlist.txt';
my $SIMPLE_PROD_PATH          = $SIMPLE_PATH."/production/simple";
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
my $opti;
my $dbg_lvl_f;
#compiler optimisation varaibales (This goes into the Makefile_macros file) this options are 
#for [linux]
#linking options (this goes into the f90g95_local file)

if( $FCOMPILER =~ /gfortran/ ) {
    if ( $PLATFORM == 0 ) {
	#debuging options
	if( $DEBUG eq 'yes' ) {
	    $dbg_lvl_f = "-g -Og -Wall -fbounds-check";
	    if( $DEBUG_LEVEL =~ /high/ ) {
		$dbg_lvl_f = $dbg_lvl_f." -Wextra -pedantic -fbacktrace -fmax-errors=25";
	    }
	} else { $dbg_lvl_f = ""; }
	#for [MacOSX]
	$option_in                = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore';
	#compiling options for the Makefile_macros
	$option_mkma_gcc_flags    = '-DCUBLAS_GFORTRAN -DADD_';
	$option_mkma_gpp_flags    = '-DADD_';
	$option_mkma_f90_flags    = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore';
	$option_mkma_f77_flags    = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore';
	$option_mkma_mpif90_flags = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fno-second-underscore';
	$option_mkma_mpif77_flags = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore';
    } elsif ( $PLATFORM == 1 ) {
	#debuging options
	if( $DEBUG eq 'yes' ) {
	               $dbg_lvl_f = "-g -Og -Wall -fbounds-check";
	    if( $DEBUG_LEVEL =~ /high/ ) {
		       $dbg_lvl_f = $dbg_lvl_f." -Wextra -pedantic -fdump-parse-tree -fbacktrace -fdump-core";
		}
	} else { $dbg_lvl_f = ""; }
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
	#debuging options
	if( $DEBUG eq 'yes' ) {
	               $dbg_lvl_f = "-g -O0 -check bounds";
	    if( $DEBUG_LEVEL =~ /high/ ) {
		       $dbg_lvl_f = $dbg_lvl_f." -traceback -check all -check pointers -debug all";
	    }
	} else { $dbg_lvl_f = ""; }
	#for [MacOSX]
	$option_in                = '-fimplicit-none -fall-intrinsics -free-form -cpp -fpic -assume no2underscores';
	#compiling options for the Makefile_macros
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
	} else { $dbg_lvl_f = ""; }
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
chdir($SIMPLE_SRC_PATH);
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
################################################################################
# Generate the Makefile_macros compilation script using the variables          #
# inputted by the user (in simple_user_input.pm)                               #
################################################################################
# goto root directory of the repo
chdir($SIMPLE_PATH);
my $filename_mkma = 'Makefile_macros';
open(my $mkma, '>', $filename_mkma) or die "Could not open file '$filename_mkma' $!";
make_Makefile_macros();
close $mkma;
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
print "$op_sys, Platform = $PLATFORM\n";
print "Architecture: $architecture\n";
print "\n";
# goto the production folder to compile the production code
chdir($SIMPLE_PROD_PATH);
# recursively walk the directory tree to produce the hash
# $ffiles{simple_prog} = /absolute/path/simple_prog.f90
dir_walk('.',\&ffiles2hash);
# extract all module names
@prgnames_all = keys %ffiles;
# substitute the initial "./" with the absolute path
my $subst = $SIMPLE_PROD_PATH.'/';
print "*********************************************************\n";
print "* These are all programs in                             *\n";
print "* Simple-v2.0                                           *\n";
print "*********************************************************\n";
foreach my $i (0 .. $#prgnames_all){
    $ffiles{$prgnames_all[$i]} =~ s/^\.\//$subst/;
	my $j = $i+1;
    print "$j: $prgnames_all[$i]\n";
	#system("ls -1 $ffiles{$prgnames_all[$i]}");
    $prod_dirs[$i]    = $prgnames_all[$i];
	$prgnames_all[$i] = $ffiles{$prgnames_all[$i]};
}
################################################################################
# Compiling  production codes                                                  #
################################################################################
# Generate the database of all modules in the library
chdir($SIMPLE_PATH);
# First, forget all about the old %ffiles
undef %ffiles;
# Then, recursively walk the directory tree to produce the hash
# $ffiles{simple_prog} = /absolute/path/simple_prog.f90
dir_walk('.',\&ffiles2hash);
# extract all file names in the library
@modnames_all = keys %ffiles;              
# substitute the initial "./" with the absolute path
$subst = $SIMPLE_PATH.'/';
#print "*********************************************************\n";
#print "* These are all *.f90 files in                          *\n";
#print "* Simple-v2.0                                           *\n";
#print "*********************************************************\n";
foreach my $i (0 .. $#modnames_all){
	$ffiles{$modnames_all[$i]} =~ s/^\.\//$subst/;
	my $j = $i+1;
	#print "$j: $modnames_all[$i]\n";
	#system("ls -1 $ffiles{$modnames_all[$i]}");
}
# Generate the compile_and_link scripts at the right locations
my $fh;          # filehandle
my $toBcompiled; # file to be compiled
my $filename;    # filename of compile_and_link script
foreach my $j(0 .. $#prod_dirs){
	# goto the production directory
	chdir($SIMPLE_PROD_PATH);
	# this is the program file:
    $toBcompiled = "$SIMPLE_PROD_PATH/$prod_dirs[$j]/$prod_dirs[$j]\.f90";
	# Forget all about the old @modfiles_glob
	undef @modfiles_glob;
    # recursively process the source files to include all dependencies
	process_fsource($toBcompiled);
	# goto the program directory
	chdir($prod_dirs[$j]);
    # produce the script
	$filename = './compile_and_link_local.csh';
	open($fh, '>', $filename) or die "Could not open file '$filename' $!";
	if( $PLATFORM == 0 ){                                  # MacOSX
	    gen_compile_and_link_OSX10_9_5($prod_dirs[$j]);
	} elsif( $PLATFORM == 1 ){                             # LINUX
	    gen_compile_and_link_LINUX($prod_dirs[$j]); # for Linux: Ubuntu_14.10
	}
	close $fh;
	system("chmod a+x $filename");
	print ">>> COMPILING & LINKING: $prod_dirs[$j]\n";
	system("$filename");
}
################################################################################
# Subroutines                                                                  #
################################################################################
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
                open(SOURCEFILE, "$ffiles{$module}") or die "Can't find source file $ffiles{$module}: $!\n";
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
        # exhausted source files
        print STDERR "Couldn't find source file for module $module\n";
    }
}
# Subroutine to set the compilation options
sub setCompiling_options {
    my $DPLAT;
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
	       next if $file eq '.' || $file eq '..';
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
    my $fname = shift;
    my $white = '                                              ';
    #print $fh qq[FLNAME=`basename \$1 .f90`\n];
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
    print $fh qq[\$GFORTRAN \$OPTION -I \$OBJDIR -I \$MODDIR -o \$FLNAME \\\n];
    print $fh qq[$white\$FLNAME.f90 \\\n];
    foreach my $name (@modfiles_glob){
        my @arr = split('/', $name);
        $arr[-1] =~ s/$ending/o/;
        print $fh qq[$white\$OBJDIR/$arr[-1] \\\n];
    }
    print $fh qq[$white\$OBJDIR/get_systemQuery_cpu.o \\\n];
    print $fh qq[$white\$OBJDIR/get_systemQuery_cpu_c.o \\\n];
    print $fh qq[$white\$OBJDIR/fortran.o \\\n];
#    print $fh qq[$white\$OBJDIR/fft123D_cpu.o \\\n];
    print $fh qq[$white\$OBJDIR/get_fft123D_cpu_c.o \\\n];
    print $fh qq[$white\$OBJDIR/get_fft123D_cpu.o \\\n];
    if( $DCUDA =~ /DCUDA/ ) {
	print $fh qq[$white\$OBJDIR/get_deviceQuery_gpu.o \\\n];
	print $fh qq[$white\$OBJDIR/get_deviceQuery_gpu_c.o \\\n];
#        print $fh qq[$white\$OBJDIR/fft123D_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/get_fft123D_gpu_c.o \\\n];
        print $fh qq[$white\$OBJDIR/get_fft123D_gpu.o \\\n];
    }
    print $fh qq[$white\$OBJDIR/get_cpu_time_c.o \\\n];
    print $fh qq[$white\$OBJDIR/timestamp.o \\\n];
    print $fh qq[$white\$OBJDIR/timming_c.o \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3 \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3f \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3_threads \\\n];
    print $fh qq[$white-lblas \\\n];
    print $fh qq[$white-llapack \\\n];
    print $fh qq[$white-lstdc++ \\\n];
    #If we use GPU
    if( $USE_GPU == 1 ){
	if( $GPU_MODEL == 0 ){
	    #for CUDA GPU driver invoke
	    if( $DCUDA =~ /DCUDA/ ) {
		print $fh qq[$white-lcublas \\\n];
		print $fh qq[$white-lcudart \\\n];
		print $fh qq[$white-lcufft \\\n];
		print $fh qq[$white-L \$CUDADIR/lib \\\n];
	    }
	} elsif( $GPU_MODEL == 1 ){
            #for OpenCL GPU driver invoke
	    if( $DOPENCL =~ /DOPENCL/ ) {
		print $fh qq[$white-lOpenCL \\\n];
		print $fh qq[$white-L \$SHARE_LIB -lGLEW \\\n];
	    }
	}
	if( $DMAGMA =~ /DMAGMA/ ) {
	    print $fh qq[$white-L -lmagma \\\n];
	    print $fh qq[$white-L \$MAGMADIR_LIB \\\n];
	}
    }
	print $fh qq[$white-lpthread \\\n];
	print $fh qq[$white-lm\n];
	print $fh qq[mv \$FLNAME \$SOURCE_DIR/bin\n];
}
# Subroutine to generate the LINUX compile_and_link compiler script
sub gen_compile_and_link_LINUX{
    my $fname = shift;
    my $white = '                                           ';
    #print $fh qq[FLNAME=`basename \$1 .f90`\n];
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
    print $fh qq[\$GFORTRAN \$OPTION -I \$OBJDIR -I \$MODDIR -o \$FLNAME \\\n];
    print $fh qq[$white\$FLNAME.f90 \\\n];
    foreach my $name (@modfiles_glob){
        my @arr = split('/', $name);
        $arr[-1] =~ s/$ending/o/;
        print $fh qq[$white\$OBJDIR/$arr[-1] \\\n];
    }
    print $fh qq[$white\$OBJDIR/get_systemQuery_cpu.o \\\n];
    print $fh qq[$white\$OBJDIR/get_systemQuery_cpu_c.o \\\n];
    print $fh qq[$white\$OBJDIR/fortran.o \\\n];
#    print $fh qq[$white\$OBJDIR/fft123D_cpu.o \\\n];
    print $fh qq[$white\$OBJDIR/get_fft123D_cpu_c.o \\\n];
    print $fh qq[$white\$OBJDIR/get_fft123D_cpu.o \\\n];
    if( $DCUDA =~ /DCUDA/ ) {
	print $fh qq[$white\$OBJDIR/get_deviceQuery_gpu.o \\\n];
	print $fh qq[$white\$OBJDIR/get_deviceQuery_gpu_c.o \\\n];
#        print $fh qq[$white\$OBJDIR/fft123D_gpu.o \\\n];
        print $fh qq[$white\$OBJDIR/get_fft123D_gpu_c.o \\\n];
        print $fh qq[$white\$OBJDIR/get_fft123D_gpu.o \\\n];
    }
    print $fh qq[$white\$OBJDIR/get_cpu_time_c.o \\\n];
    print $fh qq[$white\$OBJDIR/timestamp.o \\\n];
    print $fh qq[$white\$OBJDIR/timming_c.o \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3 \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3f \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3_threads \\\n];
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
	    print $fh qq[$white-L \$MAGMADIR_LIB \\\n];
	}
    }
    print $fh qq[$white-lm\n];
    print $fh qq[\n];
    print $fh qq[mv \$FLNAME \$SOURCE_DIR/bin\n];
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
    print $mkma qq[NVCCFLAGS= --ptxas-options=-v \$(PUSHMEM_GPU) -O3 -DADD_ $DBENCH $DCUDA $DMAGMA\n];
    print $mkma qq[NVCCCLIB=\$(NVCC) -c \$(NVCCFLAGS) -I \$(MODDIR)                               \\\n];
    print $mkma qq[                                 -I \$(CUDADIR)/include                      \\\n];
    if( $DMAGMA =~ /DMAGMA/ ) {
    print $mkma qq[                                 -I \$(MAGMADIR_INCLUDE)                     \\\n];
    print $mkma qq[                                 -I \$(MAGMADIR_CONTROL)                     \\\n];
    }
    print $mkma qq[                                 -I \$(Simple_source)/include/cuda           \\\n];
    print $mkma qq[                                 -I \$(Simple_source)/include                \\\n];
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
