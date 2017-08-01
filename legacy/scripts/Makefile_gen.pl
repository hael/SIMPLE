#!/usr/bin/perl
use lib './';
#BEGIN { push @INC, './'}
use warnings;
use strict;
use Cwd qw(getcwd);
use Env;
use Config;

################################################################################
# Calling the input file simple_user_input.pm                                  #
#                                                                              #
#                                                                              #
############################inserted here#######################################

use simple_user_input;

################################################################################

if(scalar(@ARGV) == 0){die "need at one argument: mainprogname.f90\n"};

my @modnames_all;
my @modfiles_glob;
my %ffiles;
my $ending     = "f90";
my $execdir    = getcwd();
my $mainprog   = $ARGV[0];

# set compiler options
my $option;
setCompiling_options();

# goto root directory of the repo
chdir($SIMPLE_PATH);
# recursively walk the directory tree to produce the hash
# $ffiles{simple_prog} = /absolute/path/simple_prog.f90
dir_walk('.',\&ffiles2hash);
# extract all module names
@modnames_all = keys %ffiles;
# substitute the initial "./" with the absolute path
my $subst = $SIMPLE_PATH.'/';
foreach my $i (0 .. $#modnames_all){
    $ffiles{$modnames_all[$i]} =~ s/^\.\//$subst/;  
}
# go back to execution directory
chdir($execdir);
print "\n";
print "$execdir\n";
# recursively process the source files to include all dependencies
process_fsource($mainprog);
#Opening the files for the 
my $filename = 'f90g95_local.csh';
open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";

# create the local linker script
if( $PLATFORM == 0 ){
#for MacOSX
    gen_f90g95_local_OSX10_9_5();
} elsif( $PLATFORM == 1 ){
#for Linux: Ubuntu_14.10
    gen_f90g95_local_LINUX();
}
print qq[\n];
print qq[f90g95_local.csh has been generated\n];
close $fh;
copy_f90g95_local_2prod();
print qq[\n];
print qq[done\n];

################################################################################
# Soubroutines for file handling.                                              #
#                                                                              #
#                                                                              #
################################################################################
#
##Subroutine to generate the Makefile_macros
sub generate_Makefile_macros {

}
sub setCompiling_options {

if( $FCOMPILER eq 'gfortran' ){
    if( $PLATFORM == 0 ){
	if( $DEBUG eq 'yes' ){
	    $option = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore -DBENCH -DMAGMA -fopenmp -g -Og -Wall -fbounds-check -Wextra -pedantic -fbacktrace -fmax-errors=25 ';
	}else{
	    $option = '-fimplicit-none -fall-intrinsics -ffree-form -cpp -fpic -fno-second-underscore -DBENCH -DMAGMA -O3 -fopenmp ';
	}
    } elsif( $PLATFORM == 1 ){
	if( $DEBUG eq 'yes' ){
	    $option = '-ffree-form -cpp -O3 -fPIC -fno-second-underscore -DBENCH -DMAGMA -Og -Wall -fbounds-check -Wextra -pedantic -fdump-parse-tree -ffpe-trap=list -fbacktrace -fdump-core ';
	}else{
	    $option = '-ffree-form -cpp -O3 -fPIC -fno-second-underscore -DBENCH -DMAGMA -fopenmp';
	}
    }
}elsif( $FCOMPILER eq 'ifort' ){
    if( $PLATFORM == 0 ){
	if( $DEBUG eq 'yes' ){
	    $option = '-fimplicit-none -fall-intrinsics -free-form -cpp -fpic -assume no2underscores -DBENCH -DMAGMA -fopenmp -O0 -g -traceback -check all -check bounds -check pointers -debug all ';
	}else{
	    $option = '-fimplicit-none -fall-intrinsics -free-form -cpp -fpic -assume no2underscores -DBENCH -DMAGMA -fast -fopenmp';
	}	
    } elsif( $PLATFORM == 1 ){
	if( $DEBUG eq 'yes' ){
	    $option = '-implicitnone -cpp -O3 -fPIC -nuse -DBENCH -DMAGMA -g';
	}else{
	    $option = '-implicitnone -132 -ffree-form -cpp -O3 -fPIC -nus -DBENCH -DMAGMA -fopenmp';
	}
    }
}else{
    die "Unsupported compiler!";
}


}


sub copy_f90g95_local_2prod {
    print qq[The linking path: $ARGV[0]\n];
    my @Fld = split('/',$ARGV[0]);
    print qq[$Fld[0]\n];

    print qq[has ($#Fld) elements\n];
    print qq[\n];
    my $copy2path = $Fld[0];
    for ( my $iFld = 1 ; $iFld < $#Fld ; $iFld++ ) {
	print qq[$Fld[$iFld]\n];
	$copy2path = $copy2path."/".$Fld[$iFld];  
    }
    print qq[\n];
    print qq[copying file $filename to $copy2path\n];
    system("chmod a+x $filename");
    system("mv $filename $copy2path");
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
	       next if $file eq '.' || $file eq '..';
	       push @results, dir_walk("$top/$file", $filefunc, $dirfunc); # recursion
	    }
	    return $dirfunc ? $dirfunc->($top, @results) : () ;
	}else{
	    return $filefunc ? $filefunc->($top) : () ; 
	}
}

sub file_size{ -s $_[0] }

sub empty{}

sub dir_size{
    my $dir = shift;
    my $total = -s $dir;
    my $n=0;
    for $n (@_){ $total += $n }
    return $total;   
}

sub print_ffile{
    chomp $_[0];
    if( $_[0] =~ /\.f90$/ ){
        print $_[0], "\n";
    }
    return;
}

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
        foreach my $name (@modnames_all){      # looping over all modules available in the library
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

sub gen_f90g95_local_OSX10_9_5{
    my $white = '                                              ';
    print $fh qq[FLNAME=`basename \$1 .f90`\n];
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
    print $fh qq[$white-L \$FFTW_LIB -lfftw3f \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3_threads \\\n];
    print $fh qq[$white-lblas \\\n];
    print $fh qq[$white-llapack \\\n];
    print $fh qq[$white-lstdc++ \\\n];
    if( $USE_GPU == 1 ){
#If we use GPU
	if( $GPU_MODEL == 0 ){
#for CUDA GPU driver invoke
	    print $fh qq[$white-lcublas \\\n];
	    print $fh qq[$white-lcudart \\\n];
	    print $fh qq[$white-L \$CUDADIR/lib \\\n];
	} elsif( $GPU_MODEL == 1 ){
#for OpenCL GPU driver invoke
	    print $fh qq[$white-L \$SHARE_LIB -lGLEW \\\n];
	}
    }
    print $fh qq[$white-lpthread \\\n];
    print $fh qq[$white-lm\n];
    print $fh qq[mv \$FLNAME \$SOURCE_DIR/bin\n];
}

#The gpu generating routine for the Linux OS
sub gen_f90g95_local_LINUX{
    my $white = '                                           ';
    print $fh qq[FLNAME=`basename \$1 .f90`\n];
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
    print $fh qq[$white-L \$FFTW_LIB -lfftw3f \\\n];
    print $fh qq[$white-L \$FFTW_LIB -lfftw3_threads \\\n];
    print $fh qq[$white-lblas \\\n];
    print $fh qq[$white-llapack \\\n];
    print $fh qq[$white-lstdc++ \\\n];
    print $fh qq[$white-lrt \\\n];
    print $fh qq[$white-lpthread \\\n];
    if( $USE_GPU == 1 ){
#If we use GPU
	if( $GPU_MODEL == 0 ){
#for CUDA GPU driver invoke
	    print $fh qq[$white-lcuda \\\n];
	    print $fh qq[$white-lcublas \\\n];
	    print $fh qq[$white-lcudart \\\n];
	    print $fh qq[$white-L \$CUDADIR/lib64 \\\n];
	} elsif( $GPU_MODEL == 1 ){
#for OpenCL GPU driver invoke
	    print $fh qq[$white-lOpenCL \\\n];
	    print $fh qq[$white-L \$SHARE_LIB -lGLEW_x86_64 \\\n];
	}
    }
    print $fh qq[$white-lm\n];
    print $fh qq[\n];
    print $fh qq[mv \$FLNAME \$SOURCE_DIR/bin\n];
    #print $fh qq[return\n];
}
