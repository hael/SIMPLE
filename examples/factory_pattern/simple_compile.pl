#!/usr/bin/perl
use warnings;
use Env;
#-------------------------------------------
# Parse the command line, %name_value is the hash with input data
#-------------------------------------------
if(scalar(@ARGV) == 0){die "need at least one argument: prg=<simple program name>\n"};
foreach $cmd (@ARGV){
    chomp($cmd);
    @tmp = split('=', $cmd);
    $name_value{$tmp[0]} = $tmp[1];
}
#-------------------------------------------
# Check so that the environment variable pointing to FFTW is set
#-------------------------------------------
if( $ENV{FFTWLIB} ){
}else{
    die "The environment variable FFTWLIB pointing to the FFTW library needs to be set!";   
}

#-------------------------------------------
# Check the command line and set defaults
#-------------------------------------------
if( !defined($name_value{'compiler'}) ){
    $name_value{'compiler'} = 'gfortran';
}
if( !defined($name_value{'optimize'}) ){ 
    $name_value{'optimize'} = 'yes';
}
if( !defined($name_value{'para'}) ){
    $name_value{'para'} = 'yes';
}
if( !defined($name_value{'mac'}) ){
    $name_value{'mac'} = 'yes';
}
if( !defined($name_value{'clean'}) ){
    $name_value{'clean'} = 'no';
}
#-------------------------------------------
# Set compiler directives
#-------------------------------------------
if( $name_value{'optimize'} eq 'yes' ){
    if( $name_value{'compiler'} eq 'gfortran' ){
        $name_value{'directives'} = '-fimplicit-none -fall-intrinsics -ftree-vectorize -O3 -L'.$ENV{FFTWLIB}.' -lfftw3f -lfftw3f_threads ';
        if( $name_value{'para'} eq 'yes' ){
            $name_value{'directives'} = $name_value{'directives'}.'-fopenmp ';
        }
    } elsif( $name_value{'compiler'} eq 'ifort' ){
        $name_value{'directives'} = '-L'.$ENV{FFTWLIB}.' -lfftw3f -lfftw3f_threads ';
        if( $name_value{'mac'} eq 'yes' ){
            $name_value{'directives'} = '-m64 -fast ';
        }
        if( $name_value{'para'} eq 'yes' ){
            $name_value{'directives'} = $name_value{'directives'}.'-openmp ';
        }
    }
} else {
    if( $name_value{'compiler'} eq 'gfortran' ){
        $name_value{'directives'} = '-g -fbounds-check -Wuninitialized -Wimplicit-interface -Wsurprising -Wunderflow -Woverflow -Wunused-parameter -O -ftrapv -fmax-errors=5 -fbacktrace -fimplicit-none -finit-real=nan -L'.$ENV{FFTWLIB}.' -lfftw3f -lfftw3f_threads -ffpe-trap=zero,overflow,invalid ';
        if( $name_value{'para'} eq 'yes' ){
            $name_value{'directives'} = $name_value{'directives'}.'-fopenmp ';
        } 
    } elsif( $name_value{'compiler'} eq 'ifort' ){
        $name_value{'directives'} = '-CB -CU -m64 -g -debug all -warn unused -fp-stack-check -ftrapuv -check all -traceback -heap-arrays -check pointers -check bounds -L'.$ENV{FFTWLIB}.' -lfftw3f -lfftw3f_threads ';
        if( $name_value{'para'} eq 'yes' ){
            $name_value{'directives'} = $name_value{'directives'}.'-openmp ';
        } 
    }
}
#-------------------------------------------
# Extract the program names
#-------------------------------------------
@simple_files = glob("simple_*");
chomp(@simple_files);
$i = 0;
while( $i <= $#simple_files ){
    if( $simple_files[$i] =~ /\.f90$/ or ($simple_files[$i] =~ /\.pl$/ or ($simple_files[$i] =~ /\.txt$/ or $simple_files[$i] =~ /\.log$/))){      
        splice @simple_files, $i, 1;
    } else {
        $i++;
    }
}
#-------------------------------------------
# Execute code generators
#-------------------------------------------
system("./simple_args_generator.pl");

#-------------------------------------------
# Do the compilation
#-------------------------------------------
if( $name_value{'prg'} eq 'all_fast' ){
    if( $name_value{'clean'} eq 'yes' ){ clean('simple_lib') };
    comp('simple_lib',%name_value);
    chdir("./simple_lib");
    @ofiles = glob("*.o");
    chomp(@ofiles);
    $ofiles_str = $ofiles[0];
    foreach $i (1 .. $#ofiles){
        $ofiles_str .= ' '.$ofiles[$i];
    }
    system("ar rcv libsimple.a $ofiles_str");
    chdir("../");
}
$didit = 0;
foreach $file (@simple_files){
    if( $name_value{'prg'} eq 'all_fast' ){
        if( $file ne 'simple_lib' ){ 
            if( $name_value{'clean'} eq 'yes' ){ clean($file) };
            comp_fast($file,%name_value);
            $didit = 1;
        }
    } elsif( $name_value{'prg'} eq 'clean' ){
        clean($file);
        $didit = 1;
    } elsif( $name_value{'prg'} eq 'all' ){
        if( $name_value{'clean'} eq 'yes' ){ clean($file) };
        comp($file);
        $didit = 1;
    }
}
if( $didit == 0 ){
    if( $name_value{'clean'} eq 'yes' ){ clean($name_value{'prg'}) };
    comp($name_value{'prg'});
}
#-------------------------------------------
# Subroutines
#-------------------------------------------
sub comp{
    $prg        = shift;
    $f90file    = $prg.'.f90';
    print ">>> Compiling $prg\n";
    chdir($prg);
    system("../fmkmf.pl $name_value{'compiler'} $name_value{'directives'} $f90file > makefile");
    system("make");
    system("cp $prg ../bin");
    chdir("../");
}
sub comp_fast{
    $prg        = shift;
    $prgo       = $prg.'.o';
    $f90file    = $prg.'.f90';
    print ">>> Compiling $prg\n";
    chdir($prg);
    system("$name_value{'compiler'} -c $f90file -J../simple_lib/ $name_value{'directives'}");
    system("$name_value{'compiler'} -o $prg $prgo ../simple_lib/libsimple.a $name_value{'directives'}");
    system("cp $prg ../bin");
    chdir("../");
}
sub clean{
    $prg = shift;
    print ">>> Cleaning $prg\n";
    chdir($prg);
    system("rm $prg");
    system("make clean");
    system("rm -f -r *~ *.g90 *.o *.mod *.M *.d V*.inc *.vo V*.f *.dbg album F.err libsimple.a");
    chdir("../");
}
