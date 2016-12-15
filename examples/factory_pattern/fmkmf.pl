#!/usr/bin/perl  -w
chomp(@ARGV);
$f90 = $ARGV[0];
$comp_directives = '';
foreach $i (1 .. scalar(@ARGV)-2){
    $comp_directives = $comp_directives.' '.$ARGV[$i];
}
$arg = $ARGV[-1];
$sftag="f90";
@spath=split(/:/,".:..");
@global_outlines=();
@global_objlist=();
@global_modfiles=();
$mainprogfile=$arg;
print "# Main program is $mainprogfile \n" ;
# this subroutine (def below) does most of the work.
process_fsource($mainprogfile); 
# set some makefile . 
print "\n# ------------------Macro-Defs---------------------\n";
print "F90=$f90 \n";
print "\n# -------------------End-macro-Defs---------------------------\n";
# Generate a name for the executable file
$execfile=$mainprogfile;
$execfile=~s/\.${sftag}//;
$execfile=~s|.*/||;
# Generate makefile entry for the Link step
print "\n# Here is the link step \n";
print "$execfile:@global_objlist \n";
print "\t \$(F90) -o $execfile @global_objlist $comp_directives\n";
print "\n# Here are the compile steps\n ";
print STDOUT @global_outlines;
# Add an entry for make clean at the end of the make file.  this
# removes most of the gubbage left around by most of the Unix Fortran
# 90 compilers I have tried. 
print "# This entry allows you to type \" make clean \" to get rid of\n";
print "# all object and module files \n";
print "clean:\n";
print "\trm -f -r f_{files,modd}* *~ *.g90 *.o *.mod *.M *.d V*.inc *.vo \\\n";
print "\tV*.f *.dbg album F.err";
print "\n  \n";
# End of main program 
##############################################
# Here is the subroutine that generates the compile entries in the makefile
# These end up in the global array @global_outlines. The magic part is 
# that this subroutine calls itself recursively.
##############################################
sub process_fsource {  
  my $mainprogfile=$_[0];
  print"# process_fsource called with arg $mainprogfile \n";
  open( MAINPROG, $mainprogfile) or 
    die "Can't find main program file $mainprogfile: $! \n";
  # Read through Fortran source looking for USE statements
  # There should be nothing but whitespace before the USE. Sloppily, 
  # we allow tabs, although the standard (IIRC) does not
  my @modulelist=();
  while ($line=<MAINPROG>) { 
    if ($line =~ /^[ \t]*use[ \t]+(\w+)/i ) { # line matches regexp between / /
      print "# $mainprogfile Uses Module $1\n";
      @modulelist=(@modulelist,$1);
    }
  }
  close(MAINPROG);
  #print "# Full list of modules in $mainprogfile: @modulelist \n";
  print "# Full list of modules in $mainprogfile: @modulelist \n";
  # Find which file each module is in.
 my @modfiles=();
 MODLOOP:foreach $module (@modulelist){
    foreach $directory (@spath){
      # print "# Looking in directory $directory\n";
      opendir( DIRHANDLE, $directory) or die 
	"Can't open directory $directory : $! \n";
      @sourcefiles=grep /\.${sftag}\Z/, sort(readdir(DIRHANDLE));
    foreach $sourcefile (@sourcefiles){
      $pathsourcefile="$directory/$sourcefile";
      #print "\# Checking $pathsourcefile\n";
      open( SOURCEFILE, "$pathsourcefile") or 
	die "Can't find source file $pathsourcefile: $! \n";
      while ($line=<SOURCEFILE>){
	if ($line =~ /^[ \t]*module[ \t]+(\w+)/i ){
	  if($1 =~ /^$module$/i){
	    print "# Uses $module which is in $pathsourcefile\n";
	    @modfiles=(@modfiles,$pathsourcefile);
	    if (grep (/$pathsourcefile/,@global_modfiles )){
	      print "# $pathsourcefile already in list\n";
	    } 
	    else {
	      @global_modfiles=(@global_modfiles,$pathsourcefile);
	      process_fsource($pathsourcefile);
	    }
	    # We found this module -- go on to the next one
	    close (SOURCEFILE);
	    next MODLOOP;	    
	  }
	}
      }
      close( SOURCEFILE );
    }
  }
  # exhausted source files
  print STDERR "Couldn't find source file for module $module\n";
}
# name of file we want to make
$objfile=$mainprogfile;
# replace source file name with .o
$objfile=~s/\.${sftag}/\.o/;
# strip path so object files go in current dir
$objfile=~s|.*/||;
@global_objlist=(@global_objlist,$objfile);
# list of dependencies
@objlist=();
foreach  $mf (@modfiles) { 
  $obj=$mf;
  # replace source file name with .o
  $obj=~s/\.${sftag}/\.o/;
  # strip path so object files go in current dir
  $obj=~s|.*/||;
  @objlist=(@objlist,$obj);
}
@global_outlines=(@global_outlines,"\n$objfile:$mainprogfile @objlist \n");
@global_outlines=(@global_outlines,"\t \$(F90) -c $comp_directives $mainprogfile\n"); # -g
}