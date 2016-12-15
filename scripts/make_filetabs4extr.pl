#!/usr/bin/perl
@movfiles = glob("*.mrc");
chomp(@movfiles);
foreach $file (@movfiles){
	$boxfile = $file;
	$boxfile =~ s/mrc/box/;
	if( -e $boxfile ){
		# all good
	}else{
		# need to make a fake boxfile
		open(FAKE, ">$boxfile") or die "Cannot open fake boxfile for writing: $!\n";
		print FAKE "; this is a fake boxfile, just to make the extraction smooth\n";
		close(FAKE);
	}
}
@movfiles = glob("*.mrc");
chomp(@movfiles);
@boxfiles = @movfiles;
foreach $boxfile (@boxfiles){
	$boxfile =~ s/mrc/box/;
}
if( $#movfiles != $#boxfiles ){
	die 'number of entries in movfiles vs. boxfiles arrays do not match'
}
open(BOX, ">boxfiles.txt") or die "Cannot open boxfiles.txt for writing: $!\n";
open(MOV, ">movfiles.txt") or die "Cannot open movfiles.txt for writing: $!\n";
foreach $i (0 .. $#movfiles ){
	print BOX "$boxfiles[$i]\n";
	print MOV "$movfiles[$i]\n";
}
close(BOX);
close(MOV);