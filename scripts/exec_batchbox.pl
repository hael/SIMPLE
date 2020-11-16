#!/usr/bin/perl
@files = glob("*.mrc");
chomp(@files);
$ploc    = 0.2;
$lclim   = 0.08;
$hclim   = 1.0;
$mindist = 120;
foreach $file (@files){
	$boxfile = $file;
	$boxfile =~ s/mrc/box/;
	system("batchboxer input=$file refimg=best.hed auto=$ploc,$hclim,$lclim mindist=$mindist dbout=$boxfile");
}
