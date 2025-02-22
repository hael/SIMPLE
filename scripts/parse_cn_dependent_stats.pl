#!/usr/bin/perl
use strict;
use warnings;
my $file = $ARGV[0] or die "Need to get CSV file on the command line\n";
open(my $data, '<', $file) or die "Could not open '$file' $!\n";
my $cnt = 0;
my @cn_natms;
while (my $line = <$data>) {
    $cnt++;
    if( $cnt > 1 ){
        chomp $line; 
        my @fields = split "," , $line;
        $cn_natms[$fields[0]] = $fields[1];
    }
}
my $natms_core;
my $natms_all;
foreach my $i (0 .. $#cn_natms){
    if( exists $cn_natms[$i] ){
	$natms_all += $cn_natms[$i];
	if( $i > 8 ){
	    $natms_core += $cn_natms[$i];
	}
    }
}
print "\%CN\>8 ", 100*($natms_core/$natms_all), "\n";


