#!/usr/bin/perl

use warnings;
use strict;

# check so that at least one argument
if(scalar(@ARGV) < 2){die "need at least two argument: oritab=<frealign doc> oriinclude=<yes|no>\n"};

# parse the command line
# $cmd_string is the command line
# %name_value is the hash with input data
my $cmd_string = '';
my %name_value;
foreach my $cmd (@ARGV){
    chomp($cmd);
    my @tmp = split('=', $cmd);
    if( $tmp[0] ne 'prg' ){
        $cmd_string = $cmd_string.$cmd.' ';
    }
    $name_value{$tmp[0]} = $tmp[1];
}
$cmd_string =~ s/\s+$//; # removes trailing whitespace

open(FREALIGN, "<", $name_value{'oritab'}) or die "Cannot open oritab!\n";
while(<FREALIGN>){
    chomp($_);
    my@row = split(/\s+/,$_);
    if( $name_value{'oriinclude'} eq 'yes' ){
        print "e1=$row[2] e2=$row[3] e3=$row[4] x=$row[5] y=$row[6] ";  
    }
    my $dfx = $row[9]/10000.;
    my $dfy = $row[10]/10000.;
    print "dfx=$dfx dfy=$dfy angast=$row[11]\n";
}
close(FREALIGN);