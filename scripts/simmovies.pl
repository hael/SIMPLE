#!/usr/bin/perl

use warnings;
use strict;

my $nmovies = 100;
my $filebody = 'FoilHole';

mkdir 'Movies';
chdir 'Movies';
foreach my $i (1 .. $nmovies){
    my $movie = $filebody.zero_pad_intg($i,$nmovies).'.mrc';
    system("touch $movie");
    sleep(5);
}

sub zero_pad_intg{
    my$intg   = shift;
    my$numlen = shift;
    while( length($intg) < length($numlen) ){
        $intg = '0'.$intg;
    }
    return $intg;
}
