#!/usr/bin/perl

use warnings;
use strict;

my $nrestarts   = 100;
my $restart_cmd = '';
my $dirbody     = 'ini3D_from_cavgs_run';

foreach my $i (1 .. $nrestarts){
    my $dir = $dirbody.zero_pad_intg($i,$nrestarts);
    mkdir $dir;
    chdir $dir;
    system($restart_cmd);
    chdir "../";
}

sub zero_pad_intg{
    my$intg   = shift;
    my$numlen = shift;
    while( length($intg) < length($numlen) ){
        $intg = '0'.$intg;
    }
    return $intg;
}
