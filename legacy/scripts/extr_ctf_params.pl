#!/usr/bin/perl

use warnings;
use strict;

while(<STDIN>){
    chomp($_);
    my@row = split(/\s+/,$_);
    my$str='';
    foreach my $elem (@row){
        if( $elem =~ /dfx/ or $elem =~ /dfy/ or $elem =~ /angast/ or $elem =~ /ctfres/ ){
            $str .= ' '.$elem;
        }
    }
    $str =~ s/^\s+//;
    print $str, "\n";
}