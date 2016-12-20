#!/usr/bin/perl
use warnings;
use strict;

print rounder($ARGV[0]), "\n";

sub rounder{
    my$number     = shift;
    my$number_int = int($number);
    my$diff       = $number-$number_int;
    if( $diff > 0.5 or $diff == 0.5 ){
        return $number_int+1;
    }else{
        return $number_int;
    }
}