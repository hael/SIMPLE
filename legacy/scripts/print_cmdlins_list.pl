#!/usr/bin/perl

use strict;
use warnings;

my%key_vals;
my@bins = glob("simple*");
foreach my$bin (@bins){
    my$cmd = `$bin`;
    $cmd =~ s/\*\*less commonly used\*\*\n//;
    my@values = split('=', $cmd);
    my@keys;
    foreach my$i (0 .. $#values-1){
        if( $values[$i] =~ /\s(\w+)$/ or $values[$i] =~ /\[(\w+)$/ ){
            push(@keys,$1)
        }
    }
    chomp(@values);
    foreach my$i (1 .. $#values){
        if( $i <= $#keys ){
             $values[$i] =~ s/$keys[$i]//;
        }
        $values[$i] =~ s/\<//;
        $values[$i] =~ s/\>//;
        $values[$i] =~ s/\[//;
        $values[$i] =~ s/\]//;
        $values[$i] =~ s/^\s+|\s+$//g
    }
    shift @values;
    foreach my$i (0 .. $#keys){
        $key_vals{$keys[$i]} = $values[$i];
    }
}
foreach my$key (sort (keys(%key_vals))) {
   print "$key=$key_vals{$key}\n";
}

