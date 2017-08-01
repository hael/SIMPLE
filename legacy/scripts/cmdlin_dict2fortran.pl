#!/usr/bin/perl

use warnings;
use strict;

my$single_quote = 
my$cnt = 0;
while(my$line=<>){
    $cnt++;
    chomp($line);
    my@key_val = split('=',$line);
    my$str;
    $str = 'keys('.$cnt.')  = \''.$key_val[0].'\'';
    print $str, "\n";
    $str = 'descr('.$cnt.') = \''.$key_val[1].'\'';
    print $str, "\n";
}