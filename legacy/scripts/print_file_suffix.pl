#!/usr/bin/perl
use warnings;
use strict;
my @pieces = split(/\./,$ARGV[0]);
my $ext = '.'.$pieces[-1];
print "$ext\n";