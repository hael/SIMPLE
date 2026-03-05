#!/usr/bin/perl
use strict;
use warnings;

# Get seconds from command line argument
my $seconds = shift or die "Usage: $0 <seconds>\n";

my $hours   = int($seconds / 3600);
my $minutes = int(($seconds % 3600) / 60);
my $secs    = $seconds % 60;

printf "%d seconds = %d hours, %d minutes, %d seconds\n",
       $seconds, $hours, $minutes, $secs;
       