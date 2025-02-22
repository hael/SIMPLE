#!/usr/bin/perl
use strict;
use warnings;
my $file = $ARGV[0] or die "Need to get CSV file on the command line\n";
open(my $data, '<', $file) or die "Could not open '$file' $!\n";
while (my $line = <$data>) {
  chomp $line;
  my @fields = split "," , $line;
  print $fields[0], $fields[2], $fields[9], $fields[11], $fields[57], "\n";
}
