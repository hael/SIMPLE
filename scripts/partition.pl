#!/usr/bin/perl
use warnings;
use strict;
if( scalar(@ARGV) < 1 ){
    die "Script needs (1) filebody and (2) chunksize as input to create partitions\n";
}
my @files = glob("$ARGV[0]*");
chomp(@files);
my $chunksz = $ARGV[1];
my $cnt     = 0;
my $part    = 0;
while($cnt < scalar(@files)){
  $part++;
  my $start = $cnt+1;
  $cnt      = $cnt+$chunksz;
  my $end   = $cnt;
  mkdir $part;
  foreach my$i ($start..$end){
    if( defined($files[$i-1]) ){ 
      system("mv $files[$i-1] ./$part");
    }
  }
}
