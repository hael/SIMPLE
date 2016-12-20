#!/usr/bin/perl
use warnings;
use strict;
use Tie::File;
if( scalar(@ARGV) < 2 ){
    die "Need filetable as input 1 & nparts as input 2\n";
}
my @files;
tie @files, 'Tie::File', $ARGV[0] or die "Cannot tie to file: $ARGV[0]\n";
my $npart          = $ARGV[1];
my $nlines         = scalar(@files);
my $lines_per_part = $nlines/$npart;
my$leftover        = $nlines-$lines_per_part*$npart;
my $stop           = 0;
my $start          = 0;
for(my $i=1; $i<=$npart; $i++){
    if( $i == $npart ){
        $start = $stop+1;
        $stop  = $nlines;
    }else{
        if( $leftover == 0 ){
            $start = $stop+1;
            $stop  = $start+$lines_per_part-1;
        }else{
            $stop  = $i*($lines_per_part+1);
            $start = $stop-($lines_per_part+1)+1;
            $leftover--;
        }
    }
    my $part_file = 'files_part'.$i.'.txt';
    open(PFH, ">$part_file") or die "Cannot open file: $part_file\n";
    foreach my $j ($start .. $stop){
        print PFH $files[$j-1], "\n";
    }
    close(PFH);
}
untie @files;