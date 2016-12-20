#!/usr/bin/perl
use strict;
use warnings;
if( scalar(@ARGV) < 3 ){
    die "Script needs (1) filetable of individual mrc files, (2) number of frames per movie, (3) filebody of movie\n";
}
my@frames;
open(FRMS, "<$ARGV[0]") or die "Cannot open: $ARGV[0]\n$!\n";
@frames = <FRMS>;
close(FRMS);
chomp(@frames);
my$nframes_per_movie = $ARGV[1];
my$fbody             = $ARGV[2];
my$nframes_tot       = scalar(@frames);
if( $nframes_tot % $nframes_per_movie != 0 ){
    die "Total number of frames must be divisible by number of frames!\n";
}
my$cnt = 0;
for(my $iframe=1; $iframe <= $nframes_tot; $iframe += $nframes_per_movie ){
    # increment movie counter
    $cnt++;
    # make temporary file 4 filetable input to stackops
    my@frame_chunk = splice(@frames,0,$nframes_per_movie);
    open(TMPFTAB, ">temp_filetab.txt");
    foreach my$fff (@frame_chunk) {
        print TMPFTAB $fff, "\n";
    }
    close(TMPFTAB);
    # create the output filename
    my$fname_movie = $fbody.zero_pad_intg($cnt,1000).'.mrc';
    system("simple_stackops filetab=temp_filetab.txt merge=yes outstk=$fname_movie");
}

sub zero_pad_intg{
    my$intg   = shift;
    my$numlen = shift;
    while( length($intg) < length($numlen) ){
        $intg = '0'.$intg;
    }
    return $intg;
}
