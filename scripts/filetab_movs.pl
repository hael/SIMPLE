#!/usr/bin/env perl
# Create movies.txt file from movie files in a directory (MRC/EER/TIFF)
#
# Usage: filetab_movs.pl [directory] [sample_count]
#        filetab_movs.pl [sample_count]
use strict;
use warnings;
use Cwd 'abs_path';
use List::Util 'shuffle';

if (@ARGV > 2) {
    die "Usage: $0 [directory] [sample_count]\n";
}

my ($input_dir, $sample_count);
if (@ARGV == 0) {
    $input_dir = '.';
} elsif (@ARGV == 1) {
    if (-d $ARGV[0]) {
        $input_dir = $ARGV[0];
    } elsif ($ARGV[0] =~ /^\d+$/) {
        $input_dir = '.';
        $sample_count = $ARGV[0];
    } else {
        $input_dir = $ARGV[0];
    }
} else {
    ($input_dir, $sample_count) = @ARGV;
}

die "Input is not a directory: $input_dir\n" unless -d $input_dir;
if (defined $sample_count) {
    die "Sample count must be a positive integer: $sample_count\n"
        unless $sample_count =~ /^\d+$/ && $sample_count > 0;
}

my $base_dir = abs_path($input_dir);
my $outfile = 'movies.txt';

opendir(my $dh, $base_dir) or die "Cannot open directory $base_dir: $!\n";
my @entries = readdir($dh);
closedir($dh);

my @mrc_files  = sort grep { /\.mrc$/i          && -f "$base_dir/$_" } @entries;
my @eer_files  = sort grep { /\.eer$/i          && -f "$base_dir/$_" } @entries;
my @tiff_files = sort grep { /\.(?:tif|tiff)$/i && -f "$base_dir/$_" } @entries;

my @selected;
my $kind;
if (@mrc_files) {
    @selected = @mrc_files;
    $kind = 'MRC';
} elsif (@eer_files) {
    @selected = @eer_files;
    $kind = 'EER';
} elsif (@tiff_files) {
    @selected = @tiff_files;
    $kind = 'TIFF';
} else {
    die "No .mrc, .eer, .tif, or .tiff files found in $base_dir\n";
}

my $available_count = scalar(@selected);
if (defined $sample_count) {
    die "Requested $sample_count movie(s), but only $available_count $kind file(s) found in $base_dir\n"
        if $sample_count > $available_count;

    my @sampled = (shuffle @selected)[0 .. $sample_count - 1];
    @selected = sort @sampled;
}

open(my $fh, '>', $outfile) or die "Cannot open $outfile: $!";
for my $file (@selected) {
    my $abs = abs_path("$base_dir/$file");
    print $fh "$abs\n";
}
close($fh);

print "Searched directory: $base_dir\n";
print "Using $kind entries: $available_count file(s) found";
print "; randomly selected " . scalar(@selected) . " file(s)" if defined $sample_count;
print ". Written to $outfile\n";
