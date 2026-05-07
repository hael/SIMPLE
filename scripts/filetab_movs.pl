#!/usr/bin/env perl
# Create movies.txt file from movie files in a directory (MRC/EER/TIFF)
use strict;
use warnings;
use Cwd 'abs_path';

if (@ARGV > 1) {
    die "Usage: $0 [directory]\n";
}

my $input_dir = @ARGV ? $ARGV[0] : '.';
die "Input is not a directory: $input_dir\n" unless -d $input_dir;

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

open(my $fh, '>', $outfile) or die "Cannot open $outfile: $!";
for my $file (@selected) {
    my $abs = abs_path("$base_dir/$file");
    print $fh "$abs\n";
}
close($fh);

print "Searched directory: $base_dir\n";
print "Using $kind entries: " . scalar(@selected) . " file(s). Written to $outfile\n";
