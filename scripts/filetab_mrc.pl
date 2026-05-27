#!/usr/bin/env perl
# Create filetab.txt file - List of all mrc files in a directory
use strict;
use warnings;
use Cwd 'abs_path';

if (@ARGV > 1) {
    die "Usage: $0 [directory]\n";
}

my $input_dir = @ARGV ? $ARGV[0] : '.';
die "Input is not a directory: $input_dir\n" unless -d $input_dir;

my $base_dir = abs_path($input_dir);

my $outfile = 'filetab.txt';

opendir(my $dh, $base_dir) or die "Cannot open directory $base_dir: $!\n";
my @entries = readdir($dh);
closedir($dh);

my @mrc_files = sort grep { /\.mrc$/i && -f "$base_dir/$_" } @entries;

open(my $fh, '>', $outfile) or die "Cannot open $outfile: $!";
for my $file (@mrc_files) {
    my $abs = abs_path("$base_dir/$file");
    print $fh "$abs\n";
}
close($fh);

print "Searched directory: $base_dir\n";
print "Found " . scalar(@mrc_files) . " MRC file(s). Written to $outfile\n";
