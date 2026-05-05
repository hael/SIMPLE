#!/usr/bin/env perl
# Create filetab.txt file - List of all mrc files in the current directory
use strict;
use warnings;
use Cwd 'abs_path';

my $outfile = 'filetab.txt';

opendir(my $dh, '.') or die "Cannot open current directory: $!";
my @mrc_files = sort grep { /\.mrc$/i && -f $_ } readdir($dh);
closedir($dh);

open(my $fh, '>', $outfile) or die "Cannot open $outfile: $!";
for my $file (@mrc_files) {
    my $abs = abs_path($file);
    print $fh "$abs\n";
}
close($fh);

print "Found " . scalar(@mrc_files) . " MRC file(s). Written to $outfile\n";
