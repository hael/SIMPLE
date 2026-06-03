#!/usr/bin/env perl
# Create filetab.txt file - List files in a directory filtered by kind.
#
# Usage: filetab_mrc.pl [directory] [kind]
#
#   kind can be:
#     mrc          - match *.mrc  (default)
#     mrcs         - match *.mrcs
#     all          - match *.mrc and *.mrcs
#     <suffix>     - match *<suffix>  e.g. _int.mrc  _intg.mrc
#
# Examples:
#   filetab_mrc.pl .                  -> all *.mrc in current dir
#   filetab_mrc.pl /data mrcs         -> all *.mrcs in /data
#   filetab_mrc.pl /data _int.mrc     -> all *_int.mrc in /data
#   filetab_mrc.pl /data all          -> all *.mrc and *.mrcs in /data

use strict;
use warnings;
use Cwd 'abs_path';

if (@ARGV > 2) {
    die "Usage: $0 [directory] [mrc|mrcs|all|<suffix>]\n";
}

my $input_dir = (@ARGV >= 1 && -d $ARGV[0]) ? $ARGV[0] : '.';
my $kind      = (@ARGV == 2) ? $ARGV[1]
              : (@ARGV == 1 && !-d $ARGV[0]) ? $ARGV[0]
              : 'mrc';

die "Input is not a directory: $input_dir\n" unless -d $input_dir;

my $base_dir = abs_path($input_dir);
my $outfile  = 'filetab.txt';

# Build the match sub based on kind
my $match;
if    ($kind eq 'mrc')  { $match = sub { /\.mrc$/i  } }
elsif ($kind eq 'mrcs') { $match = sub { /\.mrcs$/i } }
elsif ($kind eq 'all')  { $match = sub { /\.mrcs?$/i } }
else                    { $match = sub { /\Q$kind\E$/i } }  # e.g. _int.mrc

opendir(my $dh, $base_dir) or die "Cannot open directory $base_dir: $!\n";
my @entries = readdir($dh);
closedir($dh);

my @mrc_files = sort grep { $match->() && -f "$base_dir/$_" } @entries;

open(my $fh, '>', $outfile) or die "Cannot open $outfile: $!";
for my $file (@mrc_files) {
    my $abs = abs_path("$base_dir/$file");
    print $fh "$abs\n";
}
close($fh);

print "Searched directory : $base_dir\n";
print "Filter kind        : $kind\n";
print "Found " . scalar(@mrc_files) . " file(s). Written to $outfile\n";
