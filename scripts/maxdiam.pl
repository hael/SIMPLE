#!/usr/bin/env perl
use strict;
use warnings;

sub usage {
    die "Usage: $0 <pdb_file>\n";
}

my $pdb_file = shift @ARGV;
usage() unless defined $pdb_file;

open(my $fh, '<', $pdb_file) or die "Error: Cannot open $pdb_file: $!\n";

my ($x_min, $y_min, $z_min);
my ($x_max, $y_max, $z_max);
my $found = 0;

while (my $line = <$fh>) {
    next unless $line =~ /^(ATOM  |HETATM)/;
    next if length($line) < 54;

    my $x = substr($line, 30, 8);
    my $y = substr($line, 38, 8);
    my $z = substr($line, 46, 8);

    $x =~ s/^\s+|\s+$//g;
    $y =~ s/^\s+|\s+$//g;
    $z =~ s/^\s+|\s+$//g;

    next unless length($x) && length($y) && length($z);
    next unless $x =~ /^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$/;
    next unless $y =~ /^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$/;
    next unless $z =~ /^[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?$/;

    $x += 0.0;
    $y += 0.0;
    $z += 0.0;

    if (!$found) {
        ($x_min, $y_min, $z_min) = ($x, $y, $z);
        ($x_max, $y_max, $z_max) = ($x, $y, $z);
        $found = 1;
    } else {
        $x_min = $x if $x < $x_min;
        $y_min = $y if $y < $y_min;
        $z_min = $z if $z < $z_min;

        $x_max = $x if $x > $x_max;
        $y_max = $y if $y > $y_max;
        $z_max = $z if $z > $z_max;
    }
}

close($fh);

if (!$found) {
    print "No atomic coordinates found in the PDB file.\n";
    exit 0;
}

my $x_dim = $x_max - $x_min;
my $y_dim = $y_max - $y_min;
my $z_dim = $z_max - $z_min;

my $max_dim = $x_dim;
$max_dim = $y_dim if $y_dim > $max_dim;
$max_dim = $z_dim if $z_dim > $max_dim;

printf "Maximum dimension of the molecule in %s: %.4f A\n", $pdb_file, $max_dim;
