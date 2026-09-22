#!/usr/bin/perl
use strict;
use warnings;

my ($input, $chunks) = @ARGV;
die "Usage: $0 <input_file> <num_chunks>\n" unless $chunks;

open my $fh, '<', $input or die "Cannot open $input: $!";

my @lines = <$fh>;
close $fh;

my $total = @lines;

my $base = int($total / $chunks);
my $extra = $total % $chunks;

my $offset = 0;

for my $i (1 .. $chunks) {
    my $size = $base + ($i < $extra ? 1 : 0);

    open my $out, '>', sprintf("chunk_%02d.txt", $i)
        or die "Cannot write chunk_$i: $!";

    print $out @lines[$offset .. $offset + $size - 1];

    close $out;

    $offset += $size;
}
