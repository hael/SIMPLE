#!/usr/bin/env perl
use strict;
use warnings;
use File::Glob qw(bsd_glob);

# Usage:
#   perl parse_refine3d_bench_combined.pl [glob]
# Example:
#   perl parse_refine3d_bench_combined.pl "REFINE3D_BENCH_ITER*.txt"
#
# Output:
#   timings_combined.csv

my $pattern = shift(@ARGV) // 'REFINE3D_BENCH_ITER*.txt';
my @files = bsd_glob($pattern);
die "No files matched pattern: $pattern\n" unless @files;

my (%abs, %rel);
my (%metrics); # union of all metric keys seen in either section

for my $file (@files) {
    my ($iter) = ($file =~ /ITER(\d+)/i);
    next unless defined $iter;
    $iter = int($iter); # "001" -> 1

    open my $fh, '<', $file or die "Cannot open $file: $!\n";

    my $section = ''; # '', 'abs', 'rel'
    while (my $line = <$fh>) {
        chomp $line;

        if ($line =~ /^\*{3}\s*TIMINGS\s*\(s\)\s*\*{3}/) {
            $section = 'abs';
            next;
        }
        if ($line =~ /^\*{3}\s*RELATIVE\s*TIMINGS\s*\(%\)\s*\*{3}/) {
            $section = 'rel';
            next;
        }

        # Parse "key : value" lines inside the active section
        if ($section && $line =~ /^\s*(.+?)\s*:\s*([0-9]*\.?[0-9]+)\s*$/) {
            my ($key, $val) = ($1, $2);

            $key =~ s/^\s+//;
            $key =~ s/\s+$//;

            $metrics{$key} = 1;

            if ($section eq 'abs') {
                $abs{$iter}{$key} = $val + 0;
            } elsif ($section eq 'rel') {
                $rel{$iter}{$key} = $val + 0;
            }
        }
    }

    close $fh;
}

# Metrics sorted (swap this for a fixed order if you prefer)
my @cols = sort keys %metrics;

# Iterations sorted numerically (union of those seen in abs/rel)
my %all_iters = map { $_ => 1 } (keys %abs, keys %rel);
my @iters = sort { $a <=> $b } keys %all_iters;

my $out = 'timings_combined.csv';
open my $ofh, '>', $out or die "Cannot write $out: $!\n";

# Header: iteration, sec:<metric>..., pct:<metric>...
my @header = ('iteration',
              (map { "sec:$_" } @cols),
              (map { "pct:$_" } @cols));
print $ofh join(',', map { csv_escape($_) } @header), "\n";

for my $iter (@iters) {
    my @row = ($iter);

    for my $m (@cols) {
        push @row, (exists $abs{$iter}{$m} ? $abs{$iter}{$m} : '');
    }
    for my $m (@cols) {
        push @row, (exists $rel{$iter}{$m} ? $rel{$iter}{$m} : '');
    }

    print $ofh join(',', map { csv_escape($_) } @row), "\n";
}

close $ofh;
print "Wrote $out\n";

sub csv_escape {
    my ($s) = @_;
    $s = '' unless defined $s;
    if ($s =~ /[,"\n]/) {
        $s =~ s/"/""/g;
        return qq("$s");
    }
    return $s;
}
