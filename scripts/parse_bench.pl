#!/usr/bin/env perl
use strict;
use warnings;
use File::Glob qw(bsd_glob);
use File::Path qw(make_path);

# Parse benchmark text files (*_BENCH_ITER*.txt) and emit ONE "matrix" CSV:
#   - matrix_seconds.csv  (rows = iteration, cols = tags/metrics, values = seconds)
#
# Also: "probability table" seconds are added into "3D alignment" seconds,
# and "probability table" is NOT emitted as a separate column.
#
# Usage:
#   perl parse_bench.pl [glob] [outdir]
#
# Examples:
#   perl parse_bench.pl "CLUSTER2D_BENCH_ITER*.txt"
#   perl parse_bench.pl "*_BENCH_ITER*.txt" out_csv
#
# Output file:
#   outdir/matrix_seconds.csv

my $pattern = shift(@ARGV) // '*_BENCH_ITER*.txt';
my $outdir  = shift(@ARGV) // '.';

my @files = bsd_glob($pattern);
die "No files matched pattern: $pattern\n" unless @files;

make_path($outdir) unless -d $outdir;

# Data:
#   sec{iter}{metric} = seconds
my (%sec);
my (%metrics, %iters);

for my $file (@files) {
    my ($iter) = ($file =~ /ITER(\d+)/i);
    next unless defined $iter;
    $iter = int($iter);
    $iters{$iter} = 1;

    open my $fh, '<', $file or die "Cannot open $file: $!\n";

    my $in_seconds = 0;
    while (my $line = <$fh>) {
        chomp $line;

        if ($line =~ /^\s*\*{3}\s*TIMINGS\s*\(s\)\s*\*{3}\s*$/i) {
            $in_seconds = 1;
            next;
        }

        # Ignore percentage section entirely (no parsing, no data capture)
        if ($line =~ /^\s*\*{3}\s*RELATIVE\s*TIMINGS\s*\(%\)\s*\*{3}\s*$/i) {
            $in_seconds = 0;
            next;
        }

        # Parse only seconds metrics in TIMINGS (s) section:
        #   stochastic alignment : 0.50
        #   3D alignment         : 1.23
        if ($in_seconds && $line =~ /^\s*(.+?)\s*:\s*([0-9]+(?:\.[0-9]+)?)\s*$/) {
            my ($metric, $val) = ($1, $2);

            $metric =~ s/^\s+//;
            $metric =~ s/\s+$//;
            $metric =~ s/\s{2,}/ /g;  # normalize internal whitespace

            $metrics{$metric} = 1;
            $sec{$iter}{$metric} = $val + 0;
        }
    }

    close $fh;
}

# --- Merge "probability table" into "3D alignment" and drop it as a separate metric ---
my $m3d = find_metric_key_ci(\%metrics, '3D alignment');
my $mpt = find_metric_key_ci(\%metrics, 'probability table');

if (defined $mpt) {
    $m3d //= '3D alignment';
    $metrics{$m3d} = 1;

    for my $iter (keys %iters) {
        if (exists $sec{$iter}{$mpt}) {
            $sec{$iter}{$m3d} = ($sec{$iter}{$m3d} // 0) + $sec{$iter}{$mpt};
            delete $sec{$iter}{$mpt};
        }
    }

    delete $metrics{$mpt};
}

my @iter_list = sort { $a <=> $b } keys %iters;

# Sort metrics; push summary-ish rows to end if present
my @metric_list = sort {
    metric_rank($a) <=> metric_rank($b) || lc($a) cmp lc($b)
} keys %metrics;

write_matrix("$outdir/matrix_seconds.csv", \@iter_list, \@metric_list, \%sec);

print "Matched " . scalar(@files) . " files from pattern: $pattern\n";
print "Wrote:\n  $outdir/matrix_seconds.csv\n";

sub write_matrix {
    my ($out, $iters_ref, $metrics_ref, $data_ref) = @_;

    open my $ofh, '>', $out or die "Cannot write $out: $!\n";

    # Header: iteration,<metric1>,<metric2>,...
    print $ofh join(',', map { csv_escape($_) } ('iteration', @$metrics_ref)), "\n";

    for my $iter (@$iters_ref) {
        my @row = ($iter);

        for my $m (@$metrics_ref) {
            my $v = '';
            if (exists $data_ref->{$iter} && exists $data_ref->{$iter}{$m}) {
                $v = $data_ref->{$iter}{$m};
            }
            push @row, $v;
        }

        print $ofh join(',', map { csv_escape($_) } @row), "\n";
    }

    close $ofh;
}

sub csv_escape {
    my ($s) = @_;
    $s = '' unless defined $s;
    if ($s =~ /[,"\n]/) {
        $s =~ s/"/""/g;
        return qq("$s");
    }
    return $s;
}

sub metric_rank {
    my ($m) = @_;
    my $lc = lc($m);
    return 1000 if $lc eq 'total time';
    return 0;
}

sub find_metric_key_ci {
    my ($href, $target) = @_;
    my $t = lc($target // '');
    for my $k (keys %$href) {
        return $k if lc($k) eq $t;
    }
    return undef;
}