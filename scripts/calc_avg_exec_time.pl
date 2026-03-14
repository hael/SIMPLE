#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(floor);

# Usage:
#   perl abinitio_stats.pl ABINITIO3D_OUTPUT_RESTART*
#   perl abinitio_stats.pl 'ABINITIO3D_OUTPUT_RESTART*'
#   perl abinitio_stats.pl file1 file2 ...

my @files = @ARGV;
die "Usage: $0 <files...>\\n" if !@files;

my (@times_sec, @max_iters, @missing_time, @missing_iter, @missing_files);

for my $f (@files) {
    if (!-e $f) { push @missing_files, $f; next; }
    open my $fh, '<', $f or do { warn "Could not open $f: $!\\n"; push @missing_files, $f; next; };

    my $last_time;
    my $max_iter;

    while (my $line = <$fh>) {
        # Execution time: allow with/without trailing '.', keep last match
        if ($line =~ /Execution\s+time:\s*([0-9]+(?:\.[0-9]+)?)\s*seconds\b/i) {
            $last_time = $1 + 0.0;
        }

        # Iteration: track maximum ITERATION value seen
        # MODIFIED: The regex is updated to match the format in the provided file.
        if ($line =~ /^>>>\s+ITERATION\s+(\d+)/) {
            my $it = $1 + 0;
            $max_iter = $it if !defined($max_iter) || $it > $max_iter;
        }
    }
    close $fh;

    if (defined $last_time) { push @times_sec, $last_time; }
    else { warn "No execution time found in $f\\n"; push @missing_time, $f; }

    if (defined $max_iter) { push @max_iters, $max_iter; }
    else { warn "No REFINE3D ITERATION found in $f\\n"; push @missing_iter, $f; }
}

# ---- Print stats ----
print_stats_time(\@times_sec) if @times_sec;
print "\\n";
print_stats_iters(\@max_iters) if @max_iters;

# ---- Report missing ----
my @any_missing = (@missing_files, @missing_time, @missing_iter);
if (@any_missing) {
    print "\\nIssues encountered:\\n";
    if (@missing_files) {
        print "  Missing/unreadable files:\\n";
        print "    $_\n" for @missing_files;
    }
    if (@missing_time) {
        print "  Files with no execution time found:\\n";
        print "    $_\n" for @missing_time;
    }
    if (@missing_iter) {
        print "  Files with no iteration lines found:\\n";
        print "    $_\n" for @missing_iter;
    }
}

exit 0;

# -------------------- helpers --------------------

sub mean_and_pop_stdev {
    my ($arr) = @_;
    my $n = scalar @$arr;
    die "No values.\\n" if $n == 0;
    my $sum = 0.0;
    $sum += $_ for @$arr;
    my $mean = $sum / $n;
    my $ss = 0.0;
    for my $x (@$arr) {
        my $d = $x - $mean;
        $ss += $d * $d;
    }
    my $stdev = sqrt($ss / $n); # population
    return ($n, $mean, $stdev);
}

sub fmt_hms {
    my ($sec) = @_;
    $sec = 0 if !defined $sec;
    $sec = 0 if $sec < 0;
    my $h = floor($sec / 3600);
    my $m = floor(($sec - 3600*$h) / 60);
    my $s = $sec - 3600*$h - 60*$m;
    my $is_int = (abs($s - int($s)) < 1e-9) ? 1 : 0;
    return $is_int
        ? sprintf("%d:%02d:%02d", $h, $m, int($s))
        : sprintf("%d:%02d:%05.1f", $h, $m, $s);
}

sub fmt_split_time {
    my ($sec) = @_;
    $sec = 0 if !defined $sec;
    $sec = 0 if $sec < 0;
    my $h = floor($sec / 3600);
    my $m = floor(($sec - 3600*$h) / 60);
    my $s = $sec - 3600*$h - 60*$m;
    my $is_int = (abs($s - int($s)) < 1e-9) ? 1 : 0;
    return $is_int
        ? sprintf("hours=%d minutes=%d seconds=%d", $h, $m, int($s))
        : sprintf("hours=%d minutes=%d seconds=%.1f", $h, $m, $s);
}

sub print_stats_time {
    my ($arr) = @_;
    my ($n, $mean, $stdev) = mean_and_pop_stdev($arr);
    print "Execution time stats (N=$n):\\n";
    print "  Mean:\\n";
    print "    seconds: $mean\\n";
    print "    h:m:s : ", fmt_hms($mean), "\\n";
    print "    split : ", fmt_split_time($mean), "\\n";
    print "  Population stdev:\\n";
    print "    seconds: $stdev\\n";
    print "    h:m:s : ", fmt_hms($stdev), "\\n";
    print "    split : ", fmt_split_time($stdev), "\\n";
}

sub print_stats_iters {
    my ($arr) = @_;
    my ($n, $mean, $stdev) = mean_and_pop_stdev($arr);
    print "Convergence iteration stats (max ITERATION per file, N=$n):\\n";
    print "  Mean max-iteration: $mean\\n";
    print "  Population stdev  : $stdev\\n";
}
