#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(basename);
use File::Glob qw(bsd_glob);
use File::Path qw(make_path);
use File::Temp qw(tempdir);

# Parse benchmark text files (*_BENCH_ITER*.txt) and emit matrix CSVs.
#
# Benchmark files may contain:
#   *** BENCHMARK CONTEXT ***
#   *** TIMINGS (s) ***
#   *** RELATIVE TIMINGS (%) ***
#
# Any header ending in TIMINGS (s) is treated as a seconds section. Percentage
# metrics in that section, such as "% accounted for", are kept in the percent
# matrix. Missing percentages are derived from the parsed total-time entry.
#
# A matrix row is one benchmark family at one iteration (columns "benchmark" and
# "iteration"), so families matched by one glob never overwrite each other.
# Per-partition reports (*_PARTppp.txt, e.g. REFINE3D_BENCH_ITERnnn_PARTppp.txt)
# are represented by partition 1, as in plot_refine3d_bench.py; the other
# partitions are skipped and counted in the summary.
#
# Usage:
#   perl parse_bench.pl [glob] [outdir]
#   perl parse_bench.pl --self-test
#
# Examples:
#   perl parse_bench.pl "REFINE3D*_BENCH_ITER*.txt"
#   perl parse_bench.pl "*_BENCH_ITER*.txt" out_csv
#
# Output files:
#   outdir/matrix_seconds.csv
#   outdir/matrix_percent.csv
#   outdir/matrix_context.csv
#
# --self-test writes fixture files for two partitions of one iteration plus a
# volassemble report into a temporary directory, parses them and checks that
# the matrices hold partition 1 and keep the families apart. Exit status 0 on
# success, 1 on failure.

if (@ARGV && $ARGV[0] eq '--self-test') {
    exit(self_test());
}

my $pattern = shift(@ARGV) // '*_BENCH_ITER*.txt';
my $outdir  = shift(@ARGV) // '.';
my ($nused, $nskipped) = run_parse($pattern, $outdir);
print "Matched " . ($nused + $nskipped) . " files from pattern: $pattern";
print " ($nskipped partition reports other than partition 1 skipped)" if $nskipped;
print "\n";
print "Wrote:\n";
print "  $outdir/matrix_seconds.csv\n";
print "  $outdir/matrix_percent.csv\n";
print "  $outdir/matrix_context.csv\n";
exit 0;

sub run_parse {
    my ($pattern, $outdir) = @_;
    my @files = bsd_glob($pattern);
    die "No files matched pattern: $pattern\n" unless @files;

    make_path($outdir) unless -d $outdir;

    my (%sec, %pct);
    my (%sec_metrics, %pct_metrics);
    my @context_rows;
    my %context_keys;
    my %rows;
    my ($nused, $nskipped) = (0, 0);

    for my $file (@files) {
        my ($iter) = ($file =~ /ITER(\d+)/i);
        next unless defined $iter;
        $iter = int($iter);
        my ($part) = (basename($file) =~ /_PART(\d+)/i);
        if (defined $part && int($part) != 1) {
            $nskipped++;
            next;
        }
        $nused++;
        my $bench = benchmark_name($file);
        my $row   = "$bench\t$iter";
        $rows{$row} = [$bench, $iter];

        my ($context, $file_sec, $file_pct) = parse_bench_file($file);

        for my $metric (keys %{$file_sec}) {
            $sec{$row}{$metric} = $file_sec->{$metric};
            $sec_metrics{$metric} = 1;
        }
        for my $metric (keys %{$file_pct}) {
            $pct{$row}{$metric} = $file_pct->{$metric};
            $pct_metrics{$metric} = 1;
        }

        my %context_row = (
            file      => basename($file),
            benchmark => $bench,
            iteration => $iter,
            partition => (defined $part ? int($part) : ''),
            %{$context},
        );
        push @context_rows, \%context_row;
        $context_keys{$_} = 1 for keys %context_row;
    }

    my @row_list = sort {
        $rows{$a}[0] cmp $rows{$b}[0] || $rows{$a}[1] <=> $rows{$b}[1]
    } keys %rows;
    my @sec_metric_list = sort {
        metric_rank($a) <=> metric_rank($b) || lc($a) cmp lc($b)
    } keys %sec_metrics;
    my @pct_metric_list = sort {
        metric_rank($a) <=> metric_rank($b) || lc($a) cmp lc($b)
    } keys %pct_metrics;

    write_matrix("$outdir/matrix_seconds.csv", \@row_list, \%rows, \@sec_metric_list, \%sec);
    write_matrix("$outdir/matrix_percent.csv", \@row_list, \%rows, \@pct_metric_list, \%pct);
    write_context("$outdir/matrix_context.csv", \@context_rows, \%context_keys);
    return ($nused, $nskipped);
}

sub self_test {
    my $dir = tempdir(CLEANUP => 1);
    my %fixtures = (
        'REFINE3D_BENCH_ITER001_PART001.txt' => 11.0,
        'REFINE3D_BENCH_ITER001_PART002.txt' => 99.0,
        'VOLASSEMBLE_BENCH_ITER001.txt'      => 5.0,
    );
    for my $name (keys %fixtures) {
        open my $fh, '>', "$dir/$name" or die "Cannot write fixture $name: $!\n";
        print $fh "*** TIMINGS (s) ***\n";
        print $fh "matching : $fixtures{$name}\n";
        print $fh "total time : " . ($fixtures{$name} * 2) . "\n";
        close $fh;
    }
    my ($nused, $nskipped) = run_parse("$dir/*_BENCH_ITER*.txt", $dir);
    my @fail;
    push @fail, "expected 2 files used and 1 skipped, got $nused and $nskipped"
        unless $nused == 2 && $nskipped == 1;
    open my $fh, '<', "$dir/matrix_seconds.csv" or die "Cannot read matrix: $!\n";
    my @lines = <$fh>;
    close $fh;
    chomp @lines;
    my @header = split /,/, shift @lines;
    my ($imatch) = grep { $header[$_] eq 'matching' } 0..$#header;
    push @fail, 'no matching column' unless defined $imatch;
    my %got;
    for my $line (@lines) {
        my @f = split /,/, $line;
        $got{"$f[0] $f[1]"} = $f[$imatch] if defined $imatch;
    }
    push @fail, 'REFINE3D_BENCH iteration 1 must hold partition 1 (11)'
        unless defined $got{'REFINE3D_BENCH 1'} && $got{'REFINE3D_BENCH 1'} == 11;
    push @fail, 'VOLASSEMBLE_BENCH iteration 1 must be its own row (5)'
        unless defined $got{'VOLASSEMBLE_BENCH 1'} && $got{'VOLASSEMBLE_BENCH 1'} == 5;
    push @fail, 'expected exactly two matrix rows' unless scalar(@lines) == 2;
    if (@fail) {
        print "parse_bench self-test FAILED:\n";
        print "  $_\n" for @fail;
        return 1;
    }
    print "parse_bench self-test passed\n";
    return 0;
}

sub parse_bench_file {
    my ($file) = @_;

    open my $fh, '<', $file or die "Cannot open $file: $!\n";

    my %context;
    my (%sec, %pct);
    my $section = '';

    while (my $line = <$fh>) {
        chomp $line;

        if ($line =~ /^\s*\*{3}\s*(.*?)\s*\*{3}\s*$/) {
            my $header = normalize_header($1);
            if ($header eq 'BENCHMARK CONTEXT') {
                $section = 'context';
            } elsif ($header =~ /TIMINGS \(S\)$/) {
                $section = 'sec';
            } elsif ($header =~ /RELATIVE TIMINGS \(%\)$/) {
                $section = 'pct';
            } else {
                $section = '';
            }
            next;
        }

        if ($section eq 'context') {
            if ($line =~ /^\s*(.+?)\s*:\s*(.*?)\s*$/) {
                my ($key, $value) = (normalize_key($1), $2);
                $value =~ s/^\s+//;
                $value =~ s/\s+$//;
                $context{$key} = $value;
            }
            next;
        }

        next unless $section;

        my ($metric, $value) = parse_metric_value($line);
        next unless defined $metric;

        if (is_percent_metric($metric)) {
            $pct{$metric} = $value;
        } elsif ($section eq 'sec') {
            $sec{$metric} = $value;
        } elsif ($section eq 'pct') {
            $pct{$metric} = $value;
        }
    }

    close $fh;

    my %derived_pct = derive_percentages(\%sec);
    for my $metric (keys %derived_pct) {
        $pct{$metric} = $derived_pct{$metric} unless exists $pct{$metric};
    }

    return (\%context, \%sec, \%pct);
}

sub derive_percentages {
    my ($sec_ref) = @_;
    my %pct;
    my $total_metric = find_total_metric($sec_ref);
    return %pct unless defined $total_metric;

    my $total = $sec_ref->{$total_metric};
    return %pct unless defined $total && $total > 0;

    for my $metric (keys %{$sec_ref}) {
        next if $metric eq $total_metric;
        next if is_percent_metric($metric);
        $pct{$metric} = ($sec_ref->{$metric} / $total) * 100.0;
    }

    my $accounted = 0.0;
    for my $value (values %pct) {
        $accounted += $value;
    }
    my $prefix = $total_metric;
    $prefix =~ s/\s*total time\s*$//i;
    $prefix =~ s/\s+$//;
    my $accounted_metric = length($prefix) ? "$prefix % accounted for" : '% accounted for';
    $pct{$accounted_metric} = $accounted;

    return %pct;
}

sub find_total_metric {
    my ($sec_ref) = @_;
    for my $metric (sort keys %{$sec_ref}) {
        return $metric if lc($metric) =~ /\btotal time$/;
    }
    return;
}

sub parse_metric_value {
    my ($line) = @_;
    my $num = qr/[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?/;
    return unless $line =~ /^\s*(.+?)\s*:\s*($num)\s*$/;

    my ($metric, $value) = (normalize_key($1), $2 + 0);
    return ($metric, $value);
}

sub is_percent_metric {
    my ($metric) = @_;
    return lc($metric) =~ /%/;
}

sub normalize_header {
    my ($header) = @_;
    $header =~ s/^\s+//;
    $header =~ s/\s+$//;
    $header =~ s/\s+/ /g;
    return uc($header);
}

sub normalize_key {
    my ($key) = @_;
    $key =~ s/^\s+//;
    $key =~ s/\s+$//;
    $key =~ s/\s{2,}/ /g;
    return $key;
}

sub benchmark_name {
    my ($file) = @_;
    my $name = basename($file);
    $name =~ s/_?ITER\d+.*$//i;
    return $name;
}

sub write_matrix {
    my ($out, $rows_ref, $ids_ref, $metrics_ref, $data_ref) = @_;

    open my $ofh, '>', $out or die "Cannot write $out: $!\n";
    print $ofh join(',', map { csv_escape($_) } ('benchmark', 'iteration', @{$metrics_ref})), "\n";

    for my $row (@{$rows_ref}) {
        my @out = @{$ids_ref->{$row}};
        for my $metric (@{$metrics_ref}) {
            my $value = '';
            if (exists $data_ref->{$row} && exists $data_ref->{$row}{$metric}) {
                $value = $data_ref->{$row}{$metric};
            }
            push @out, $value;
        }
        print $ofh join(',', map { csv_escape($_) } @out), "\n";
    }

    close $ofh;
}

sub write_context {
    my ($out, $rows_ref, $keys_ref) = @_;

    my @prefix = qw(file benchmark iteration partition);
    my %is_prefix = map { $_ => 1 } @prefix;
    my @keys = (
        @prefix,
        sort { lc($a) cmp lc($b) } grep { !$is_prefix{$_} } keys %{$keys_ref},
    );

    open my $ofh, '>', $out or die "Cannot write $out: $!\n";
    print $ofh join(',', map { csv_escape($_) } @keys), "\n";

    for my $row (@{$rows_ref}) {
        print $ofh join(',', map { csv_escape($row->{$_} // '') } @keys), "\n";
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
    my ($metric) = @_;
    my $lc = lc($metric);
    return 1000 if $lc =~ /\btotal time$/;
    return 1001 if $lc =~ /% accounted for$/;
    return 0;
}
