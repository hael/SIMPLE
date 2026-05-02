#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(basename);
use File::Glob qw(bsd_glob);
use File::Path qw(make_path);

# Parse benchmark text files (*_BENCH_ITER*.txt) and emit matrix CSVs.
#
# New sectioned benchmark files may contain:
#   *** BENCHMARK CONTEXT ***
#   *** COMPARABLE TIMINGS (s) ***
#   *** COMPARABLE DETAIL TIMINGS (s) ***
#   *** COMPARABLE RELATIVE TIMINGS (%) ***
#
# The parser prefers the COMPARABLE sections when present and falls back to the
# older TIMINGS / RELATIVE TIMINGS sections for files that have not been updated.
#
# Usage:
#   perl parse_bench.pl [glob] [outdir]
#
# Examples:
#   perl parse_bench.pl "REFINE3D*_BENCH_ITER*.txt"
#   perl parse_bench.pl "*_BENCH_ITER*.txt" out_csv
#
# Output files:
#   outdir/matrix_seconds.csv
#   outdir/matrix_percent.csv
#   outdir/matrix_context.csv

my $pattern = shift(@ARGV) // '*_BENCH_ITER*.txt';
my $outdir  = shift(@ARGV) // '.';

my @files = bsd_glob($pattern);
die "No files matched pattern: $pattern\n" unless @files;

make_path($outdir) unless -d $outdir;

my (%sec, %pct);
my (%sec_metrics, %pct_metrics);
my @context_rows;
my %context_keys;
my %iters;

for my $file (@files) {
    my ($iter) = ($file =~ /ITER(\d+)/i);
    next unless defined $iter;
    $iter = int($iter);
    $iters{$iter} = 1;

    my ($context, $file_sec, $file_pct) = parse_bench_file($file);

    for my $metric (keys %{$file_sec}) {
        $sec{$iter}{$metric} = $file_sec->{$metric};
        $sec_metrics{$metric} = 1;
    }
    for my $metric (keys %{$file_pct}) {
        $pct{$iter}{$metric} = $file_pct->{$metric};
        $pct_metrics{$metric} = 1;
    }

    my %context_row = (
        file      => basename($file),
        benchmark => benchmark_name($file),
        iteration => $iter,
        %{$context},
    );
    push @context_rows, \%context_row;
    $context_keys{$_} = 1 for keys %context_row;
}

my @iter_list = sort { $a <=> $b } keys %iters;
my @sec_metric_list = sort {
    metric_rank($a) <=> metric_rank($b) || lc($a) cmp lc($b)
} keys %sec_metrics;
my @pct_metric_list = sort {
    metric_rank($a) <=> metric_rank($b) || lc($a) cmp lc($b)
} keys %pct_metrics;

write_matrix("$outdir/matrix_seconds.csv", \@iter_list, \@sec_metric_list, \%sec);
write_matrix("$outdir/matrix_percent.csv", \@iter_list, \@pct_metric_list, \%pct);
write_context("$outdir/matrix_context.csv", \@context_rows, \%context_keys);

print "Matched " . scalar(@files) . " files from pattern: $pattern\n";
print "Wrote:\n";
print "  $outdir/matrix_seconds.csv\n";
print "  $outdir/matrix_percent.csv\n";
print "  $outdir/matrix_context.csv\n";

sub parse_bench_file {
    my ($file) = @_;

    open my $fh, '<', $file or die "Cannot open $file: $!\n";

    my %context;
    my (%sec_comparable, %sec_detail, %sec_legacy);
    my (%pct_comparable, %pct_legacy);
    my $section = '';

    while (my $line = <$fh>) {
        chomp $line;

        if ($line =~ /^\s*\*{3}\s*(.*?)\s*\*{3}\s*$/) {
            my $header = normalize_header($1);
            if ($header eq 'BENCHMARK CONTEXT') {
                $section = 'context';
            } elsif ($header eq 'COMPARABLE TIMINGS (S)') {
                $section = 'sec_comparable';
            } elsif ($header eq 'COMPARABLE DETAIL TIMINGS (S)') {
                $section = 'sec_detail';
            } elsif ($header eq 'TIMINGS (S)') {
                $section = 'sec_legacy';
            } elsif ($header eq 'COMPARABLE RELATIVE TIMINGS (%)') {
                $section = 'pct_comparable';
            } elsif ($header eq 'RELATIVE TIMINGS (%)') {
                $section = 'pct_legacy';
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

        if ($section eq 'sec_comparable') {
            $sec_comparable{$metric} = $value;
        } elsif ($section eq 'sec_detail') {
            $sec_detail{$metric} = $value;
        } elsif ($section eq 'sec_legacy') {
            $sec_legacy{$metric} = $value;
        } elsif ($section eq 'pct_comparable') {
            $pct_comparable{$metric} = $value;
        } elsif ($section eq 'pct_legacy') {
            $pct_legacy{$metric} = $value;
        }
    }

    close $fh;

    my %file_sec = ((keys %sec_comparable) || (keys %sec_detail))
        ? (%sec_comparable, %sec_detail)
        : %sec_legacy;
    my %file_pct = (keys %pct_comparable) ? %pct_comparable : %pct_legacy;

    return (\%context, \%file_sec, \%file_pct);
}

sub parse_metric_value {
    my ($line) = @_;
    my $num = qr/[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?/;
    return unless $line =~ /^\s*(.+?)\s*:\s*($num)\s*$/;

    my ($metric, $value) = (normalize_key($1), $2 + 0);
    return ($metric, $value);
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
    my ($out, $iters_ref, $metrics_ref, $data_ref) = @_;

    open my $ofh, '>', $out or die "Cannot write $out: $!\n";
    print $ofh join(',', map { csv_escape($_) } ('iteration', @{$metrics_ref})), "\n";

    for my $iter (@{$iters_ref}) {
        my @row = ($iter);
        for my $metric (@{$metrics_ref}) {
            my $value = '';
            if (exists $data_ref->{$iter} && exists $data_ref->{$iter}{$metric}) {
                $value = $data_ref->{$iter}{$metric};
            }
            push @row, $value;
        }
        print $ofh join(',', map { csv_escape($_) } @row), "\n";
    }

    close $ofh;
}

sub write_context {
    my ($out, $rows_ref, $keys_ref) = @_;

    my @prefix = qw(file benchmark iteration);
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
