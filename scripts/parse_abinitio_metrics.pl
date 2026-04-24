#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $infile  = '';
my $outfile = '';
my $format  = 'tsv';   # tsv|csv

GetOptions(
  'in=s'     => \$infile,
  'out=s'    => \$outfile,
  'format=s' => \$format,
) or die usage();

die usage() if !$infile;
$format = lc($format);
die "ERROR: --format must be tsv or csv\n" if $format ne 'tsv' && $format ne 'csv';
my $sep = ($format eq 'csv') ? ',' : "\t";

open(my $IN, '<', $infile) or die "ERROR: cannot open --in '$infile': $!\n";

my $OUT;
if ($outfile) {
  open($OUT, '>', $outfile) or die "ERROR: cannot open --out '$outfile': $!\n";
} else {
  $OUT = *STDOUT;
}

# Current context
my ($stage, $lp) = (undef, undef);
my $cur = undef;          # hashref for current iteration record
my @rows;

# For wrapped "SHIFT INCR ARG" line
my $pending_shift = 0;

# Track detected states for per-state resolution columns
my %detected_states = ();  # $state => 1 if seen

sub new_record {
  my ($iter) = @_;
  my $r = {
    stage => defined($stage) ? $stage : 'NA',
    lp    => defined($lp)    ? $lp    : 'NA',
    iteration => $iter,

    overlap => 'NA',
    pct_sampled => 'NA',
    pct_updated => 'NA',
    pct_used_model => 'NA',
    pct_greedy => 'NA',

    dist_best_avg => 'NA', dist_best_sdev => 'NA', dist_best_min => 'NA', dist_best_max => 'NA',
    inplane_avg   => 'NA', inplane_sdev   => 'NA', inplane_min   => 'NA', inplane_max   => 'NA',
    shift_avg     => 'NA', shift_sdev     => 'NA', shift_min     => 'NA', shift_max     => 'NA',

    scanned_avg   => 'NA', scanned_sdev   => 'NA', scanned_min   => 'NA', scanned_max   => 'NA',
    pwt_avg       => 'NA', pwt_sdev       => 'NA', pwt_min       => 'NA', pwt_max       => 'NA',
    matchlp_avg   => 'NA', matchlp_sdev   => 'NA', matchlp_min   => 'NA', matchlp_max   => 'NA',
    estlp_avg     => 'NA', estlp_sdev     => 'NA', estlp_min     => 'NA', estlp_max     => 'NA',
    fsc143_avg    => 'NA', fsc143_sdev    => 'NA', fsc143_min    => 'NA', fsc143_max    => 'NA',
    score_avg     => 'NA', score_sdev     => 'NA', score_min     => 'NA', score_max     => 'NA',

    refinement_mode => 'NA',
    trailing_rec_update_fraction => 'NA',
    converged => 'NA',
  };
  # Initialize per-state resolution fields for all previously detected states
  foreach my $state (sort {$a <=> $b} keys %detected_states) {
    my $key = 'res_state_' . sprintf('%02d', $state);
    $r->{$key} = 'NA';
  }
  return $r;
}

sub finalize_record {
  my ($r) = @_;
  return if !$r;
  push @rows, $r;
}

# Parse a "AVG/SDEV/MIN/MAX" numeric quadruple at end of a line
sub parse_quad {
  my ($line) = @_;
  # captures last 4 numeric-ish tokens
  if ($line =~ /:\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$/) {
    return ($1,$2,$3,$4);
  }
  return;
}

while (my $line = <$IN>) {
  chomp $line;

  # Stage marker
  # Example: >>> STAGE   1 WITH LP   9.8
  if ($line =~ /^>>>\s*STAGE\s+(\d+)\s+WITH\s+LP\s+([0-9.]+)/) {
    $stage = $1;
    $lp    = $2;
    next;
  }

  # Iteration marker
  # Example: >>> ITERATION      12
  if ($line =~ /^>>>\s*ITERATION\s+(\d+)/) {
    finalize_record($cur);
    $cur = new_record($1);
    $pending_shift = 0;
    next;
  }

  next if !$cur; # ignore header noise before first iteration

  # Handle wrapped shift line:
  # >>> SHIFT INCR ARG           AVG/SDEV/MIN/
  # MAX:        0.000 ...
  if ($line =~ /^>>>\s*SHIFT\s+INCR\s+ARG\s+AVG\/SDEV\/MIN\/\s*$/) {
    $pending_shift = 1;
    next;
  }
  if ($pending_shift && $line =~ /^MAX:\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$/) {
    ($cur->{shift_avg}, $cur->{shift_sdev}, $cur->{shift_min}, $cur->{shift_max}) = ($1,$2,$3,$4);
    $pending_shift = 0;
    next;
  }
  # If something else happens, stop waiting
  $pending_shift = 0 if $pending_shift && $line !~ /^MAX:/;

  # Simple scalar metrics
  if ($line =~ /^>>>\s*ORIENTATION\s+OVERLAP:\s+([-\d.]+)/) {
    $cur->{overlap} = $1; next;
  }
  if ($line =~ /^>>>\s*%\s*PARTICLES\s+SAMPLED\s+THIS\s+ITERATION\s+([-\d.]+)/) {
    $cur->{pct_sampled} = $1; next;
  }
  if ($line =~ /^>>>\s*%\s*PARTICLES\s+UPDATED\s+SO\s+FAR\s+([-\d.]+)/) {
    $cur->{pct_updated} = $1; next;
  }
  if ($line =~ /^>>>\s*%\s*PARTICLES\s+USED\s+FOR\s+UPDATING\s+THE\s+MODEL\s+([-\d.]+)/) {
    $cur->{pct_used_model} = $1; next;
  }
  if ($line =~ /^>>>\s*%\s*GREEDY\s+SEARCHES\s+([-\d.]+)/) {
    $cur->{pct_greedy} = $1; next;
  }

  # Quad metrics
  if ($line =~ /^>>>\s*DIST\s+BTW\s+BEST\s+ORIS.*AVG\/SDEV\/MIN\/MAX:/) {
    my @q = parse_quad($line);
    if (@q) { @{$cur}{qw(dist_best_avg dist_best_sdev dist_best_min dist_best_max)} = @q; }
    next;
  }
  if ($line =~ /^>>>\s*IN-PLANE\s+DIST.*AVG\/SDEV\/MIN\/MAX:/) {
    my @q = parse_quad($line);
    if (@q) { @{$cur}{qw(inplane_avg inplane_sdev inplane_min inplane_max)} = @q; }
    next;
  }
  if ($line =~ /^>>>\s*%\s*SEARCH\s+SPACE\s+SCANNED.*AVG\/SDEV\/MIN\/MAX:/) {
    my @q = parse_quad($line);
    if (@q) { @{$cur}{qw(scanned_avg scanned_sdev scanned_min scanned_max)} = @q; }
    next;
  }
  if ($line =~ /^>>>\s*PARTICLE\s+WEIGHTS.*AVG\/SDEV\/MIN\/MAX:/) {
    my @q = parse_quad($line);
    if (@q) { @{$cur}{qw(pwt_avg pwt_sdev pwt_min pwt_max)} = @q; }
    next;
  }
  if ($line =~ /^>>>\s*MATCHING\s+LOW-PASS\s+LIMIT.*AVG\/SDEV\/MIN\/MAX:/) {
    my @q = parse_quad($line);
    if (@q) { @{$cur}{qw(matchlp_avg matchlp_sdev matchlp_min matchlp_max)} = @q; }
    next;
  }
  if ($line =~ /^>>>\s*ESTIMATED\s+LOW-PASS\s+LIMIT.*AVG\/SDEV\/MIN\/MAX:/) {
    my @q = parse_quad($line);
    if (@q) { @{$cur}{qw(estlp_avg estlp_sdev estlp_min estlp_max)} = @q; }
    next;
  }
  if ($line =~ /^>>>\s*RESOLUTION\s+\@\s*FSC=0\.143.*AVG\/SDEV\/MIN\/MAX:/) {
    my @q = parse_quad($line);
    if (@q) { @{$cur}{qw(fsc143_avg fsc143_sdev fsc143_min fsc143_max)} = @q; }
    next;
  }

  # Per-state resolution: >>> RESOLUTION @ FSC=0.143 STATE 01 AVG/SDEV/MIN/MAX: ...
  if ($line =~ /^>>>\s*RESOLUTION\s+\@\s*FSC=0\.143\s+STATE\s+(\d+)\s+AVG\/SDEV\/MIN\/MAX:/) {
    my $state = $1;
    my @q = parse_quad($line);
    if (@q) {
      # Register this state if new
      if (!exists $detected_states{$state}) {
        $detected_states{$state} = 1;
        # Backfill all existing rows with this new state key
        my $key = 'res_state_' . sprintf('%02d', $state);
        foreach my $r (@rows) {
          $r->{$key} = 'NA' if !exists $r->{$key};
        }
      }
      # Store the resolution value
      my $key = 'res_state_' . sprintf('%02d', $state);
      $cur->{$key} = $q[0];  # Use AVG value
    }
    next;
  }

  if ($line =~ /^>>>\s*SCORE\s+\[0,1\].*AVG\/SDEV\/MIN\/MAX:/) {
    my @q = parse_quad($line);
    if (@q) { @{$cur}{qw(score_avg score_sdev score_min score_max)} = @q; }
    next;
  }

  # Other flags / modes
  if ($line =~ /^>>>\s*REFINEMENT\s+MODE\s+IS\s+(\S+)/) {
    $cur->{refinement_mode} = $1; next;
  }
  if ($line =~ /^>>>\s*TRAILING\s+REC\s+UPDATE\s+FRACTION:\s*([-\d.]+)/) {
    $cur->{trailing_rec_update_fraction} = $1; next;
  }
  if ($line =~ /^>>>\s*CONVERGED:\s*(\.\w+\.)/) {
    $cur->{converged} = $1; next;
  }
}

finalize_record($cur);
close $IN;

# Build dynamic state-wise resolution header columns
my @state_cols = map { 'res_state_' . sprintf('%02d', $_) } sort {$a <=> $b} keys %detected_states;

my @header = qw(
  stage lp iteration
  overlap pct_sampled pct_updated pct_used_model pct_greedy
  dist_best_avg dist_best_sdev dist_best_min dist_best_max
  inplane_avg inplane_sdev inplane_min inplane_max
  shift_avg shift_sdev shift_min shift_max
  scanned_avg scanned_sdev scanned_min scanned_max
  pwt_avg pwt_sdev pwt_min pwt_max
  matchlp_avg matchlp_sdev matchlp_min matchlp_max
  estlp_avg estlp_sdev estlp_min estlp_max
  fsc143_avg fsc143_sdev fsc143_min fsc143_max
  score_avg score_sdev score_min score_max
  refinement_mode trailing_rec_update_fraction converged
);

# Append per-state resolution columns to header
push @header, @state_cols;

print {$OUT} join($sep, @header) . "\n";
for my $r (@rows) {
  my @vals = map { exists $r->{$_} ? $r->{$_} : 'NA' } @header;
  print {$OUT} join($sep, @vals) . "\n";
}
close $OUT if $outfile;

exit 0;

sub usage {
  return <<"USAGE";
Usage:
  perl parse_abinitio_metrics.pl --in ABINITIO3D_OUTPUT_RESTART1 [--out metrics.tsv] [--format tsv|csv]

Examples:
  perl parse_abinitio_metrics.pl --in ABINITIO3D_OUTPUT_RESTART1 > metrics.tsv
  perl parse_abinitio_metrics.pl --in ABINITIO3D_OUTPUT_RESTART1 --format csv --out metrics.csv
USAGE
}
