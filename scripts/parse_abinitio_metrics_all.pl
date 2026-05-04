#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(basename);
use File::Find qw(find);
use File::Path qw(make_path);
use Cwd qw(abs_path);

our $NUM_RE = qr/[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?/;

my $indir            = '.';
my $outdir           = 'abinitio_summary';
my $format           = 'tsv';       # tsv|csv
my $pattern          = '^ABINITIO3D_OUTPUT';
my $recursive        = 0;
my $ang_threshold    = 2.5;         # ANGTHRES_MI_PROJ_DEFAULT
my $overlap_limits_s = '0.80,0.85,0.90,0.92,0.95,0.98,0.99';
my $frac_limits_s    = 'stage';     # stage or comma-separated percentages
my $res_bin_width    = 0.25;        # Angstrom bins for resolution-vs-distance averages
my $res_benefit_min  = 0.05;        # Angstrom improvement counted as useful after first halt
my $detailed         = 0;
my $help             = 0;

GetOptions(
  'indir=s'          => \$indir,
  'outdir=s'         => \$outdir,
  'format=s'         => \$format,
  'pattern=s'        => \$pattern,
  'recursive!'       => \$recursive,
  'ang-threshold=f'  => \$ang_threshold,
  'overlap-limits=s' => \$overlap_limits_s,
  'fracsrch-limits=s'=> \$frac_limits_s,
  'res-bin-width=f'  => \$res_bin_width,
  'res-benefit-min=f'=> \$res_benefit_min,
  'detailed!'        => \$detailed,
  'help|h'           => \$help,
) or die usage();

if ($help) {
  print usage();
  exit 0;
}

$format = lc($format);
die "ERROR: --format must be tsv or csv\n" if $format ne 'tsv' && $format ne 'csv';
die "ERROR: --res-bin-width must be positive\n" if $res_bin_width <= 0.0;
die "ERROR: --res-benefit-min must be non-negative\n" if $res_benefit_min < 0.0;

$indir  = expand_tilde($indir);
$outdir = expand_tilde($outdir);
die "ERROR: --indir '$indir' is not a directory\n" if !-d $indir;
make_path($outdir) if !-d $outdir;

my @overlap_limits = parse_list($overlap_limits_s);
die "ERROR: --overlap-limits did not contain any numeric values\n" if !@overlap_limits;

my @frac_limits = ();
my $use_stage_frac = lc($frac_limits_s) eq 'stage';
@frac_limits = parse_list($frac_limits_s) if !$use_stage_frac;
die "ERROR: --fracsrch-limits did not contain any numeric values\n"
  if !$use_stage_frac && !@frac_limits;

my @files = find_input_files($indir, $pattern, $recursive);
die "ERROR: no files matching pattern /$pattern/ under '$indir'\n" if !@files;

my @iterations;
my %state_cols;
my @run_meta;

for my $file (@files) {
  my ($rows, $meta) = parse_abinitio_file($file);
  push @iterations, @$rows;
  push @run_meta, $meta;
  for my $r (@$rows) {
    for my $k (keys %$r) {
      $state_cols{$k} = 1 if $k =~ /^(?:res_state_|state\d+_(?:joint_overlap|population))/;
    }
  }
}

annotate_iterations(\@iterations, $ang_threshold);

my @halting_benefit     = build_halting_benefit(\@iterations, $res_benefit_min);
my @stage_averages      = build_stage_averages_compact(\@iterations);
my @resolution_distance = build_resolution_distance(\@iterations, $res_bin_width);
my @overall_summary     = build_overall_summary(\@iterations, \@run_meta, \@halting_benefit);

my @overall_header = qw(metric value note);
my @stage_avg_header = qw(
  stage lp_mean n_runs n_iters_mean n_iters_min n_iters_max
  res_start_mean res_end_mean res_best_mean res_gain_start_to_best_mean
  res_end_minus_best_mean
  dist_start_mean dist_end_mean dist_min_mean dist_drop_start_to_end_mean
  dist_end_minus_min_mean
  overlap_end_mean fracsrch_end_mean score_end_mean
);
my @resolution_distance_header = qw(
  resolution_bin_low resolution_bin_high n_iterations n_runs n_run_stages
  resolution_mean dist_best_oris_mean dist_best_oris_sdev
  overlap_mean fracsrch_mean score_mean stage_mean
);
my @halting_header = qw(
  stage lp_mean n_halt_points n_with_later_iters mean_later_iters_observed
  mean_res_at_first_halt mean_best_later_res mean_res_gain_later pct_later_res_improved
  mean_dist_at_first_halt mean_final_dist mean_dist_drop_after_halt
  mean_overlap_at_first_halt mean_final_overlap
);

emit_table("$outdir/overall_summary.$format",     \@overall_header,             \@overall_summary,     $format);
emit_table("$outdir/stage_averages.$format",      \@stage_avg_header,           \@stage_averages,      $format);
emit_table("$outdir/resolution_distance.$format", \@resolution_distance_header, \@resolution_distance, $format);
emit_table("$outdir/halting_benefit.$format",     \@halting_header,             \@halting_benefit,     $format);

if ($detailed) {
  my @stage_summary = build_stage_summary(\@iterations);
  my @run_summary   = build_run_summary(\@iterations, \@run_meta);
  my @scan          = build_convergence_scan(\@iterations, \@overlap_limits, \@frac_limits, $use_stage_frac);
  my @near_misses   = build_near_misses(\@iterations);

  my @dynamic_state_cols = sort state_col_sort keys %state_cols;

  my @iter_header = (
    qw(run_id restart_id file execution_dir normal_stop stage lp iteration stage_iter
       code_overlap_lim code_fracsrch_lim code_would_converge converged
       overlap overlap_gap fracsrch_avg fracsrch_gap
       pct_sampled n_sampled n_active pct_updated pct_used_model pct_greedy
       dist_best_avg dist_best_sdev dist_best_min dist_best_max
       dist_ang_threshold dist_best_avg_minus_ang_threshold dist_best_avg_over_ang_threshold
       nonoverlap_fraction nonoverlap_mean_lb nonoverlap_mean_ub nonoverlap_mean_lb_minus_ang_threshold
       inplane_avg inplane_sdev inplane_min inplane_max
       shift_avg shift_sdev shift_min shift_max
       scanned_avg scanned_sdev scanned_min scanned_max
       pwt_avg pwt_sdev pwt_min pwt_max
       matchlp_avg matchlp_sdev matchlp_min matchlp_max
       estlp_avg estlp_sdev estlp_min estlp_max
       fsc143_avg fsc143_sdev fsc143_min fsc143_max
       score_avg score_sdev score_min score_max
       refinement_mode nspace gau_regularization ml_regularization ml_tau
       icm_regularization icm_lambda automasking filter_mode
       trailing_rec_update_fraction fillin_particle_sampling)
  );
  push @iter_header, @dynamic_state_cols;

  my @stage_header = qw(
    run_id restart_id file stage lp n_iters first_iteration last_iteration
    first_logged_converged_iteration first_code_converged_iteration
    code_overlap_lim code_fracsrch_lim
    max_overlap last_overlap max_fracsrch last_fracsrch
    min_dist_best_avg last_dist_best_avg min_fsc143_avg last_fsc143_avg
    max_score_avg last_score_avg
    overlap_gap_last fracsrch_gap_last
  );

  my @run_header = qw(
    run_id restart_id file execution_dir normal_stop n_iters stages_seen
    first_stage last_stage first_iteration last_iteration
    first_logged_converged_iteration first_code_converged_iteration
    final_stage final_overlap final_fracsrch final_dist_best_avg final_fsc143_avg final_score_avg
    min_fsc143_avg min_fsc143_stage min_fsc143_iteration
    max_score_avg max_score_stage max_score_iteration
  );

  my @scan_header = qw(
    run_id restart_id file stage lp overlap_limit fracsrch_limit
    first_converged_iteration first_converged_stage_iter
    code_first_converged_iteration iters_saved_vs_code
    n_iters max_overlap max_fracsrch
  );

  my @miss_header = qw(
    run_id restart_id file stage lp iteration stage_iter
    overlap overlap_limit needed_extra_overlap
    fracsrch_avg fracsrch_limit fracsrch_margin
    dist_best_avg dist_best_sdev dist_best_max
    dist_ang_threshold dist_best_avg_minus_ang_threshold dist_best_avg_over_ang_threshold
    nonoverlap_fraction nonoverlap_mean_lb nonoverlap_mean_ub nonoverlap_mean_lb_minus_ang_threshold
    refinement_mode nspace
  );

  emit_table("$outdir/iterations.$format",       \@iter_header,  \@iterations,   $format);
  emit_table("$outdir/stage_summary.$format",    \@stage_header, \@stage_summary,$format);
  emit_table("$outdir/run_summary.$format",      \@run_header,   \@run_summary,  $format);
  emit_table("$outdir/convergence_scan.$format", \@scan_header,  \@scan,         $format);
  emit_table("$outdir/near_misses.$format",      \@miss_header,  \@near_misses,  $format);
}

print "Parsed " . scalar(@files) . " files and " . scalar(@iterations) . " iteration records\n";
print "Wrote:\n";
print "  $outdir/overall_summary.$format\n";
print "  $outdir/stage_averages.$format\n";
print "  $outdir/resolution_distance.$format\n";
print "  $outdir/halting_benefit.$format\n";
if ($detailed) {
  print "  $outdir/iterations.$format\n";
  print "  $outdir/stage_summary.$format\n";
  print "  $outdir/run_summary.$format\n";
  print "  $outdir/convergence_scan.$format\n";
  print "  $outdir/near_misses.$format\n";
}

exit 0;

sub parse_abinitio_file {
  my ($file) = @_;
  open(my $IN, '<', $file) or die "ERROR: cannot open '$file': $!\n";

  my $base = basename($file);
  my ($restart_id) = ($base =~ /RESTART(\d+)/);
  $restart_id = 'NA' if !defined $restart_id;

  my $meta = {
    run_id        => $base,
    restart_id    => $restart_id,
    file          => $file,
    execution_dir => 'NA',
    normal_stop   => 0,
  };

  my @rows;
  my ($stage, $lp) = ('NA', 'NA');
  my $stage_iter = 0;
  my $cur;
  my $pending_shift = 0;
  my $line_no = 0;

  while (my $line = <$IN>) {
    $line_no++;
    chomp $line;

    if ($line =~ /^>>>\s*EXECUTION\s+DIRECTORY:\s*(.+?)\s*$/) {
      $meta->{execution_dir} = $1;
      $cur->{execution_dir} = $1 if $cur;
      next;
    }
    if ($line =~ /SIMPLE_ABINITIO3D\s+NORMAL\s+STOP/) {
      $meta->{normal_stop} = 1;
      next;
    }

    if ($line =~ /^>>>\s*STAGE\s+(\d+)\s+WITH\s+LP\s+($NUM_RE)/) {
      finalize_record(\@rows, $cur);
      $cur = undef;
      $stage = $1;
      $lp    = $2;
      $stage_iter = 0;
      next;
    }

    if ($line =~ /^>>>\s*ITERATION\s+(\d+)/) {
      finalize_record(\@rows, $cur);
      $stage_iter++;
      $cur = new_record($meta, $stage, $lp, $1, $stage_iter, $line_no);
      $pending_shift = 0;
      next;
    }

    next if !$cur;

    if ($line =~ /^>>>\s*SHIFT\s+INCR\s+ARG\s+AVG\/SDEV\/MIN\/\s*$/) {
      $pending_shift = 1;
      next;
    }
    if ($pending_shift && $line =~ /^MAX:\s*($NUM_RE)\s+($NUM_RE)\s+($NUM_RE)\s+($NUM_RE)\s*$/) {
      @{$cur}{qw(shift_avg shift_sdev shift_min shift_max)} = ($1, $2, $3, $4);
      $pending_shift = 0;
      next;
    }
    $pending_shift = 0 if $pending_shift && $line !~ /^MAX:/;

    if ($line =~ /^>>>\s*ORIENTATION\s+OVERLAP:\s+($NUM_RE)/) {
      $cur->{overlap} = $1; next;
    }
    if ($line =~ /^>>>\s*STATE\s+OVERLAP:\s+($NUM_RE)/) {
      $cur->{state_overlap} = $1; next;
    }
    if ($line =~ /^>>>\s*%\s*PARTICLES\s+SAMPLED\s+THIS\s+ITERATION\s+($NUM_RE)(?:\s+\((\d+)\/(\d+)\))?/) {
      $cur->{pct_sampled} = $1;
      $cur->{n_sampled} = $2 if defined $2;
      $cur->{n_active}  = $3 if defined $3;
      next;
    }
    if ($line =~ /^>>>\s*%\s*PARTICLES\s+UPDATED\s+SO\s+FAR\s+($NUM_RE)/) {
      $cur->{pct_updated} = $1; next;
    }
    if ($line =~ /^>>>\s*%\s*PARTICLES\s+USED\s+FOR\s+UPDATING\s+THE\s+MODEL\s+($NUM_RE)/) {
      $cur->{pct_used_model} = $1; next;
    }
    if ($line =~ /^>>>\s*%\s*GREEDY\s+SEARCHES\s+($NUM_RE)/) {
      $cur->{pct_greedy} = $1; next;
    }

    if ($line =~ /^>>>\s*DIST\s+BTW\s+BEST\s+ORIS.*AVG\/SDEV\/MIN\/MAX:/) {
      my @q = parse_quad($line);
      @{$cur}{qw(dist_best_avg dist_best_sdev dist_best_min dist_best_max)} = @q if @q;
      next;
    }
    if ($line =~ /^>>>\s*IN-PLANE\s+DIST.*AVG\/SDEV\/MIN\/MAX:/) {
      my @q = parse_quad($line);
      @{$cur}{qw(inplane_avg inplane_sdev inplane_min inplane_max)} = @q if @q;
      next;
    }
    if ($line =~ /^>>>\s*SHIFT\s+INCR\s+ARG.*AVG\/SDEV\/MIN\/MAX:/) {
      my @q = parse_quad($line);
      @{$cur}{qw(shift_avg shift_sdev shift_min shift_max)} = @q if @q;
      next;
    }
    if ($line =~ /^>>>\s*%\s*SEARCH\s+SPACE\s+SCANNED.*AVG\/SDEV\/MIN\/MAX:/) {
      my @q = parse_quad($line);
      @{$cur}{qw(scanned_avg scanned_sdev scanned_min scanned_max)} = @q if @q;
      next;
    }
    if ($line =~ /^>>>\s*PARTICLE\s+WEIGHTS.*AVG\/SDEV\/MIN\/MAX:/) {
      my @q = parse_quad($line);
      @{$cur}{qw(pwt_avg pwt_sdev pwt_min pwt_max)} = @q if @q;
      next;
    }
    if ($line =~ /^>>>\s*MATCHING\s+LOW-PASS\s+LIMIT.*AVG\/SDEV\/MIN\/MAX:/) {
      my @q = parse_quad($line);
      @{$cur}{qw(matchlp_avg matchlp_sdev matchlp_min matchlp_max)} = @q if @q;
      next;
    }
    if ($line =~ /^>>>\s*ESTIMATED\s+LOW-PASS\s+LIMIT.*AVG\/SDEV\/MIN\/MAX:/) {
      my @q = parse_quad($line);
      @{$cur}{qw(estlp_avg estlp_sdev estlp_min estlp_max)} = @q if @q;
      next;
    }
    if ($line =~ /^>>>\s*RESOLUTION\s+\@\s*FSC=0\.143\s+STATE\s+(\d+)\s+AVG\/SDEV\/MIN\/MAX:/) {
      my $state = sprintf('%02d', $1);
      my @q = parse_quad($line);
      if (@q) {
        @{$cur}{"res_state_${state}_avg", "res_state_${state}_sdev", "res_state_${state}_min", "res_state_${state}_max"} = @q;
      }
      next;
    }
    if ($line =~ /^>>>\s*RESOLUTION\s+\@\s*FSC=0\.143.*AVG\/SDEV\/MIN\/MAX:/) {
      my @q = parse_quad($line);
      @{$cur}{qw(fsc143_avg fsc143_sdev fsc143_min fsc143_max)} = @q if @q;
      next;
    }
    if ($line =~ /^>>>\s*SCORE\s+\[0,1\].*AVG\/SDEV\/MIN\/MAX:/) {
      my @q = parse_quad($line);
      @{$cur}{qw(score_avg score_sdev score_min score_max)} = @q if @q;
      next;
    }

    if ($line =~ /^>>>\s*STATE\s+(\d+)\s+JOINT\s+DISTRIBUTION\s+OVERLAP:\s*($NUM_RE)\s+POPULATION:\s*(\d+)/) {
      my $state = sprintf('%02d', $1);
      $cur->{"state${state}_joint_overlap"} = $2;
      $cur->{"state${state}_population"} = $3;
      next;
    }

    if ($line =~ /^>>>\s*\|\s*(.*?)\s*\|\s*(.*?)\s*$/) {
      parse_setting_line($cur, $1, $2);
      next;
    }
    if ($line =~ /^>>>\s*REFINEMENT\s+MODE\s+IS\s+(\S+)/) {
      $cur->{refinement_mode} = $1; next;
    }
    if ($line =~ /^>>>\s*TRAILING\s+REC\s+UPDATE\s+FRACTION:\s*($NUM_RE)/) {
      $cur->{trailing_rec_update_fraction} = $1; next;
    }
    if ($line =~ /^>>>\s*CONVERGED:\s*(\S+)/) {
      $cur->{converged} = normalise_yes_no($1); next;
    }
  }

  finalize_record(\@rows, $cur);
  close $IN;

  for my $r (@rows) {
    $r->{normal_stop} = $meta->{normal_stop};
    $r->{execution_dir} = $meta->{execution_dir} if $r->{execution_dir} eq 'NA';
  }
  return (\@rows, $meta);
}

sub new_record {
  my ($meta, $stage, $lp, $iter, $stage_iter, $line_no) = @_;
  my ($ol, $fl) = stage_default_limits($stage);
  return {
    run_id        => $meta->{run_id},
    restart_id    => $meta->{restart_id},
    file          => $meta->{file},
    execution_dir => $meta->{execution_dir},
    normal_stop   => $meta->{normal_stop},
    stage         => $stage,
    lp            => $lp,
    iteration     => $iter,
    stage_iter    => $stage_iter,
    line_start    => $line_no,
    code_overlap_lim  => $ol,
    code_fracsrch_lim => $fl,
    code_would_converge => 'NA',
    converged     => 'NA',

    overlap => 'NA',
    overlap_gap => 'NA',
    state_overlap => 'NA',
    pct_sampled => 'NA',
    n_sampled => 'NA',
    n_active => 'NA',
    pct_updated => 'NA',
    pct_used_model => 'NA',
    pct_greedy => 'NA',

    dist_best_avg => 'NA', dist_best_sdev => 'NA', dist_best_min => 'NA', dist_best_max => 'NA',
    dist_ang_threshold => 'NA', dist_best_avg_minus_ang_threshold => 'NA', dist_best_avg_over_ang_threshold => 'NA',
    nonoverlap_fraction => 'NA', nonoverlap_mean_lb => 'NA', nonoverlap_mean_ub => 'NA',
    nonoverlap_mean_lb_minus_ang_threshold => 'NA',
    inplane_avg   => 'NA', inplane_sdev   => 'NA', inplane_min   => 'NA', inplane_max   => 'NA',
    shift_avg     => 'NA', shift_sdev     => 'NA', shift_min     => 'NA', shift_max     => 'NA',
    scanned_avg   => 'NA', scanned_sdev   => 'NA', scanned_min   => 'NA', scanned_max   => 'NA',
    fracsrch_avg  => 'NA', fracsrch_gap   => 'NA',
    pwt_avg       => 'NA', pwt_sdev       => 'NA', pwt_min       => 'NA', pwt_max       => 'NA',
    matchlp_avg   => 'NA', matchlp_sdev   => 'NA', matchlp_min   => 'NA', matchlp_max   => 'NA',
    estlp_avg     => 'NA', estlp_sdev     => 'NA', estlp_min     => 'NA', estlp_max     => 'NA',
    fsc143_avg    => 'NA', fsc143_sdev    => 'NA', fsc143_min    => 'NA', fsc143_max    => 'NA',
    score_avg     => 'NA', score_sdev     => 'NA', score_min     => 'NA', score_max     => 'NA',

    refinement_mode => 'NA',
    nspace => 'NA',
    gau_regularization => 'NA',
    ml_regularization => 'NA',
    ml_tau => 'NA',
    icm_regularization => 'NA',
    icm_lambda => 'NA',
    automasking => 'NA',
    filter_mode => 'NA',
    trailing_rec_update_fraction => 'NA',
    fillin_particle_sampling => 'NA',
  };
}

sub finalize_record {
  my ($rows, $r) = @_;
  return if !$r;
  $r->{fracsrch_avg} = $r->{scanned_avg};
  push @$rows, $r;
}

sub annotate_iterations {
  my ($rows, $ang) = @_;
  for my $r (@$rows) {
    my $ov = num($r->{overlap});
    my $fr = num($r->{fracsrch_avg});
    my $ol = num($r->{code_overlap_lim});
    my $fl = num($r->{code_fracsrch_lim});
    if (defined $ov && defined $ol) {
      $r->{overlap_gap} = fmt($ov - $ol);
    }
    if (defined $fr && defined $fl) {
      $r->{fracsrch_gap} = fmt($fr - $fl);
    }
    if (defined $ov && defined $fr && defined $ol && defined $fl) {
      $r->{code_would_converge} = ($ov > $ol && $fr > $fl) ? 'yes' : 'no';
    }

    my $dist = num($r->{dist_best_avg});
    $r->{dist_ang_threshold} = fmt($ang);
    if (defined $ov && defined $dist) {
      $r->{dist_best_avg_minus_ang_threshold} = fmt($dist - $ang);
      $r->{dist_best_avg_over_ang_threshold} = fmt($dist / $ang) if $ang > 0;
      my $non = 1.0 - $ov;
      $r->{nonoverlap_fraction} = fmt($non);
      if ($non > 1.0e-9) {
        my $lb = ($dist - $ov * $ang) / $non;
        $lb = $ang if $lb < $ang;
        my $ub = $dist / $non;
        $r->{nonoverlap_mean_lb} = fmt($lb);
        $r->{nonoverlap_mean_ub} = fmt($ub);
        $r->{nonoverlap_mean_lb_minus_ang_threshold} = fmt($lb - $ang);
      }
    }
  }
}

sub build_stage_summary {
  my ($rows) = @_;
  my %groups;
  for my $r (@$rows) {
    push @{$groups{join("\t", $r->{run_id}, $r->{stage})}}, $r;
  }

  my @out;
  for my $key (sort stage_group_sort keys %groups) {
    my @g = sort { num($a->{iteration}) <=> num($b->{iteration}) } @{$groups{$key}};
    my $first = $g[0];
    my $last  = $g[-1];
    my ($best_res_r)   = sort { num_or_inf($a->{fsc143_avg}) <=> num_or_inf($b->{fsc143_avg}) } @g;
    my ($best_score_r) = sort { num_or_ninf($b->{score_avg}) <=> num_or_ninf($a->{score_avg}) } @g;
    my $logged = first_row(\@g, sub { $_[0]->{converged} eq 'yes' });
    my $code   = first_row(\@g, sub { $_[0]->{code_would_converge} eq 'yes' });

    push @out, {
      run_id => $first->{run_id}, restart_id => $first->{restart_id}, file => $first->{file},
      stage => $first->{stage}, lp => $first->{lp},
      n_iters => scalar(@g),
      first_iteration => $first->{iteration}, last_iteration => $last->{iteration},
      first_logged_converged_iteration => $logged ? $logged->{iteration} : 'NA',
      first_code_converged_iteration   => $code   ? $code->{iteration}   : 'NA',
      code_overlap_lim => $first->{code_overlap_lim}, code_fracsrch_lim => $first->{code_fracsrch_lim},
      max_overlap => max_field(\@g, 'overlap'),
      last_overlap => $last->{overlap},
      max_fracsrch => max_field(\@g, 'fracsrch_avg'),
      last_fracsrch => $last->{fracsrch_avg},
      min_dist_best_avg => min_field(\@g, 'dist_best_avg'),
      last_dist_best_avg => $last->{dist_best_avg},
      min_fsc143_avg => $best_res_r->{fsc143_avg},
      last_fsc143_avg => $last->{fsc143_avg},
      max_score_avg => $best_score_r->{score_avg},
      last_score_avg => $last->{score_avg},
      overlap_gap_last => $last->{overlap_gap},
      fracsrch_gap_last => $last->{fracsrch_gap},
    };
  }
  return @out;
}

sub build_run_summary {
  my ($rows, $metas) = @_;
  my %groups;
  for my $r (@$rows) {
    push @{$groups{$r->{run_id}}}, $r;
  }
  my %meta_by_run = map { $_->{run_id} => $_ } @$metas;

  my @out;
  for my $run (sort run_sort keys %groups) {
    my @g = sort {
      num($a->{iteration}) <=> num($b->{iteration})
        || num($a->{stage}) <=> num($b->{stage})
    } @{$groups{$run}};
    my $first = $g[0];
    my $last = $g[-1];
    my %stages = map { $_->{stage} => 1 } @g;
    my ($best_res_r)   = sort { num_or_inf($a->{fsc143_avg}) <=> num_or_inf($b->{fsc143_avg}) } @g;
    my ($best_score_r) = sort { num_or_ninf($b->{score_avg}) <=> num_or_ninf($a->{score_avg}) } @g;
    my $logged = first_row(\@g, sub { $_[0]->{converged} eq 'yes' });
    my $code   = first_row(\@g, sub { $_[0]->{code_would_converge} eq 'yes' });
    my $meta   = $meta_by_run{$run} || {};

    push @out, {
      run_id => $run,
      restart_id => $first->{restart_id},
      file => $first->{file},
      execution_dir => $meta->{execution_dir} || $first->{execution_dir},
      normal_stop => $meta->{normal_stop} || 0,
      n_iters => scalar(@g),
      stages_seen => join(',', sort { num($a) <=> num($b) } keys %stages),
      first_stage => $first->{stage},
      last_stage => $last->{stage},
      first_iteration => $first->{iteration},
      last_iteration => $last->{iteration},
      first_logged_converged_iteration => $logged ? $logged->{iteration} : 'NA',
      first_code_converged_iteration   => $code   ? $code->{iteration}   : 'NA',
      final_stage => $last->{stage},
      final_overlap => $last->{overlap},
      final_fracsrch => $last->{fracsrch_avg},
      final_dist_best_avg => $last->{dist_best_avg},
      final_fsc143_avg => $last->{fsc143_avg},
      final_score_avg => $last->{score_avg},
      min_fsc143_avg => $best_res_r->{fsc143_avg},
      min_fsc143_stage => $best_res_r->{stage},
      min_fsc143_iteration => $best_res_r->{iteration},
      max_score_avg => $best_score_r->{score_avg},
      max_score_stage => $best_score_r->{stage},
      max_score_iteration => $best_score_r->{iteration},
    };
  }
  return @out;
}

sub build_convergence_scan {
  my ($rows, $overlap_limits, $frac_limits, $use_stage_frac) = @_;
  my %groups;
  for my $r (@$rows) {
    push @{$groups{join("\t", $r->{run_id}, $r->{stage})}}, $r;
  }
  my @out;
  for my $key (sort stage_group_sort keys %groups) {
    my @g = sort { num($a->{iteration}) <=> num($b->{iteration}) } @{$groups{$key}};
    my $first = $g[0];
    my $code_first = first_row(\@g, sub { $_[0]->{code_would_converge} eq 'yes' });
    my @flims = $use_stage_frac ? ($first->{code_fracsrch_lim}) : @$frac_limits;
    for my $ol (@$overlap_limits) {
      for my $fl (@flims) {
        my $hit = first_row(\@g, sub {
          my $ov = num($_[0]->{overlap});
          my $fr = num($_[0]->{fracsrch_avg});
          defined $ov && defined $fr && $ov > $ol && $fr > $fl;
        });
        my $saved = 'NA';
        if ($hit && $code_first) {
          $saved = num($code_first->{iteration}) - num($hit->{iteration});
        }
        push @out, {
          run_id => $first->{run_id}, restart_id => $first->{restart_id}, file => $first->{file},
          stage => $first->{stage}, lp => $first->{lp},
          overlap_limit => fmt($ol),
          fracsrch_limit => fmt($fl),
          first_converged_iteration => $hit ? $hit->{iteration} : 'NA',
          first_converged_stage_iter => $hit ? $hit->{stage_iter} : 'NA',
          code_first_converged_iteration => $code_first ? $code_first->{iteration} : 'NA',
          iters_saved_vs_code => $saved,
          n_iters => scalar(@g),
          max_overlap => max_field(\@g, 'overlap'),
          max_fracsrch => max_field(\@g, 'fracsrch_avg'),
        };
      }
    }
  }
  return @out;
}

sub build_near_misses {
  my ($rows) = @_;
  my @out;
  for my $r (@$rows) {
    my $ov = num($r->{overlap});
    my $fr = num($r->{fracsrch_avg});
    my $ol = num($r->{code_overlap_lim});
    my $fl = num($r->{code_fracsrch_lim});
    next unless defined $ov && defined $fr && defined $ol && defined $fl;
    next unless $fr > $fl && $ov <= $ol;
    push @out, {
      run_id => $r->{run_id}, restart_id => $r->{restart_id}, file => $r->{file},
      stage => $r->{stage}, lp => $r->{lp}, iteration => $r->{iteration}, stage_iter => $r->{stage_iter},
      overlap => $r->{overlap}, overlap_limit => $r->{code_overlap_lim},
      needed_extra_overlap => fmt($ol - $ov),
      fracsrch_avg => $r->{fracsrch_avg}, fracsrch_limit => $r->{code_fracsrch_lim},
      fracsrch_margin => $r->{fracsrch_gap},
      dist_best_avg => $r->{dist_best_avg}, dist_best_sdev => $r->{dist_best_sdev}, dist_best_max => $r->{dist_best_max},
      dist_ang_threshold => $r->{dist_ang_threshold},
      dist_best_avg_minus_ang_threshold => $r->{dist_best_avg_minus_ang_threshold},
      dist_best_avg_over_ang_threshold => $r->{dist_best_avg_over_ang_threshold},
      nonoverlap_fraction => $r->{nonoverlap_fraction},
      nonoverlap_mean_lb => $r->{nonoverlap_mean_lb},
      nonoverlap_mean_ub => $r->{nonoverlap_mean_ub},
      nonoverlap_mean_lb_minus_ang_threshold => $r->{nonoverlap_mean_lb_minus_ang_threshold},
      refinement_mode => $r->{refinement_mode}, nspace => $r->{nspace},
    };
  }
  return sort {
    num($a->{needed_extra_overlap}) <=> num($b->{needed_extra_overlap})
      || run_sort_key($a->{run_id}) cmp run_sort_key($b->{run_id})
      || num($a->{iteration}) <=> num($b->{iteration})
  } @out;
}

sub build_stage_averages_compact {
  my ($rows) = @_;
  my @summaries = run_stage_summaries($rows);
  my %by_stage;
  for my $s (@summaries) {
    push @{$by_stage{$s->{stage}}}, $s;
  }

  my @out;
  for my $stage (sort { num($a) <=> num($b) } keys %by_stage) {
    my @s = @{$by_stage{$stage}};
    push @out, {
      stage => $stage,
      lp_mean => fmt(mean_summary_field(\@s, 'lp')),
      n_runs => scalar(@s),
      n_iters_mean => fmt(mean_summary_field(\@s, 'n_iters')),
      n_iters_min => fmt(min_summary_field(\@s, 'n_iters')),
      n_iters_max => fmt(max_summary_field(\@s, 'n_iters')),
      res_start_mean => fmt(mean_summary_field(\@s, 'res_start')),
      res_end_mean => fmt(mean_summary_field(\@s, 'res_end')),
      res_best_mean => fmt(mean_summary_field(\@s, 'res_best')),
      res_gain_start_to_best_mean => fmt(mean_summary_field(\@s, 'res_gain_start_to_best')),
      res_end_minus_best_mean => fmt(mean_summary_field(\@s, 'res_end_minus_best')),
      dist_start_mean => fmt(mean_summary_field(\@s, 'dist_start')),
      dist_end_mean => fmt(mean_summary_field(\@s, 'dist_end')),
      dist_min_mean => fmt(mean_summary_field(\@s, 'dist_min')),
      dist_drop_start_to_end_mean => fmt(mean_summary_field(\@s, 'dist_drop_start_to_end')),
      dist_end_minus_min_mean => fmt(mean_summary_field(\@s, 'dist_end_minus_min')),
      overlap_end_mean => fmt(mean_summary_field(\@s, 'overlap_end')),
      fracsrch_end_mean => fmt(mean_summary_field(\@s, 'fracsrch_end')),
      score_end_mean => fmt(mean_summary_field(\@s, 'score_end')),
    };
  }
  return @out;
}

sub build_resolution_distance {
  my ($rows, $bin_width) = @_;
  my %bins;
  for my $r (@$rows) {
    my $res = num($r->{fsc143_avg});
    my $dist = num($r->{dist_best_avg});
    next if !defined($res) || !defined($dist);
    my $idx = int($res / $bin_width);
    push @{$bins{$idx}}, $r;
  }

  my @out;
  for my $idx (sort { $a <=> $b } keys %bins) {
    my @g = @{$bins{$idx}};
    my %runs = map { $_->{run_id} => 1 } @g;
    my %run_stages = map { $_->{run_id} . "\t" . $_->{stage} => 1 } @g;
    my @dist_vals = map { num($_->{dist_best_avg}) } @g;
    push @out, {
      resolution_bin_low => fmt($idx * $bin_width),
      resolution_bin_high => fmt(($idx + 1) * $bin_width),
      n_iterations => scalar(@g),
      n_runs => scalar(keys %runs),
      n_run_stages => scalar(keys %run_stages),
      resolution_mean => fmt(mean_field(\@g, 'fsc143_avg')),
      dist_best_oris_mean => fmt(mean_field(\@g, 'dist_best_avg')),
      dist_best_oris_sdev => fmt(sdev_vals(@dist_vals)),
      overlap_mean => fmt(mean_field(\@g, 'overlap')),
      fracsrch_mean => fmt(mean_field(\@g, 'fracsrch_avg')),
      score_mean => fmt(mean_field(\@g, 'score_avg')),
      stage_mean => fmt(mean_field(\@g, 'stage')),
    };
  }
  return @out;
}

sub build_halting_benefit {
  my ($rows, $benefit_min) = @_;
  my %by_stage;

  for my $g (run_stage_groups($rows)) {
    my @rows = @$g;
    my $first = $rows[0];
    my $last = $rows[-1];
    my $halt = first_row(\@rows, sub { $_[0]->{code_would_converge} eq 'yes' });
    next if !$halt;

    my @later = grep {
      defined(num($_->{stage_iter})) && defined(num($halt->{stage_iter}))
        && num($_->{stage_iter}) > num($halt->{stage_iter})
    } @rows;
    my $best_later = min_row_by_field(\@later, 'fsc143_avg');

    my $res_halt = num($halt->{fsc143_avg});
    my $res_later = $best_later ? num($best_later->{fsc143_avg}) : undef;
    my $res_gain = (defined($res_halt) && defined($res_later)) ? $res_halt - $res_later : undef;
    my $dist_halt = num($halt->{dist_best_avg});
    my $dist_final = num($last->{dist_best_avg});
    my $dist_drop = (defined($dist_halt) && defined($dist_final)) ? $dist_halt - $dist_final : undef;

    push @{$by_stage{$first->{stage}}}, {
      stage => $first->{stage},
      lp => num($first->{lp}),
      later_iters => scalar(@later),
      res_at_halt => $res_halt,
      best_later_res => $res_later,
      res_gain_later => $res_gain,
      later_res_improved => (defined($res_gain) && $res_gain > $benefit_min) ? 1 : 0,
      dist_at_halt => $dist_halt,
      final_dist => $dist_final,
      dist_drop_after_halt => $dist_drop,
      overlap_at_halt => num($halt->{overlap}),
      final_overlap => num($last->{overlap}),
    };
  }

  my @out;
  for my $stage (sort { num($a) <=> num($b) } keys %by_stage) {
    my @s = @{$by_stage{$stage}};
    my @with_later = grep { $_->{later_iters} > 0 } @s;
    my @with_gain = grep { defined($_->{res_gain_later}) } @with_later;
    my $pct_improved = @with_gain
      ? 100.0 * scalar(grep { $_->{later_res_improved} } @with_gain) / scalar(@with_gain)
      : undef;
    push @out, {
      stage => $stage,
      lp_mean => fmt(mean_summary_field(\@s, 'lp')),
      n_halt_points => scalar(@s),
      n_with_later_iters => scalar(@with_later),
      mean_later_iters_observed => fmt(mean_summary_field(\@s, 'later_iters')),
      mean_res_at_first_halt => fmt(mean_summary_field(\@s, 'res_at_halt')),
      mean_best_later_res => fmt(mean_summary_field(\@with_gain, 'best_later_res')),
      mean_res_gain_later => fmt(mean_summary_field(\@with_gain, 'res_gain_later')),
      pct_later_res_improved => fmt($pct_improved),
      mean_dist_at_first_halt => fmt(mean_summary_field(\@s, 'dist_at_halt')),
      mean_final_dist => fmt(mean_summary_field(\@s, 'final_dist')),
      mean_dist_drop_after_halt => fmt(mean_summary_field(\@s, 'dist_drop_after_halt')),
      mean_overlap_at_first_halt => fmt(mean_summary_field(\@s, 'overlap_at_halt')),
      mean_final_overlap => fmt(mean_summary_field(\@s, 'final_overlap')),
    };
  }
  return @out;
}

sub build_overall_summary {
  my ($rows, $metas, $halting) = @_;
  my @run_groups = run_groups($rows);
  my @final_rows = map { $_->[-1] } @run_groups;
  my @best_res_rows = grep { defined($_) } map { min_row_by_field($_, 'fsc143_avg') } @run_groups;
  my %stages = map { $_->{stage} => 1 } grep { defined(num($_->{stage})) } @$rows;
  my $normal_stops = scalar(grep { $_->{normal_stop} } @$metas);
  my $pair_stats = pair_stats($rows, 'fsc143_avg', 'dist_best_avg');
  my $best_stage = best_halting_stage($halting);

  my @out;
  push @out, row_metric('input_files', scalar(@$metas), 'ABINITIO3D_OUTPUT files parsed');
  push @out, row_metric('iteration_records', scalar(@$rows), 'logged iterations parsed');
  push @out, row_metric('normal_stop_runs', $normal_stops, 'runs ending with SIMPLE_ABINITIO3D NORMAL STOP');
  push @out, row_metric('normal_stop_fraction', fmt(scalar(@$metas) ? $normal_stops / scalar(@$metas) : undef), '1.0 means every parsed run reached normal stop');
  push @out, row_metric('stages_seen', join(',', sort { num($a) <=> num($b) } keys %stages), 'stage ids observed in the logs');
  push @out, row_metric('mean_iterations_per_run', fmt(mean_vals(map { scalar(@$_) } @run_groups)), 'averaged over runs');
  push @out, row_metric('mean_final_resolution', fmt(mean_field(\@final_rows, 'fsc143_avg')), 'final logged FSC=0.143 resolution; lower is better');
  push @out, row_metric('mean_best_resolution', fmt(mean_field(\@best_res_rows, 'fsc143_avg')), 'best logged FSC=0.143 resolution per run; lower is better');
  push @out, row_metric('mean_final_minus_best_resolution', fmt(mean_final_minus_best_resolution(\@run_groups)), 'near zero means runs ended at their best observed resolution');
  push @out, row_metric('mean_final_dist_best_oris', fmt(mean_field(\@final_rows, 'dist_best_avg')), 'final DIST BTW BEST ORIS average in degrees');
  push @out, row_metric('mean_final_overlap', fmt(mean_field(\@final_rows, 'overlap')), 'final binary orientation overlap');
  push @out, row_metric('mean_final_fracsrch', fmt(mean_field(\@final_rows, 'fracsrch_avg')), 'final percent search space scanned');
  push @out, row_metric('resolution_dist_pairs', $pair_stats->{n}, 'iterations with both resolution and DIST BTW BEST ORIS');
  push @out, row_metric('resolution_dist_pearson_r', fmt($pair_stats->{r}), 'correlation of resolution Angstrom vs orientation distance degrees');
  push @out, row_metric('resolution_dist_slope_deg_per_angstrom', fmt($pair_stats->{slope}), 'linear slope of distance against resolution');
  push @out, row_metric('best_stage_for_extra_iterations', $best_stage ? $best_stage->{stage} : 'NA', 'largest positive mean_res_gain_later in halting_benefit');
  push @out, row_metric('best_stage_mean_res_gain_after_halt', $best_stage ? $best_stage->{mean_res_gain_later} : 'NA', 'positive means later logged iterations improved resolution after first halt point');
  push @out, row_metric('best_stage_pct_later_res_improved', $best_stage ? $best_stage->{pct_later_res_improved} : 'NA', 'percent of observed later-iteration cases improving by --res-benefit-min');
  return @out;
}

sub run_stage_summaries {
  my ($rows) = @_;
  my @out;
  for my $g (run_stage_groups($rows)) {
    my @rows = @$g;
    my $first = $rows[0];
    my $last = $rows[-1];
    my $best_res = min_row_by_field(\@rows, 'fsc143_avg');
    my $min_dist = min_row_by_field(\@rows, 'dist_best_avg');
    my $res_start = num($first->{fsc143_avg});
    my $res_end = num($last->{fsc143_avg});
    my $res_best = $best_res ? num($best_res->{fsc143_avg}) : undef;
    my $dist_start = num($first->{dist_best_avg});
    my $dist_end = num($last->{dist_best_avg});
    my $dist_min = $min_dist ? num($min_dist->{dist_best_avg}) : undef;
    push @out, {
      run_id => $first->{run_id},
      stage => $first->{stage},
      lp => num($first->{lp}),
      n_iters => scalar(@rows),
      res_start => $res_start,
      res_end => $res_end,
      res_best => $res_best,
      res_gain_start_to_best => (defined($res_start) && defined($res_best)) ? $res_start - $res_best : undef,
      res_end_minus_best => (defined($res_end) && defined($res_best)) ? $res_end - $res_best : undef,
      dist_start => $dist_start,
      dist_end => $dist_end,
      dist_min => $dist_min,
      dist_drop_start_to_end => (defined($dist_start) && defined($dist_end)) ? $dist_start - $dist_end : undef,
      dist_end_minus_min => (defined($dist_end) && defined($dist_min)) ? $dist_end - $dist_min : undef,
      overlap_end => num($last->{overlap}),
      fracsrch_end => num($last->{fracsrch_avg}),
      score_end => num($last->{score_avg}),
    };
  }
  return @out;
}

sub run_stage_groups {
  my ($rows) = @_;
  my %groups;
  for my $r (@$rows) {
    next if !defined(num($r->{stage}));
    push @{$groups{join("\t", $r->{run_id}, $r->{stage})}}, $r;
  }
  return map {
    [ sort {
        num($a->{stage_iter}) <=> num($b->{stage_iter})
          || num($a->{iteration}) <=> num($b->{iteration})
      } @{$groups{$_}} ]
  } sort stage_group_sort keys %groups;
}

sub run_groups {
  my ($rows) = @_;
  my %groups;
  for my $r (@$rows) {
    push @{$groups{$r->{run_id}}}, $r;
  }
  return map {
    [ sort {
        num($a->{iteration}) <=> num($b->{iteration})
          || num($a->{stage}) <=> num($b->{stage})
      } @{$groups{$_}} ]
  } sort run_sort keys %groups;
}

sub min_row_by_field {
  my ($rows, $field) = @_;
  my $best;
  my $best_val;
  for my $r (@$rows) {
    my $v = num($r->{$field});
    next if !defined $v;
    if (!defined($best_val) || $v < $best_val) {
      $best = $r;
      $best_val = $v;
    }
  }
  return $best;
}

sub best_halting_stage {
  my ($rows) = @_;
  my $best;
  my $best_gain;
  for my $r (@$rows) {
    next if num($r->{n_with_later_iters}) <= 0;
    my $gain = num($r->{mean_res_gain_later});
    next if !defined $gain;
    if (!defined($best_gain) || $gain > $best_gain) {
      $best = $r;
      $best_gain = $gain;
    }
  }
  return $best;
}

sub mean_final_minus_best_resolution {
  my ($run_groups) = @_;
  my @gaps;
  for my $g (@$run_groups) {
    my $final = $g->[-1];
    my $best = min_row_by_field($g, 'fsc143_avg');
    next if !$best;
    my $final_res = num($final->{fsc143_avg});
    my $best_res = num($best->{fsc143_avg});
    push @gaps, $final_res - $best_res if defined($final_res) && defined($best_res);
  }
  return mean_vals(@gaps);
}

sub row_metric {
  my ($metric, $value, $note) = @_;
  return { metric => $metric, value => (defined($value) ? $value : 'NA'), note => $note };
}

sub parse_setting_line {
  my ($r, $key, $val) = @_;
  $key =~ s/^\s+|\s+$//g;
  $val =~ s/^\s+|\s+$//g;
  my $norm = uc($key);
  $norm =~ s/\s+/ /g;
  if ($norm eq 'REFINEMENT MODE') {
    $r->{refinement_mode} = $val;
  } elsif ($norm eq 'NSPACE PROJECTION DIRECTIONS') {
    $r->{nspace} = $val;
  } elsif ($norm eq 'GAU REGULARIZATION') {
    $r->{gau_regularization} = $val;
  } elsif ($norm eq 'ML REGULARIZATION') {
    if ($val =~ /^(on|off)\b/i) { $r->{ml_regularization} = lc($1); }
    if ($val =~ /TAU:\s*($NUM_RE)/i) { $r->{ml_tau} = $1; }
  } elsif ($norm eq 'ICM REGULARIZATION') {
    if ($val =~ /^(on|off)\b/i) { $r->{icm_regularization} = lc($1); }
    if ($val =~ /LAMBDA:\s*($NUM_RE)/i) { $r->{icm_lambda} = $1; }
  } elsif ($norm eq 'AUTOMASKING') {
    $r->{automasking} = $val;
  } elsif ($norm eq 'FILTER MODE') {
    $r->{filter_mode} = $val;
  } elsif ($norm eq 'TRAILING REC UPDATE FRACTION') {
    $r->{trailing_rec_update_fraction} = $val;
  } elsif ($norm eq 'FILLIN PARTICLE SAMPLING') {
    $r->{fillin_particle_sampling} = $val;
  }
}

sub parse_quad {
  my ($line) = @_;
  if ($line =~ /:\s*($NUM_RE)\s+($NUM_RE)\s+($NUM_RE)\s+($NUM_RE)\s*$/) {
    return ($1, $2, $3, $4);
  }
  return;
}

sub stage_default_limits {
  my ($stage) = @_;
  return ('NA', 'NA') if !defined($stage) || $stage eq 'NA';
  return (0.99, 99) if $stage <= 3;
  return (0.90, 90) if $stage <= 6;
  return (0.95, 90);
}

sub normalise_yes_no {
  my ($x) = @_;
  $x = lc($x // '');
  $x =~ s/\.//g;
  return 'yes' if $x eq 'yes' || $x eq 'true' || $x eq 't';
  return 'no'  if $x eq 'no'  || $x eq 'false' || $x eq 'f';
  return $x || 'NA';
}

sub find_input_files {
  my ($dir, $regex, $rec) = @_;
  my @files;
  my $re = qr/$regex/;
  if ($rec) {
    find({
      wanted => sub {
        return if !-f $_;
        push @files, $File::Find::name if basename($_) =~ $re;
      },
      no_chdir => 1,
    }, $dir);
  } else {
    opendir(my $DH, $dir) or die "ERROR: cannot open directory '$dir': $!\n";
    while (my $f = readdir($DH)) {
      next if $f =~ /^\./;
      my $p = "$dir/$f";
      push @files, $p if -f $p && $f =~ $re;
    }
    closedir($DH);
  }
  return sort { run_sort_key(basename($a)) cmp run_sort_key(basename($b)) } map { abs_path($_) || $_ } @files;
}

sub emit_table {
  my ($path, $header, $rows, $fmt) = @_;
  open(my $OUT, '>', $path) or die "ERROR: cannot write '$path': $!\n";
  print {$OUT} join_sep($fmt, @$header) . "\n";
  for my $r (@$rows) {
    print {$OUT} join_sep($fmt, map { exists $r->{$_} ? $r->{$_} : 'NA' } @$header) . "\n";
  }
  close $OUT;
}

sub join_sep {
  my ($fmt, @vals) = @_;
  if ($fmt eq 'csv') {
    return join(',', map { csv_escape($_) } @vals);
  }
  return join("\t", map { defined($_) ? $_ : 'NA' } @vals);
}

sub csv_escape {
  my ($v) = @_;
  $v = 'NA' if !defined $v;
  if ($v =~ /[",\n\r]/) {
    $v =~ s/"/""/g;
    return qq("$v");
  }
  return $v;
}

sub expand_tilde {
  my ($p) = @_;
  $p =~ s/^~(?=\/|$)/$ENV{HOME}/e;
  return $p;
}

sub parse_list {
  my ($s) = @_;
  return sort { $a <=> $b } grep { defined } map {
    s/^\s+|\s+$//gr =~ /^$NUM_RE$/ ? 0 + $_ : undef
  } split(/,/, $s);
}

sub num {
  my ($x) = @_;
  return undef if !defined($x) || $x eq 'NA' || $x eq '';
  return 0 + $x if $x =~ /^$NUM_RE$/;
  return undef;
}

sub num_or_inf {
  my ($x) = @_;
  my $n = num($x);
  return defined($n) ? $n : 9.0e99;
}

sub num_or_ninf {
  my ($x) = @_;
  my $n = num($x);
  return defined($n) ? $n : -9.0e99;
}

sub fmt {
  my ($x) = @_;
  return 'NA' if !defined $x;
  return sprintf('%.6g', $x);
}

sub mean_field {
  my ($rows, $field) = @_;
  return mean_vals(map { num($_->{$field}) } @$rows);
}

sub mean_summary_field {
  my ($rows, $field) = @_;
  return mean_vals(map { num($_->{$field}) } @$rows);
}

sub min_summary_field {
  my ($rows, $field) = @_;
  my @vals = grep { defined } map { num($_->{$field}) } @$rows;
  return undef if !@vals;
  my $best = $vals[0];
  for my $v (@vals) {
    $best = $v if $v < $best;
  }
  return $best;
}

sub max_summary_field {
  my ($rows, $field) = @_;
  my @vals = grep { defined } map { num($_->{$field}) } @$rows;
  return undef if !@vals;
  my $best = $vals[0];
  for my $v (@vals) {
    $best = $v if $v > $best;
  }
  return $best;
}

sub mean_vals {
  my @vals = grep { defined } @_;
  return undef if !@vals;
  my $sum = 0.0;
  $sum += $_ for @vals;
  return $sum / scalar(@vals);
}

sub sdev_vals {
  my @vals = grep { defined } @_;
  return undef if @vals < 2;
  my $mean = mean_vals(@vals);
  my $ss = 0.0;
  $ss += ($_ - $mean) * ($_ - $mean) for @vals;
  return sqrt($ss / (scalar(@vals) - 1));
}

sub pair_stats {
  my ($rows, $xfield, $yfield) = @_;
  my ($n, $sx, $sy, $sxx, $syy, $sxy) = (0, 0.0, 0.0, 0.0, 0.0, 0.0);
  for my $r (@$rows) {
    my $x = num($r->{$xfield});
    my $y = num($r->{$yfield});
    next if !defined($x) || !defined($y);
    $n++;
    $sx += $x;
    $sy += $y;
    $sxx += $x * $x;
    $syy += $y * $y;
    $sxy += $x * $y;
  }
  return { n => 0, r => undef, slope => undef } if $n < 2;
  my $cov = $n * $sxy - $sx * $sy;
  my $vx = $n * $sxx - $sx * $sx;
  my $vy = $n * $syy - $sy * $sy;
  my $r = ($vx > 0.0 && $vy > 0.0) ? $cov / sqrt($vx * $vy) : undef;
  my $slope = $vx > 0.0 ? $cov / $vx : undef;
  return { n => $n, r => $r, slope => $slope };
}

sub max_field {
  my ($rows, $field) = @_;
  my $best;
  for my $r (@$rows) {
    my $n = num($r->{$field});
    next if !defined $n;
    $best = $n if !defined($best) || $n > $best;
  }
  return defined($best) ? fmt($best) : 'NA';
}

sub min_field {
  my ($rows, $field) = @_;
  my $best;
  for my $r (@$rows) {
    my $n = num($r->{$field});
    next if !defined $n;
    $best = $n if !defined($best) || $n < $best;
  }
  return defined($best) ? fmt($best) : 'NA';
}

sub first_row {
  my ($rows, $pred) = @_;
  for my $r (@$rows) {
    return $r if $pred->($r);
  }
  return undef;
}

sub run_sort_key {
  my ($s) = @_;
  my ($n) = ($s =~ /RESTART(\d+)/);
  return defined($n) ? sprintf('%09d_%s', $n, $s) : $s;
}

sub run_sort {
  return run_sort_key($a) cmp run_sort_key($b);
}

sub stage_group_sort {
  my ($ra, $sa) = split(/\t/, $a);
  my ($rb, $sb) = split(/\t/, $b);
  return run_sort_key($ra) cmp run_sort_key($rb) || num($sa) <=> num($sb);
}

sub state_col_sort {
  my ($ap) = ($a =~ /(\d+)/);
  my ($bp) = ($b =~ /(\d+)/);
  return ($ap || 0) <=> ($bp || 0) || $a cmp $b;
}

sub usage {
  return <<"USAGE";
Usage:
  parse_abinitio_metrics_all.pl [options]
  perl scripts/parse_abinitio_metrics_all.pl [options]

Options:
  --help, -h                   Show this help and exit
  --indir DIR                  Directory containing ABINITIO3D_OUTPUT files
                               (default: current directory)
  --outdir DIR                 Directory for parsed tables
                               (default: ./abinitio_summary)
  --format tsv|csv              Output format for all tables (default: tsv)
  --pattern REGEX               Input basename regex (default: ^ABINITIO3D_OUTPUT)
  --recursive / --no-recursive  Recurse below --indir (default: no)
  --ang-threshold FLOAT         Projection-direction overlap angle threshold in degrees (default: 2.5)
  --res-bin-width FLOAT         Resolution bin width for resolution-vs-distance table
                                (default: 0.25 Angstrom)
  --res-benefit-min FLOAT       Resolution gain counted as a meaningful post-halt benefit
                                (default: 0.05 Angstrom)
  --overlap-limits LIST         Comma-separated convergence overlap limits for --detailed scan
                                (default: 0.80,0.85,0.90,0.92,0.95,0.98,0.99)
  --fracsrch-limits stage|LIST  Use stage policy or comma-separated search-space limits
                                for --detailed scan
                                (default: stage)
  --detailed / --no-detailed    Also write raw parsed iteration/run/stage/debug tables
                                (default: no)

Outputs:
  overall_summary.tsv/csv      compact key/value summary of the full batch
  stage_averages.tsv/csv       average start/end/best resolution and orientation distance by stage
  resolution_distance.tsv/csv  binned relationship between FSC resolution and DIST BTW BEST ORIS
  halting_benefit.tsv/csv      whether logged iterations after first convergence improved resolution

Detailed outputs with --detailed:
  iterations.tsv/csv           one row per parsed iteration
  stage_summary.tsv/csv        one row per run/stage
  run_summary.tsv/csv          one row per input file
  convergence_scan.tsv/csv     first convergence under candidate overlap/search limits
  near_misses.tsv/csv          rows where search passed but overlap missed current stage criterion

Examples:
  cd ~/abinitio_outputs
  /Users/elmlundho/src/SIMPLE/scripts/parse_abinitio_metrics_all.pl

  /Users/elmlundho/src/SIMPLE/scripts/parse_abinitio_metrics_all.pl \\
    --res-bin-width 0.5 --res-benefit-min 0.1

  /Users/elmlundho/src/SIMPLE/scripts/parse_abinitio_metrics_all.pl \\
    --detailed --overlap-limits 0.85,0.90,0.95 --fracsrch-limits stage

Notes:
  The script can scan overlap-limit convergence exactly from the logged metrics.
  DIST BTW BEST ORIS is the aggregate of the per-particle 'dist' angular
  distance stored by assign_ori. The script uses these aggregate angular-
  distance stats directly. It still cannot recompute exact overlap under
  alternate angular thresholds from logs alone because the per-particle
  distance CDF is not printed. The nonoverlap_* columns provide aggregate
  bounds for judging whether particles that failed the 2.5-degree binary test
  are probably close to the threshold.
USAGE
}
