# perl move_constructors.pl \
#   --src simple_user_interface.f90 \
#   --map /mnt/data/instance_assigments.txt \
#   --ctors /mnt/data/ui_procedures.txt \
#   --dry-run
#
# perl move_constructors.pl \
#   --src simple_user_interface.f90 \
#   --map /mnt/data/instance_assigments.txt \
#   --ctors /mnt/data/ui_procedures.txt \
#   --apply
#
#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Copy qw(copy);

# --------------------------
# CLI
# --------------------------
my $src      = 'simple_user_interface.f90';
my $mapfile  = 'instance_assigments.txt';   # note spelling matches your upload path/name
my $ctors    = 'ui_procedures.txt';
my $dry_run  = 0;
my $apply    = 0;

GetOptions(
  'src=s'   => \$src,
  'map=s'   => \$mapfile,
  'ctors=s' => \$ctors,
  'dry-run' => \$dry_run,
  'apply'   => \$apply,
) or die "Usage: $0 --src simple_user_interface.f90 --map instance_assigments.txt --ctors ui_procedures.txt [--dry-run|--apply]\n";

if (!$dry_run && !$apply) {
  $dry_run = 1; # default safe behavior
}

# --------------------------
# Helpers
# --------------------------
sub slurp {
  my ($path) = @_;
  open my $fh, '<', $path or die "Cannot open $path: $!";
  local $/;
  my $txt = <$fh>;
  close $fh;
  return $txt;
}

sub write_file {
  my ($path, $txt) = @_;
  open my $fh, '>', $path or die "Cannot write $path: $!";
  print $fh $txt;
  close $fh;
}

sub backup_once {
  my ($path) = @_;
  my $bak = "$path.bak";
  return if -e $bak;
  copy($path, $bak) or die "Backup failed $path -> $bak: $!";
}

sub norm {
  my ($s) = @_;
  $s =~ s/^\s+|\s+$//g;
  return lc($s);
}

# --------------------------
# Parse mapping: proc -> file
# format like:
#   simple_ui_api_image.f90:
#   binarize
#   convert
#   ...
# --------------------------
my %proc2file;
{
  open my $fh, '<', $mapfile or die "Cannot open mapping file $mapfile: $!";
  my $current;
  while (my $line = <$fh>) {
    chomp $line;
    $line =~ s/\r$//;
    next if $line =~ /^\s*$/;

    if ($line =~ /^\s*([^\s:]+\.f90)\s*:\s*$/i) {
      $current = $1;
      next;
    }
    next unless defined $current;

    my $p = norm($line);
    next if $p eq '';
    $proc2file{$p} = $current;
  }
  close $fh;
}

# --------------------------
# Read constructors list: new_xxx
# --------------------------
my @ctor_names;
{
  open my $fh, '<', $ctors or die "Cannot open ctors file $ctors: $!";
  while (my $line = <$fh>) {
    chomp $line;
    $line =~ s/\r$//;
    $line =~ s/^\s+|\s+$//g;
    next if $line eq '';
    push @ctor_names, $line;
  }
  close $fh;
}

# --------------------------
# Extract procedure blocks from source
# We do a case-insensitive scan for:
#   subroutine NAME ...  / function NAME ...
# and end at:
#   end subroutine [NAME] / end function [NAME]
# --------------------------
my $src_txt = slurp($src);

# We'll collect extracted blocks keyed by constructor name (lowercase)
my %extracted_block;
my %extracted_kind;  # 'subroutine' or 'function'

for my $ctor (@ctor_names) {
  my $ctor_lc = norm($ctor);

  # Regex for start:
  # allow leading attributes e.g. "pure subroutine", "elemental function", etc.
  # capture kind and name
  my $start_re = qr{
    ^[ \t]*
    (?:(?:pure|impure|elemental|recursive|module)\s+)*   # optional prefixes
    (subroutine|function)
    \s+
    \Q$ctor\E
    \b
  }imx;

  # Find start position
  if ($src_txt !~ /$start_re/) {
    warn "WARN: constructor not found in $src: $ctor\n";
    next;
  }

  my $kind = lc($1);
  my $start_pos = $-[0];

  # Find the matching end line after start_pos
  # Prefer "end subroutine name" but allow "end subroutine"
  my $after = substr($src_txt, $start_pos);

  my $end_re_named = qr{
    ^[ \t]*end[ \t]+$kind[ \t]+\Q$ctor\E\b.*\n
  }imx;

  my $end_re_unnamed = qr{
    ^[ \t]*end[ \t]+$kind\b.*\n
  }imx;

  my ($end_pos_rel, $end_len);
  if ($after =~ /$end_re_named/) {
    $end_pos_rel = $-[0];
    $end_len     = $+[0] - $-[0];
  } elsif ($after =~ /$end_re_unnamed/) {
    $end_pos_rel = $-[0];
    $end_len     = $+[0] - $-[0];
  } else {
    die "ERROR: Found start of $ctor ($kind) but could not find end $kind in $src\n";
  }

  my $block_len = $end_pos_rel + $end_len;
  my $block = substr($after, 0, $block_len);

  $extracted_block{$ctor_lc} = $block;
  $extracted_kind{$ctor_lc}  = $kind;

  # Remove from source by replacing with blank line(s)
  substr($src_txt, $start_pos, $block_len) = "\n! moved: $ctor\n";
}

# --------------------------
# Route blocks to target files
# ctor new_X goes to file where X lives
# --------------------------
my %file2blocks;
for my $ctor (@ctor_names) {
  my $ctor_lc = norm($ctor);
  next unless exists $extracted_block{$ctor_lc};

  (my $base = $ctor_lc) =~ s/^new_//;  # map new_stack -> stack

  my $tgt = $proc2file{$base};
  if (!$tgt) {
    warn "WARN: No mapping for base proc '$base' (from $ctor). Leaving extracted block unused.\n";
    next;
  }

  push @{ $file2blocks{$tgt} }, $extracted_block{$ctor_lc};
}

# --------------------------
# Apply or dry-run
# --------------------------
print "Source file: $src\n";
print "Mapping file: $mapfile\n";
print "Constructors file: $ctors\n";
print($dry_run ? "Mode: DRY-RUN (no files will be modified)\n" : "Mode: APPLY (files will be modified; .bak created)\n");

# Update target files
for my $tgt (sort keys %file2blocks) {
  my $tgt_txt = slurp($tgt);
  my $insert = "";

  $insert .= "\n! ============================================================\n";
  $insert .= "! Constructors moved from $src\n";
  $insert .= "! ============================================================\n";
  for my $blk (@{ $file2blocks{$tgt} }) {
    $insert .= "\n$blk\n";
  }

  # Ensure there's a CONTAINS. If none, inject before end module.
  if ($tgt_txt !~ /^\s*contains\b/im) {
    if ($tgt_txt =~ /^\s*end\s+module\b/im) {
      $tgt_txt =~ s/^\s*end\s+module\b/contains\n\nend module/im;
    } else {
      warn "WARN: $tgt has no CONTAINS and no END MODULE; inserting anyway at end of file.\n";
      $tgt_txt .= "\ncontains\n";
    }
  }

  # Insert before END MODULE
  if ($tgt_txt =~ /^\s*end\s+module\b/im) {
    $tgt_txt =~ s/^\s*end\s+module\b/$insert\nend module/im;
  } else {
    warn "WARN: Could not find END MODULE in $tgt; appending at end.\n";
    $tgt_txt .= $insert;
  }

  print "Would move into $tgt: " . scalar(@{ $file2blocks{$tgt} }) . " procedure(s)\n";

  if ($apply) {
    backup_once($tgt);
    write_file($tgt, $tgt_txt);
  }
}

# Update source file (with removed blocks)
print "Would update source $src (constructors removed/replaced with '! moved: ...')\n";
if ($apply) {
  backup_once($src);
  write_file($src, $src_txt);
}

print "Done.\n";
