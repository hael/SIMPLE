#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use File::Find qw(find);
use File::Spec;
use File::Basename qw(basename dirname);
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);

# One-stop generator:
# - Walk repo tree
# - Directory descriptions: read from <dir>/<dir>.inf (one-line description)
# - Fortran file descriptions: read from !@desc: within first N lines
# - Emit a single COMPLETE_CODEBASE_MAP.md with tree + descriptions
#
# Usage:
#   perl generate_complete_codebase_map.pl --root . --out docs/COMPLETE_CODEBASE_MAP.md
#
# Options:
#   --file-scan-lines 40
#   --ext f90 --ext F90 --ext f.90   (repeatable)
#   --include-nonfortran
#   --exclude-dir build --exclude-dir dist  (repeatable)

my $root               = '..';
my $out                = '../doc/code_overview/code_base_map.md';
my $file_scan_lines    = 40;
my @exts               = ('f90', 'F90');
my $include_nonfortran = 0;
my @exclude_dirs_user  = ('extlibs');
# Only include these top-level directories (relative to repo root)
my %ALLOWED_TOP_DIRS = map { $_ => 1 } qw(src production scripts);

GetOptions(
  'root=s'              => \$root,
  'out=s'               => \$out,
  'file-scan-lines=i'   => \$file_scan_lines,
  'ext=s@'              => \@exts,
  'include-nonfortran!' => \$include_nonfortran,
  'exclude-dir=s@'      => \@exclude_dirs_user,
) or die "Error parsing arguments\n";

$root = abs_path($root) // die "Cannot resolve root path\n";

my %DEFAULT_EXCLUDES = map { $_ => 1 } qw(
  .git .svn .hg
  __pycache__ .pytest_cache
  build dist
  .idea .vscode
);

for my $d (@exclude_dirs_user) {
  $DEFAULT_EXCLUDES{$d} = 1;
}

# Normalize extensions: store as lowercase strings without leading dot
my %EXTS = ();
for my $e (@exts) {
  $e =~ s/^\.+//;
  $EXTS{lc($e)} = 1;
}

sub is_blank {
  my ($s) = @_;
  return $s =~ /^\s*$/;
}

sub slurp_first_n_lines {
  my ($path, $n) = @_;
  open my $fh, '<:raw', $path or return ();
  my @lines;
  for (1..$n) {
    my $line = <$fh>;
    last unless defined $line;
    push @lines, $line;
  }
  close $fh;
  return @lines;
}

sub read_one_line_desc_file {
  my ($path) = @_;
  return '' unless -f $path;

  open my $fh, '<:raw', $path or return '';
  my @lines = <$fh>;
  close $fh;

  for my $ln (@lines) {
    chomp $ln;
    next if is_blank($ln);
    $ln =~ s/\s+/ /g;

    # If it uses !@descr: prefix, strip it
    if ($ln =~ /^\s*!\s*\@descr\s*:\s*(.+?)\s*$/i) {
      return $1;
    }
    return $ln; # plain text
  }
  return '';
}

sub read_dir_desc_from_inf {
  my ($dir) = @_;
  my $dname = basename($dir);
  $dname = '' unless defined $dname;

  # For root directory, basename might be something odd; handle gracefully.
  # We'll look for <root>/<basename(root)>.inf if basename is non-empty,
  # otherwise we skip.
  return '' if $dname eq '';

  my $inf = File::Spec->catfile($dir, $dname . '.inf');
  return read_one_line_desc_file($inf);
}

sub is_fortran_file {
  my ($path) = @_;
  my $lower = lc($path);
  for my $e (keys %EXTS) {
    my $suffix = '.' . $e;
    return 1 if substr($lower, -length($suffix)) eq $suffix;
  }
  return 0;
}

sub read_fortran_desc {
  my ($path) = @_;
  my @lines = slurp_first_n_lines($path, $file_scan_lines);
  for my $ln (@lines) {
    chomp $ln;
    if ($ln =~ /^\s*!\s*\@descr\s*:\s*(.+?)\s*$/i) {
      my $d = $1;
      $d =~ s/\s+/ /g;
      return $d;
    }
  }
  return '';
}

# Collect directories and files
my %dirs = ();            # dir => 1
my %files_by_dir = ();    # dir => [files...]

find(
  {
    no_chdir => 1,
    wanted => sub {
      my $path = $File::Find::name;

      # Normalize path
      my $abs = abs_path($path) // return;

      # Skip everything outside allowed top-level dirs
      if ($abs ne $root) {
        my $rel = File::Spec->abs2rel($abs, $root);
        $rel =~ s{\\}{/}g;

        my ($top) = split('/', $rel, 2);
        unless ($ALLOWED_TOP_DIRS{$top}) {
          if (-d $abs) {
            $File::Find::prune = 1;
          }
          return;
        }
      }

      # Directory handling
      if (-d $abs) {
        my $name = basename($abs);
        if ($DEFAULT_EXCLUDES{$name}) {
          $File::Find::prune = 1;
          return;
        }
        $dirs{$abs} = 1;
        return;
      }

      # File handling
      return unless -f $abs;

      my $dir = dirname($abs);
      if ($include_nonfortran) {
        push @{ $files_by_dir{$dir} }, $abs;
      } else {
        return unless is_fortran_file($abs);
        push @{ $files_by_dir{$dir} }, $abs;
      }
    },
  },
  $root
);

$dirs{$root} = 1;

# Sort dirs by depth then name
my @dir_list = sort {
  my $da = () = File::Spec->splitdir($a);
  my $db = () = File::Spec->splitdir($b);
  $da <=> $db || lc($a) cmp lc($b)
} keys %dirs;

# Sort files per dir
for my $d (keys %files_by_dir) {
  my @sorted = sort { lc(basename($a)) cmp lc(basename($b)) } @{ $files_by_dir{$d} };
  $files_by_dir{$d} = \@sorted;
}

# Build children map
my %children = ();
for my $d (@dir_list) {
  next if $d eq $root;
  my $parent = dirname($d);
  push @{ $children{$parent} }, $d;
}
for my $p (keys %children) {
  my @sorted = sort { lc(basename($a)) cmp lc(basename($b)) } @{ $children{$p} };
  $children{$p} = \@sorted;
}

# Ensure output directory exists
my $out_dir = dirname($out);
if ($out_dir && $out_dir ne '.' && !-d $out_dir) {
  my @parts = File::Spec->splitdir($out_dir);
  my $acc = '';
  for my $p (@parts) {
    $acc = $acc eq '' ? $p : File::Spec->catdir($acc, $p);
    mkdir $acc unless -d $acc;
  }
}

open my $out_fh, '>:encoding(UTF-8)', $out or die "Cannot write $out: $!\n";

print $out_fh "# SIMPLE Codebase Map\n\n";
#print $out_fh "Tree view with:\n\n";
#print $out_fh "- one-line directory descriptions from `<dir>/<dir>.inf` (first non-empty line)\n";
#print $out_fh "- one-line Fortran file descriptions from `!\@desc:` near the top of each file\n\n";

# Iterative DFS stack: [dir, depth]
my @stack = ([$root, 0]);
my %visited = ();

while (@stack) {
  my $item = pop @stack;
  my ($dir, $depth) = @$item;

  next if $visited{$dir}++;
  my $indent = '  ' x $depth;

  my $ddesc = read_dir_desc_from_inf($dir);
  my $dname = basename($dir);
  $dname = $dname eq '' ? '/' : $dname;

  if ($ddesc ne '') {
    print $out_fh "${indent}- **$dname/** — $ddesc\n";
  } else {
    print $out_fh "${indent}- **$dname/**\n";
  }

  # Files
  if (exists $files_by_dir{$dir}) {
    for my $f (@{ $files_by_dir{$dir} }) {
      my $fname = basename($f);
      my $fdesc = $include_nonfortran ? '' : read_fortran_desc($f);

      if ($fdesc ne '') {
        print $out_fh "${indent}  - `$fname` — $fdesc\n";
      } else {
        print $out_fh "${indent}  - `$fname`\n";
      }
    }
  }

  # Subdirectories: push reverse for stable DFS
  my $subs = $children{$dir} // [];
  for my $sd (reverse @$subs) {
    push @stack, [$sd, $depth + 1];
  }
}

close $out_fh;

my $dir_count  = scalar(keys %dirs);
my $file_count = 0;
$file_count += scalar(@{ $files_by_dir{$_} }) for keys %files_by_dir;

print "Wrote: $out\n";
print "Directories: $dir_count\n";
print "Listed files: $file_count\n";

