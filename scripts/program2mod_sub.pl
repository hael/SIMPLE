#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

# Convert Fortran PROGRAM blocks into contained SUBROUTINE blocks
# - Removes IMPLICIT NONE
# - Emits "subroutine name" (no parentheses)
# - Strips leading "simple_" from program names
# - Uses 4-space indentation for the subroutine body
# - Handles multiple PROGRAM blocks in one file/stream

my $prefix     = "";
my $indent     = "  ";     # indent for CONTAINS level
my $bodyindent = "    ";   # always 4 spaces
my $help       = 0;

GetOptions(
  "prefix=s" => \$prefix,
  "indent=s" => \$indent,
  "help"     => \$help,
) or die "Use --help\n";

if ($help) {
  print <<"HELP";
Usage:
  prog2contains.pl [--prefix=ut_] [--indent="  "] [file.f90]
  cat file.f90 | prog2contains.pl > contains.inc

Behavior:
  - program simple_foo   -> subroutine foo
  - end program          -> end subroutine foo
  - removes 'implicit none'
  - body indented by 4 spaces
  - supports multiple PROGRAM blocks
HELP
  exit 0;
}

my @lines;
if (@ARGV) {
  my $file = shift @ARGV;
  open my $fh, "<", $file or die "Cannot open $file: $!\n";
  @lines = <$fh>;
  close $fh;
} else {
  @lines = <STDIN>;
}

my $in_prog   = 0;
my $prog_name = "";
my @buf;

sub strip_simple_prefix {
  my ($name) = @_;
  $name =~ s/^simple_//i;
  return $name;
}

sub flush_prog {
  my ($name, $bufref) = @_;
  return if !$name;

  my $clean = strip_simple_prefix($name);
  my $sname = $prefix . $clean;

  my $head = $indent;
  my $body = $indent . $bodyindent;

  print $head . "subroutine $sname\n";

  for my $l (@$bufref) {

    next if $l =~ /^\s*program\b/i;
    next if $l =~ /^\s*implicit\s+none\b/i;

    if ($l =~ /^\s*end\s+program\b/i) {
      print $head . "end subroutine $sname\n\n";
      return;
    }

    if ($l =~ /^\s*$/) {
      print "\n";
    } else {
      print $body . $l;
    }
  }

  print $head . "end subroutine $sname\n\n";
}

for my $line (@lines) {

  if ($line =~ /^\s*program\s+([A-Za-z_]\w*)\b/i) {

    if ($in_prog) {
      flush_prog($prog_name, \@buf);
      @buf = ();
    }

    $in_prog   = 1;
    $prog_name = $1;
    push @buf, $line;
    next;
  }

  if ($in_prog) {
    push @buf, $line;

    if ($line =~ /^\s*end\s+program\b/i) {
      flush_prog($prog_name, \@buf);
      @buf = ();
      $in_prog   = 0;
      $prog_name = "";
    }
  }
}

flush_prog($prog_name, \@buf) if $in_prog;

