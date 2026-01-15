#!/usr/bin/env perl
use strict;
use warnings;
use File::Find;

my $root_dir = shift // '.';

# Fortran source extensions
my %is_fortran = map { $_ => 1 } qw(.f .f08 .f90 .f95 .F .F90 .F95);

# Regex for lines to remove
my $remove_re = qr/^\s*use\s+(?:simple_defs
                           |simple_type_defs
                           |simple_defs_fname
                           |simple_defs_stream
                           |simple_defs_conv)\b/ix;

# Regex to detect the trigger
my $trigger_re = qr/^\s*use\s+simple_core_module_api\b/i;

find(
    {
        wanted => sub {
            return unless -f $_;

            my ($ext) = $_ =~ /(\.[^.]+)$/;
            return unless $ext && $is_fortran{$ext};

            open my $fh, '<', $_ or die "Cannot open $_: $!";
            my @lines = <$fh>;
            close $fh;

            # Check if trigger is present
            my $has_trigger = 0;
            for (@lines) {
                if (/$trigger_re/) {
                    $has_trigger = 1;
                    last;
                }
            }
            return unless $has_trigger;

            # Filter lines
            my @new_lines;
            my $changed = 0;
            for my $line (@lines) {
                if ($line =~ $remove_re) {
                    $changed = 1;
                    next;
                }
                push @new_lines, $line;
            }

            return unless $changed;

            # Backup original
            rename $_, "$_.bak" or die "Backup failed for $_: $!";

            open my $out, '>', $_ or die "Cannot write $_: $!";
            print {$out} @new_lines;
            close $out;

            print "Updated $_\n";
        },
        no_chdir => 1,
    },
    $root_dir
);
