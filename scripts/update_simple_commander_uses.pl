#!/usr/bin/env perl
use strict;
use warnings;
use File::Find;

my $root_dir = shift // '.';

# Files must be named simple_commanders*
sub is_target_file {
    my ($name) = @_;
    return $name =~ /^simple_commanders.*\.(f|f90|f95)$/i;
}

# Trigger if ANY of these appear
my $trigger_re = qr/^\s*use\s+(?:simple_core_module_api
                             |simple_commanders_euclid
                             |simple_default_clines
                             |simple_nice
                             |simple_qsys_funs
                             |simple_binoris_io
                             |simple_builder)\b/ix;

# Remove ANY form of these use-statements
my $remove_re = $trigger_re;

find(
    {
        wanted => sub {
            return unless -f $_;
            return unless is_target_file($_);

            open my $fh, '<', $_ or die "Cannot open $_: $!";
            my @lines = <$fh>;
            close $fh;

            # Check trigger
            my $triggered = 0;
            for (@lines) {
                if (/$trigger_re/) {
                    $triggered = 1;
                    last;
                }
            }
            return unless $triggered;

            my @out;
            my $changed  = 0;
            my $inserted = 0;

            for my $line (@lines) {

                # Insert single replacement after module declaration
                if (!$inserted && $line =~ /^\s*module\s+simple_commanders\b/i) {
                    push @out, $line;
                    push @out, "use simple_commander_module_api\n";
                    $inserted = 1;
                    $changed  = 1;
                    next;
                }

                # Remove old use statements
                if ($line =~ $remove_re) {
                    $changed = 1;
                    next;
                }

                push @out, $line;
            }

            return unless $changed;

            rename $_, "$_.bak" or die "Backup failed for $_: $!";
            open my $out_fh, '>', $_ or die "Cannot write $_: $!";
            print {$out_fh} @out;
            close $out_fh;

            print "Updated $_\n";
        },
        no_chdir => 1,
    },
    $root_dir
);
