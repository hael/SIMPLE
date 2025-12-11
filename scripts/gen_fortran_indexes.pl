#!/usr/bin/env perl
use strict;
use warnings;
use File::Find;
use Getopt::Long;
use File::Path qw(make_path);
use File::Spec;

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

my @FORTRAN_EXT = qw(.f .F .f90 .F90 .f95 .F95);

# ------------------------------------------------------------
# CLI args
# ------------------------------------------------------------

my $root_dir;
my $out_dir = '.';

GetOptions(
    'root=s' => \$root_dir,
    'out=s'  => \$out_dir,
) or die "Usage: $0 --root /path/to/src [--out docs]\n";

die "Usage: $0 --root /path/to/src [--out docs]\n"
    unless defined $root_dir && -d $root_dir;

make_path($out_dir) unless -d $out_dir;

# ------------------------------------------------------------
# Data structures
# ------------------------------------------------------------
# %modules:
#   key   = module name (lowercase)
#   value = {
#       name        => original name as seen,
#       files       => { file_path => 1, ... },
#       procs       => [ { name, kind, file, line, visibility }, ... ],
#       visibility  => { symbol_name => 'public'|'private' },
#       default_vis => 'public'|'private'|undef
#   }

my %modules;

# List of all procedures for global index
# @procedures: { name, kind, module, file, line, visibility }

my @procedures;

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

sub is_fortran_file {
    my ($file) = @_;
    foreach my $ext (@FORTRAN_EXT) {
        return 1 if $file =~ /\Q$ext\E$/;
    }
    return 0;
}

sub normalize_name {
    my ($name) = @_;
    $name =~ s/^\s+//;
    $name =~ s/\s+$//;
    return lc $name;
}

sub strip_inline_comment {
    my ($line) = @_;
    # crude: strip anything after ! unless it's in a string (ignore that complexity)
    $line =~ s/!.*$//;
    return $line;
}

# ------------------------------------------------------------
# Parse a single file
# ------------------------------------------------------------

sub process_fortran_file {
    my ($file) = @_;

    open my $fh, '<', $file or do {
        warn "Cannot open $file: $!";
        return;
    };

    my $current_module;
    my $line_no = 0;

    while (my $line = <$fh>) {
        $line_no++;
        my $raw_line = $line;
        $line = strip_inline_comment($line);
        next if $line =~ /^\s*$/; # skip empty

        # Detect "module" (but not "module procedure")
        if ($line =~ /^\s*module\s+([a-zA-Z0-9_]+)\b(?!\s*procedure)/i) {
            my $modname_orig = $1;
            my $modkey = normalize_name($modname_orig);

            $modules{$modkey} ||= {
                name        => $modname_orig,
                files       => {},
                procs       => [],
                visibility  => {},
                default_vis => undef, # unknown until we see PRIVATE without list
            };
            $modules{$modkey}{files}{$file} = 1;
            $current_module = $modkey;
            next;
        }

        # End module
        if ($line =~ /^\s*end\s+module\b/i) {
            $current_module = undef;
            next;
        }

        # If inside a module, try to pick up PUBLIC/PRIVATE statements
        if (defined $current_module) {
            my $mod = $modules{$current_module};

            # PRIVATE with no list => default private
            if ($line =~ /^\s*private\s*$/i) {
                $mod->{default_vis} = 'private';
                next;
            }

            # PUBLIC with no list => default public
            if ($line =~ /^\s*public\s*$/i) {
                $mod->{default_vis} = 'public';
                next;
            }

            # PUBLIC :: a, b
            if ($line =~ /^\s*public\s*::\s*(.+)$/i) {
                my $list = $1;
                my @names = map { normalize_name($_) } split /,/, $list;
                $mod->{visibility}{$_} = 'public' for @names;
                next;
            }

            # PRIVATE :: a, b
            if ($line =~ /^\s*private\s*::\s*(.+)$/i) {
                my $list = $1;
                my @names = map { normalize_name($_) } split /,/, $list;
                $mod->{visibility}{$_} = 'private' for @names;
                next;
            }
        }

        # Detect subroutine
        # allow prefixes: recursive, pure, elemental, module ...
        if ($line =~ /^\s*(?:(?:recursive|pure|elemental|impure|module)\s+)*subroutine\s+([a-zA-Z0-9_]+)\b/i) {
            my $subname_orig = $1;
            my $subkey = normalize_name($subname_orig);

            my $kind = 'subroutine';
            my $mod_name = defined $current_module ? $modules{$current_module}{name} : '';
            my $vis = 'unknown';

            if (defined $current_module) {
                my $mod = $modules{$current_module};
                if (exists $mod->{visibility}{$subkey}) {
                    $vis = $mod->{visibility}{$subkey};
                } elsif (defined $mod->{default_vis}) {
                    $vis = $mod->{default_vis};
                } else {
                    # Fortran default: public
                    $vis = 'public';
                }
            }

            my $proc = {
                name       => $subname_orig,
                kind       => $kind,
                module     => $mod_name,
                file       => $file,
                line       => $line_no,
                visibility => $vis,
            };

            push @procedures, $proc;

            if (defined $current_module) {
                push @{$modules{$current_module}{procs}}, $proc;
            }

            next;
        }

        # Detect function
        # allow prefixes: recursive, pure, elemental, impure, module ...
        if ($line =~ /^\s*(?:(?:recursive|pure|elemental|impure|module)\s+)*function\s+([a-zA-Z0-9_]+)\b/i) {
            my $fname_orig = $1;
            my $fkey = normalize_name($fname_orig);

            my $kind = 'function';
            my $mod_name = defined $current_module ? $modules{$current_module}{name} : '';
            my $vis = 'unknown';

            if (defined $current_module) {
                my $mod = $modules{$current_module};
                if (exists $mod->{visibility}{$fkey}) {
                    $vis = $mod->{visibility}{$fkey};
                } elsif (defined $mod->{default_vis}) {
                    $vis = $mod->{default_vis};
                } else {
                    $vis = 'public';
                }
            }

            my $proc = {
                name       => $fname_orig,
                kind       => $kind,
                module     => $mod_name,
                file       => $file,
                line       => $line_no,
                visibility => $vis,
            };

            push @procedures, $proc;

            if (defined $current_module) {
                push @{$modules{$current_module}{procs}}, $proc;
            }

            next;
        }
    }

    close $fh;
}

# ------------------------------------------------------------
# Walk the tree
# ------------------------------------------------------------

my @files;

find(
    {
        wanted => sub {
            return unless -f $_;
            return unless is_fortran_file($_);
            my $full = $File::Find::name;
            push @files, File::Spec->rel2abs($full);
        },
        no_chdir => 1,
    },
    $root_dir
);

foreach my $f (@files) {
    process_fortran_file($f);
}

# ------------------------------------------------------------
# Emit modules.md
# ------------------------------------------------------------

sub write_modules_md {
    my ($path) = @_;

    open my $out, '>', $path or die "Cannot write $path: $!";

    print $out "# Module Index\n\n";
    print $out "| Module | Files | # Procedures |\n";
    print $out "|--------|-------|--------------|\n";

    foreach my $modkey (sort keys %modules) {
        my $m = $modules{$modkey};
        my $name = $m->{name};
        my @files = sort keys %{$m->{files}};
        my $files_str = join("<br>", @files);
        my $n_procs = scalar @{$m->{procs}};

        print $out "| `$name` | $files_str | $n_procs |\n";
    }

    close $out;
}

# ------------------------------------------------------------
# Emit api_index.md
# ------------------------------------------------------------

sub write_api_index_md {
    my ($path) = @_;

    open my $out, '>', $path or die "Cannot write $path: $!";

    print $out "# API Index\n\n";
    print $out "| Procedure | Kind | Module | File | Line | Visibility |\n";
    print $out "|-----------|------|--------|------|------|------------|\n";

    # sort by module, then name
    my @sorted = sort {
           ($a->{module} cmp $b->{module})
        || ($a->{name}   cmp $b->{name})
    } @procedures;

    foreach my $p (@sorted) {
        my $pname = $p->{name};
        my $kind  = $p->{kind};
        my $mod   = $p->{module} || '';
        my $file  = $p->{file};
        my $line  = $p->{line};
        my $vis   = $p->{visibility};

        print $out "| `$pname` | $kind | `$mod` | $file | $line | $vis |\n";
    }

    close $out;
}

# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------

my $modules_md   = File::Spec->catfile($out_dir, 'modules.md');
my $api_index_md = File::Spec->catfile($out_dir, 'api_index.md');

write_modules_md($modules_md);
write_api_index_md($api_index_md);

print "Wrote:\n  $modules_md\n  $api_index_md\n";
