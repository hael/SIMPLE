#!/usr/bin/env perl
use strict;
use warnings;
use File::Find;
use Getopt::Long;
use File::Path qw(make_path);
use File::Spec;
use Cwd qw(abs_path);

# ------------------------------------------------------------
# Config
# ------------------------------------------------------------

# Only these statement starters can affect the index. Include every return type
# and qualifier accepted by procedure_prefix so typed functions still reach it.
my %INDEX_STATEMENTS = map { $_ => 1 } qw(
    module submodule program end endmodule endsubmodule endinterface endtype
    endsubroutine endfunction endprocedure endprogram
    abstract interface contains use public private
    subroutine function recursive non_recursive pure elemental impure
    integer real complex logical character double type class
);

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

my $root_abs = abs_path($root_dir);
make_path($out_dir) unless -d $out_dir;

# ------------------------------------------------------------
# Data structures
# ------------------------------------------------------------
# %modules:
#   key   = module name (lowercase)
#   value = {
#       name        => original name as seen,
#       files       => { file_path => 1, ... },
#       procedure_count => number of indexed procedures,
#       visibility  => { symbol_name => 'public'|'private' },
#       default_vis => 'public'|'private'|undef,
#       uses        => { module_name => 1 },
#       symbols     => [ { name, kind, file, line, visibility }, ... ]
#   }

my %modules;

# List of all procedures for global index
# @procedures: { name, kind, module, file, line, visibility, procedure_scope }
# procedure_scope is the enclosing program/procedure name (or names joined by ::).

my @procedures;
my @symbols;
my %submodule_parents;
my %procedure_declarations;

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

sub normalize_name {
    my ($name) = @_;
    $name =~ s/^\s+//;
    $name =~ s/\s+$//;
    return lc $name;
}

sub strip_inline_comment {
    my ($line, $quote_state) = @_;
    return $line if !defined $quote_state && index($line, '!') < 0;

    pos($line) = 0;
    if (defined($quote_state) && length($$quote_state)) {
        # Finish a literal carried over from the preceding physical line.
        my $closed = $$quote_state eq "'"
            ? $line =~ /\G[^']*(?:''[^']*)*'(?!')/gc
            : $line =~ /\G[^"]*(?:""[^"]*)*"(?!")/gc;
        return $line unless $closed;
    }

    # Skip text and complete literals in bulk, preserving doubled quotes.
    # The first unconsumed character is a comment marker or an open quote.
    $line =~ /\G [^'"!]*
        (?: (?: '[^']*(?:''[^']*)*'(?!') | "[^"]*(?:""[^"]*)*"(?!") )
            [^'"!]* )*
    /gcx;
    my $end = pos($line);
    my $next = substr($line, $end, 1);
    $$quote_state = $next eq '!' ? '' : $next if defined $quote_state;
    return $next eq '!' ? substr($line, 0, $end) : $line;
}

sub join_continued_statement {
    my ($fh, $line, $physical_line_no) = @_;
    my $quote = '';
    $line = strip_inline_comment($line, \$quote);
    my @parts;
    while (1) {
        # A blank fragment can expose a trailing & in the preceding fragment.
        while (@parts && $line =~ /^\s*$/) {
            $line = pop(@parts) . $line;
        }
        # Inspect only the current fragment; assemble the statement once at the end.
        my $continued = $line =~ s/&\s*$//;
        push @parts, $line;
        last unless $continued;
        my $next_line;
        while (defined($next_line = <$fh>)) {
            ++$$physical_line_no;
            next if $next_line =~ /^\s*(?:!|$)/;
            last;
        }
        last unless defined $next_line;
        # Preserve spaces after an optional leading &, including inside literals.
        $next_line =~ s/^\s*&//;
        $line = strip_inline_comment($next_line, \$quote);
    }
    return join('', @parts);
}

sub csv_quote {
    my ($value) = @_;
    $value //= '';
    $value =~ s/"/""/g;
    return '"' . $value . '"' if $value =~ /[,"\n\r]/;
    return $value;
}

sub relative_path {
    my ($path) = @_;
    my $rel = File::Spec->abs2rel($path, $root_abs);
    $rel =~ s{\\}{/}g;
    return $rel;
}

sub module_record {
    my ($name) = @_;
    my $key = normalize_name($name);
    return $modules{$key} ||= {
        name => $name, files => {}, procedure_count => 0, visibility => {},
        default_vis => undef, uses => {}, symbols => [],
    };
}

sub procedure_prefix {
    my ($prefix, $kind) = @_;
    my ($module, $typed) = (0, 0);
    while ($prefix =~ /\S/) {
        $prefix =~ s/^\s+//;
        if ($prefix =~ s/^(recursive|non_recursive|pure|elemental|impure|module)\b\s*//i) {
            $module = 1 if lc($1) eq 'module';
            next;
        }
        return (0, 0) unless $kind eq 'function' && !$typed &&
            $prefix =~ s/^(double\s+(?:precision|complex)|integer|real|complex|logical|character|type|class)\b\s*//i;
        $typed = 1;
        if ($prefix =~ s/^\(//) {
            # Kind/length selectors can contain nested calls and quoted parentheses.
            my $depth = 1;
            while (length($prefix) && $depth) {
                next if $prefix =~ s/^(?:'(?:[^']|'')*'|"(?:[^"]|"")*")//;
                my $char = substr($prefix, 0, 1, '');
                $depth++ if $char eq '(';
                $depth-- if $char eq ')';
            }
            return (0, 0) if $depth;
        } else {
            # Legacy intrinsic declarations include real*8 and character*(*) forms.
            $prefix =~ s/^\*\s*(?:\d+|\(\s*(?:\d+|\*)\s*\))\s*//;
        }
    }
    return (1, $module);
}

sub add_procedure {
    my ($name, $kind, $owner, $file, $line, $scope, $separate, $interface_depth, $procedure_scope) = @_;
    my $mod = defined $owner ? $modules{$owner} : undef;
    # Module accessibility is resolved after scanning all declarations.
    my $vis = length($procedure_scope) ? 'local' : !$mod ? 'unknown' : 'private';
    my $proc = {
        name => $name, kind => $kind, module => $mod ? $mod->{name} : '',
        file => $file, line => $line, visibility => $vis,
        owner => $owner, scope => $scope, separate => $separate,
        procedure_scope => $procedure_scope,
        declaration => $separate && $interface_depth && !length($procedure_scope),
        module_visibility => $mod && !defined $scope && !length($procedure_scope),
    };
    push @procedures, $proc;
    push @symbols, $proc;
    if ($mod) {
        $mod->{procedure_count}++;
        push @{$mod->{symbols}}, $proc;
        if ($proc->{declaration}) {
            $procedure_declarations{join(':', $owner, $scope // '', normalize_name($name))} = $proc;
        }
    }
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
    my $current_submodule;
    my @procedure_stack;
    my $interface_depth = 0;
    my $in_type = 0;
    my $module_spec = 0;
    my $physical_line_no = 0;
    my $free_form = $file !~ /\.f$/i;

    while (my $line = <$fh>) {
        my $line_no = ++$physical_line_no;
        # Join before filtering, retaining the first physical line for the index.
        if ($free_form && index($line, '&') >= 0) {
            $line = join_continued_statement($fh, $line, \$physical_line_no);
        }
        # Component and binding accessibility belongs to the type, not its module.
        if ($in_type) {
            $in_type = 0 if $line =~ /^\s*end\s*type\b/i;
            next;
        }

        # Skip executable statements, comments and blank lines before detailed parsing.
        my ($keyword) = $line =~ /^\s*([a-zA-Z][a-zA-Z0-9_]*)/;
        next unless defined $keyword && $INDEX_STATEMENTS{lc $keyword};
        $line = strip_inline_comment($line);
        $line =~ s/^\s+//; # Trim indentation once for the declaration checks below.

        # A module declaration ends after its name, unlike module procedure prefixes.
        if ($line =~ /^\s*module\s+([a-zA-Z][a-zA-Z0-9_]*)\s*(?:;|$)/i) {
            my $modname_orig = $1;
            my $modkey = normalize_name($modname_orig);

            my $mod = module_record($modname_orig);
            $mod->{name} = $modname_orig;
            $mod->{files}{$file} = 1;
            $current_submodule = undef;
            $interface_depth = 0;
            @procedure_stack = ();
            $current_module = $modkey;
            $module_spec = 1;
            next;
        }

        # Submodules share their ancestor's procedures, but keep local declarations scoped.
        if ($line =~ /^\s*submodule\s*\(\s*([a-zA-Z][a-zA-Z0-9_]*)\s*(?::\s*([a-zA-Z][a-zA-Z0-9_]*)\s*)?\)\s*([a-zA-Z][a-zA-Z0-9_]*)\s*(?:;|$)/i) {
            my ($ancestor, $parent, $child) = ($1, $2, $3);
            $current_module = normalize_name($ancestor);
            $current_submodule = normalize_name($child);
            $submodule_parents{"$current_module:$current_submodule"} =
                defined $parent ? normalize_name($parent) : '';
            module_record($ancestor)->{files}{$file} = 1;
            $interface_depth = 0;
            @procedure_stack = ();
            $module_spec = 0;
            next;
        }
        if ($line =~ /^\s*end\s*submodule\b/i) {
            $current_module = undef;
            $current_submodule = undef;
            $interface_depth = 0;
            @procedure_stack = ();
            $module_spec = 0;
            next;
        }
        if ($line =~ /^\s*(?:abstract\s+)?interface\b/i) {
            $interface_depth++;
            next;
        }
        if ($line =~ /^\s*end\s*interface\b/i) {
            $interface_depth-- if $interface_depth;
            next;
        }

        # End module
        if ($line =~ /^\s*end\s+module\b/i) {
            $current_module = undef;
            $module_spec = 0;
            @procedure_stack = ();
            next;
        }

        # Interface declarations do not open or close an implementation's scope.
        if (!$interface_depth) {
            if ($line =~ /^\s*end\s*(?:subroutine|function|procedure|program)\b\s*([a-zA-Z][a-zA-Z0-9_]*)?/i) {
                my $name = $1;
                pop @procedure_stack if @procedure_stack &&
                    (!defined $name || lc($name) eq lc($procedure_stack[-1]));
                next;
            }
            if ($line =~ /^\s*end\s*(?:;|$)/i) {
                pop @procedure_stack if @procedure_stack;
                next;
            }
            if ($line =~ /^\s*program\s+([a-zA-Z][a-zA-Z0-9_]*)\b/i) {
                @procedure_stack = ($1);
                next;
            }
        }
        if (!$interface_depth && $line =~ /^\s*contains\s*$/i) {
            $module_spec = 0;
            next;
        }

        # If inside a module, try to pick up PUBLIC/PRIVATE statements
        if (defined $current_module) {
            my $mod = $modules{$current_module};

            if ($line =~ /^\s*use(?:\s*,\s*(?:non_intrinsic|intrinsic)\s*::\s*|\s*::\s*|\s+)([a-zA-Z][a-zA-Z0-9_]*)\b/i) {
                my $used = normalize_name($1);
                $mod->{uses}{$used} = 1 unless $used eq $current_module;
                next;
            }

            # Access statements apply only in the module specification part.
            next if (!$module_spec || $interface_depth) && $line =~ /^\s*(?:public|private)\b/i;

            # A bare access statement sets the module default.
            if ($line =~ /^\s*(public|private)\s*$/i) {
                $mod->{default_vis} = lc($1);
                next;
            }

            # PUBLIC/PRIVATE [::] a, b overrides the default for those names.
            if ($line =~ /^\s*(public|private)(?:\s*::\s*|\s+)(.+)$/i) {
                my ($visibility, $list) = (lc($1), $2);
                my @names = map { normalize_name($_) } split /,/, $list;
                $mod->{visibility}{$_} = $visibility for @names;
                next;
            }
        }

        # A module procedure body inherits its kind from an ancestor interface.
        if (defined $current_submodule && !$interface_depth &&
            $line =~ /^\s*module\s+procedure\s+([a-zA-Z][a-zA-Z0-9_]*)\s*(?:;|$)/i) {
            my $name = $1;
            add_procedure($name, 'procedure', $current_module, $file, $line_no,
                $current_submodule, 1, 0, join('::', @procedure_stack));
            push @procedure_stack, $name;
            next;
        }

        # Find the keyword first to avoid backtracking over lines without a header.
        if ($line =~ /\b(subroutine|function)\s+([a-zA-Z][a-zA-Z0-9_]*)\b/i) {
            my ($prefix, $kind, $name) = (substr($line, 0, $-[0]), lc($1), $2);
            my ($valid, $separate) = procedure_prefix($prefix, $kind);
            if ($valid) {
                add_procedure($name, $kind, $current_module, $file, $line_no,
                    $current_submodule, $separate, $interface_depth, join('::', @procedure_stack));
                push @procedure_stack, $name unless $interface_depth;
                next;
            }
        }

        if (defined $current_module && $line =~ /^\s*type\b(?!\s*(?:\(|is\b|default\b))(?:\s*,\s*([^:]+)::\s*|\s*::\s*|\s+)([a-zA-Z][a-zA-Z0-9_]*)\b/i) {
            my ($attrs, $name) = ($1 // '', $2);
            $in_type = 1;
            my $mod = $modules{$current_module};
            my $procedure_scope = join('::', @procedure_stack);
            my $module_visibility = $module_spec && !$interface_depth && !@procedure_stack;
            if ($module_visibility && $attrs =~ /\b(public|private)\b/i) {
                $mod->{visibility}{normalize_name($name)} = lc($1);
            }
            # As for procedures, module access is finalized after the scan.
            my $vis = length($procedure_scope) ? 'local' : 'private';
            my $symbol = {
                name => $name, kind => 'type', module => $mod->{name},
                file => $file, line => $line_no, visibility => $vis,
                owner => $current_module, scope => $current_submodule,
                procedure_scope => $procedure_scope,
                module_visibility => $module_visibility,
            };
            push @symbols, $symbol;
            push @{$mod->{symbols}}, $symbol;
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
            return unless /\.(?:f|f90|f95)$/i;
            push @files, File::Spec->rel2abs($File::Find::name);
        },
        no_chdir => 1,
    },
    $root_dir
);

foreach my $f (sort @files) {
    process_fortran_file($f);
}

# Resolve after scanning: a submodule may precede its ancestor on disk.
for my $symbol (@symbols) {
    if (defined $symbol->{scope}) {
        $symbol->{module} = $modules{$symbol->{owner}}{name};
    }
    if ($symbol->{module_visibility}) {
        my $mod = $modules{$symbol->{owner}};
        $symbol->{visibility} = $mod->{visibility}{normalize_name($symbol->{name})}
            // $mod->{default_vis} // 'public';
    }
}
for my $proc (@procedures) {
    next unless defined $proc->{scope};
    next if length $proc->{procedure_scope};
    next unless $proc->{separate} && !$proc->{declaration};
    my $scope = $submodule_parents{"$proc->{owner}:$proc->{scope}"} // '';
    my %seen;
    my $declaration;
    while (!$seen{$scope}++) {
        my $key = join(':', $proc->{owner}, $scope, normalize_name($proc->{name}));
        $declaration = $procedure_declarations{$key};
        last if $declaration || $scope eq '';
        $scope = $submodule_parents{"$proc->{owner}:$scope"} // '';
    }
    $proc->{visibility} = $declaration ? $declaration->{visibility} : 'unknown';
    $proc->{kind} = $declaration->{kind} if $declaration;
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
        my @files = sort map { relative_path($_) } keys %{$m->{files}};
        my $files_str = join("<br>", @files);
        print $out "| `$name` | $files_str | $m->{procedure_count} |\n";
    }

    close $out;
}

# ------------------------------------------------------------
# Emit api_index.md
# ------------------------------------------------------------

sub write_api_index_md {
    my ($path) = @_;

    open my $out, '>', $path or die "Cannot write $path: $!";

    print $out "<!-- AUTO-GENERATED by scripts/gen_fortran_indexes.pl; do not manually edit. -->\n\n";
    print $out "# API Index\n\n";
    print $out "| Procedure | Kind | Module | File | Line | Visibility |\n";
    print $out "| --------- | ---- | ------ | ---- | ---- | ---------- |\n";

    # sort by module, then name
    my @sorted = sort {
           ($a->{module} cmp $b->{module})
        || ($a->{name}   cmp $b->{name})
    } @procedures;

    foreach my $p (@sorted) {
        my $pname = $p->{name};
        my $kind  = $p->{kind};
        my $mod   = $p->{module} || '';
        my $file  = relative_path($p->{file});
        my $line  = $p->{line};
        my $vis   = $p->{visibility};

        print $out "| `$pname` | $kind | `$mod` | $file | $line | $vis |\n";
    }

    close $out;
}

sub write_module_index_md {
    my ($path) = @_;
    open my $out, '>', $path or die "Cannot write $path: $!";
    print $out "# Module Index\n\n";

    foreach my $modkey (sort keys %modules) {
        my $m = $modules{$modkey};
        my @files = sort map { relative_path($_) } keys %{$m->{files}};
        print $out "## Module: $m->{name}\n\n";
        print $out "Files:\n";
        print $out "- `$_`\n" for @files;

        my @uses = sort keys %{$m->{uses}};
        if (@uses) {
            print $out "\nUses:\n";
            print $out "- `$_`\n" for @uses;
        }

        for my $visibility (qw(public private)) {
            my @visible = sort { lc($a->{name}) cmp lc($b->{name}) }
                grep { $_->{visibility} eq $visibility } @{$m->{symbols}};
            next unless @visible;
            print $out "\n", ucfirst($visibility), " symbols:\n";
            print $out "- `$_->{name}` — $_->{kind}\n" for @visible;
        }
        my @local = sort {
            lc($a->{procedure_scope}) cmp lc($b->{procedure_scope})
                || lc($a->{name}) cmp lc($b->{name})
        } grep { $_->{visibility} eq 'local' } @{$m->{symbols}};
        if (@local) {
            print $out "\nLocal symbols (enclosing procedure scope):\n";
            print $out "- `$_->{procedure_scope}::$_->{name}` — $_->{kind}\n" for @local;
        }
        print $out "\n---\n";
    }
    close $out;
}

sub write_symbol_index_csv {
    my ($path) = @_;
    open my $out, '>', $path or die "Cannot write $path: $!";
    print $out "module,file,symbol,kind,visibility,line\n";
    my @sorted = sort {
        lc($a->{module} // '') cmp lc($b->{module} // '')
            || lc($a->{name}) cmp lc($b->{name})
            || $a->{line} <=> $b->{line}
    } @symbols;
    for my $s (@sorted) {
        print $out join(',', map { csv_quote($_) }
            $s->{module} // '', relative_path($s->{file}), $s->{name},
            $s->{kind}, $s->{visibility}, $s->{line}), "\n";
    }
    close $out;
}

sub write_module_graph_dot {
    my ($path) = @_;
    open my $out, '>', $path or die "Cannot write $path: $!";
    print $out "digraph module_graph {\n";
    print $out "  \"$modules{$_}{name}\";\n" for sort keys %modules;
    for my $modkey (sort keys %modules) {
        for my $used (sort keys %{$modules{$modkey}{uses}}) {
            print $out "  \"$modules{$modkey}{name}\" -> \"$used\";\n";
        }
    }
    print $out "}\n";
    close $out;
}

# ------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------

my $modules_md   = File::Spec->catfile($out_dir, 'modules.md');
my $api_index_md = File::Spec->catfile($out_dir, 'api_index.md');
my $module_index_md = File::Spec->catfile($out_dir, 'module_index.md');
my $symbol_index_csv = File::Spec->catfile($out_dir, 'symbol_index.csv');
my $module_graph_dot = File::Spec->catfile($out_dir, 'module_graph.dot');

write_modules_md($modules_md);
write_api_index_md($api_index_md);
write_module_index_md($module_index_md);
write_symbol_index_csv($symbol_index_csv);
write_module_graph_dot($module_graph_dot);

print "Wrote:\n  $modules_md\n  $api_index_md\n";
print "  $module_index_md\n  $symbol_index_csv\n  $module_graph_dot\n";
