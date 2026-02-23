#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;

# Files to fix (relative to script directory)
my @files = qw(
  simple_stream.f90
  simple_test_exec.f90
  single_exec.f90
);

my $prod_dir = File::Spec->catdir($Bin, '..', 'production');

# Regex: a single Fortran line calling simple_print_git_version('...').
# Allows leading/trailing whitespace.
# Accepts hex-ish hashes or numeric IDs inside the quotes.
my $CALL_RE = qr/^\s*call\s+simple_print_git_version\s*\(\s*'([0-9A-Fa-f]+)'\s*\)\s*$/i;

sub extract_nonblank_lines {
    my ($text) = @_;
    my @lines = split(/\n/, $text, -1);      # keep empties for structure, filter later
    @lines = grep { $_ !~ /^\s*$/ } @lines;  # drop blank/whitespace-only lines
    return @lines;
}

for my $fname (@files) {
    my $path = File::Spec->catfile($prod_dir, $fname);

    if (!-e $path) { warn "SKIP (not found): $path\n"; next; }
    if (!-f $path) { warn "SKIP (not a regular file): $path\n"; next; }

    open my $in, '<', $path or do { warn "ERROR opening $path for read: $!\n"; next; };
    local $/;
    my $content = <$in>;
    close $in;

    # Fast pre-check: if file can't possibly contain a targeted exact conflict, do nothing.
    if ($content !~ /^<<<<<<< /m || $content !~ /simple_print_git_version\s*\(/i) {
        print "OK (no targeted exact conflicts): $path\n";
        next;
    }

    my @lines = split(/(?<=\n)/, $content);  # keep newlines
    my @out;

    my $changed = 0;
    my $fixed_blocks = 0;

    my $i = 0;
    while ($i < @lines) {
        my $line = $lines[$i];

        # Not a conflict start: pass through
        if ($line !~ /^<<<<<<< /) {
            push @out, $line;
            $i++;
            next;
        }

        # Parse a whole conflict block: <<<<<<< ... ======= ... >>>>>>>
        my @block;
        push @block, $lines[$i++]; # <<<<<<<

        my @first_lines;
        while ($i < @lines && $lines[$i] !~ /^=======\s*$/) {
            push @first_lines, $lines[$i];
            push @block,       $lines[$i];
            $i++;
        }

        # Malformed: missing =======
        if ($i >= @lines) {
            warn "ERROR: Unterminated conflict (missing =======) in $path (file not modified)\n";
            @out = ();
            $changed = 0;
            last;
        }

        push @block, $lines[$i++]; # =======

        my @second_lines;
        while ($i < @lines && $lines[$i] !~ /^>>>>>>> /) {
            push @second_lines, $lines[$i];
            push @block,        $lines[$i];
            $i++;
        }

        # Malformed: missing >>>>>>>
        if ($i >= @lines) {
            warn "ERROR: Unterminated conflict (missing >>>>>>>) in $path (file not modified)\n";
            @out = ();
            $changed = 0;
            last;
        }

        push @block, $lines[$i++]; # >>>>>>>

        # Decide whether THIS conflict is an "exact" simple_print_git_version conflict.
        my $first_text  = join('', @first_lines);
        my $second_text = join('', @second_lines);

        my @first_nonblank  = extract_nonblank_lines($first_text);
        my @second_nonblank = extract_nonblank_lines($second_text);

        my $exact =
            (@first_nonblank == 1 && @second_nonblank == 1
             && $first_nonblank[0]  =~ $CALL_RE
             && $second_nonblank[0] =~ $CALL_RE);

        if ($exact) {
            # Keep ONLY the first call line, preserving its original newline if present.
            # (first_lines may include blanks; we want to output just the call line as it appeared.)
            # Reconstruct by taking the original first side and filtering blanks,
            # but keep the original line content (including indentation) and add a newline.
            my ($kept_line) = grep { $_ !~ /^\s*$/ } @first_lines;

            # Ensure it ends with a newline (since we split with (?<=\n), it should already)
            $kept_line .= "\n" if $kept_line !~ /\n\z/;

            push @out, $kept_line;

            $changed = 1;
            $fixed_blocks++;
        } else {
            # Leave untouched (including conflict markers)
            push @out, @block;
        }
    }

    # If we hit a malformed conflict and cleared @out, skip writing
    if (!@out && $changed == 0 && $content =~ /^<<<<<<< /m) {
        next;
    }

    if (!$changed) {
        print "OK (no exact blocks to resolve): $path\n";
        next;
    }

    my $new_content = join('', @out);

    # Preserve mtime if no net change (extra safety)
    if ($new_content eq $content) {
        print "OK (no changes needed): $path\n";
        next;
    }

    # Safe writeback
    my $tmp = "$path.tmp.$$";
    open my $outf, '>', $tmp or do { warn "ERROR opening $tmp for write: $!\n"; next; };
    print {$outf} $new_content or do { warn "ERROR writing $tmp: $!\n"; close $outf; unlink $tmp; next; };
    close $outf or do { warn "ERROR closing $tmp: $!\n"; unlink $tmp; next; };

    rename $tmp, $path or do {
        warn "ERROR replacing $path with $tmp: $!\n";
        unlink $tmp;
        next;
    };

    print "FIXED ($fixed_blocks exact conflict block(s)): $path\n";
}