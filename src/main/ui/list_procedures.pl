#!/usr/bin/env perl
use strict;
use warnings;

my $file = shift or die "Usage: $0 <fortran-file>\n";

open(my $fh, "<", $file) or die "Cannot open $file: $!\n";

my %seen;
my $buffer = "";

while (my $line = <$fh>) {

    # Remove trailing newline
    chomp $line;

    # Remove comments (anything after ! unless inside string)
    $line =~ s/!.*$//;

    # Handle continuation lines (&)
    if ($line =~ /&\s*$/) {
        $line =~ s/&\s*$//;
        $buffer .= " $line";
        next;
    } elsif ($buffer ne "") {
        $line = $buffer . " " . $line;
        $buffer = "";
    }

    # Normalize whitespace
    $line =~ s/^\s+|\s+$//g;

    next if $line eq "";

    # Match subroutine
    if ($line =~ /^\s*
        (?:recursive|pure|elemental|module)?\s*
        subroutine\s+
        ([a-zA-Z_]\w*)
        /ix)
    {
        $seen{lc $1} = 1;
        next;
    }

    # Match function
    if ($line =~ /^\s*
        (?:recursive|pure|elemental|module)?\s*
        (?:[a-zA-Z_]\w*\s+)?   # optional type (real, integer, etc.)
        function\s+
        ([a-zA-Z_]\w*)
        /ix)
    {
        $seen{lc $1} = 1;
        next;
    }
}

close $fh;

# Print sorted list
foreach my $name (sort keys %seen) {
    print "$name\n";
}

