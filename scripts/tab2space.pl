#!/usr/bin/perl
use Tie::File;
@code = glob("*.f90");
chomp(@code);
foreach $file (@code){
	tie @lines, 'Tie::File', $file;
	foreach $line (@lines) {
		$line =~ s/\t/    /g;
	}
	untie @lines;
}
