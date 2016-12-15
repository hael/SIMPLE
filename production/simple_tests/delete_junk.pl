#!/usr/bin/perl
@folders = glob("simple_*");
chomp(@folders);
foreach my $folder (@folders){
	chdir($folder);
	print ">>> DELETING JUNK FROM: $folder\n";
	system("rm -rf *.dSYM SIMPLE_UNIT_TEST*");
	chdir("../");
}
