#!/usr/bin/perl
@folders = glob("simple_*");
chomp(@folders);
foreach my $folder (@folders){
	chdir($folder);
	print ">>> COMPILING & LINKING: $folder\n";
	system("./compile_and_link.csh");
	chdir("../");
}
