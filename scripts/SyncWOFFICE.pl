#!/usr/bin/perl
use warnings;
use strict;
my $path ='/Users/hael/src/fortran/lib/simple_wfrederic/Simple_Restruct.projet';
system("scp -r ./defs/* office:$path/defs/");
system("scp -r ./include/* office:$path/include/");
system("scp -r ./scripts/* office:$path/scripts/");
system("scp -r ./simple_utils/* office:$path/simple_utils/");
system("scp -r ./src/simple/* office:$path/src/simple/");
chdir "./production/simple";
system("./delete_junk.pl");
system("scp -r * office:$path/production/simple/");
