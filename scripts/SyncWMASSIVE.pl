#!/usr/bin/perl
use warnings;
use strict;
system("scp -r ./defs/* m2-2:/home/hael/simple/defs/");
system("scp -r ./include/* m2-2:/home/hael/simple/include/");
system("scp -r ./scripts/* m2-2:/home/hael/simple/scripts/");
system("scp -r ./simple_utils/* m2-2:/home/hael/simple/simple_utils/");
system("scp -r ./src/simple_main/* m2-2:/home/hael/simple/src/simple_main/");
chdir "./production/simple";
system("./delete_junk.pl");
system("scp -r * m2-2:/home/hael/simple/production/simple/");
chdir "./production/simple_tests";
system("./delete_junk.pl");
system("scp -r * m2-2:/home/hael/simple/production/simple_tests/");
