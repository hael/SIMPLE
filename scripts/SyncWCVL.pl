#!/usr/bin/perl
use warnings;
use strict;
system("scp -r ./defs/* cvl:/home/hael/simple/defs/");
system("scp -r ./include/* cvl:/home/hael/simple/include/");
system("scp -r ./scripts/* cvl:/home/hael/simple/scripts/");
system("scp -r ./simple_utils/* cvl:/home/hael/simple/simple_utils/");
system("scp -r ./src/simple/* cvl:/home/hael/simple/src/simple/");
chdir "./production/simple";
system("./delete_junk.pl");
system("scp -r * cvl:/home/hael/simple/production/simple/");
