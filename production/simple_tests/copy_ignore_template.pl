#!/usr/bin/perl
@folders = glob("simple_*");
chomp(@folders);
foreach $dir (@folders){
  system("cp ignore_template ./$dir/.gitignore");
}
