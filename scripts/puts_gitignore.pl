#!/usr/bin/perl
@folders = glob("simple_*");
chomp(@folders);
foreach $folder (@folders){
  print "Trying to chdir to $folder\n";
  system("pwd");
  chdir $folder;
  system("cp ../Temple_gitignore.txt .gitignore");
  chdir("../");
}
