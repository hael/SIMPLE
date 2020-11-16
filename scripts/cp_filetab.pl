#!/usr/bin/perl
while(<>){
  chomp($_);
  system("cp $_ .");
}
