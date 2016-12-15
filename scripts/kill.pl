#!/usr/bin/perl
while(<>){
  chomp($_);
  system("kill -9 $_");
}
