#!/usr/bin/perl

@files = glob('*mrcs');
foreach $file (@files){
  $file_new = $file;
  $file_new =~ s/\.mrcs/\_\.mrcs/;  
  rename($file, $file_new);
}

