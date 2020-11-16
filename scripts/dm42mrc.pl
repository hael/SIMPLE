#!/usr/bin/perl
while(<>){
    chomp($_);
    $mrc = $_;
    $mrc =~ s/dm4/mrc/;
    system("e2proc2d.py $_ $mrc");
}