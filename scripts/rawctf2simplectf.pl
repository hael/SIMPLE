#!/usr/bin/perl
while(<>){
    chomp($_);
    my@vals = split(/\s+/, $_);
    my$dfx    = $vals[1]/10000.;
    my$dfy    = $vals[2]/10000.;
    my$angast = $vals[3];
    print "dfx=$dfx dfy=$dfy angast=$angast\n";
}
