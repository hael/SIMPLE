#!/usr/bin/perl
while(<>){
  chmod($_);
  open(FOO, "<$_") or die "Cannot open $_: $!\n";
  while($line=<FOO>){
    if( $line =~ /^(.+)Final\sValues$/ ){
      @vals = split(/\s+/,$1);
      $dfx = $vals[1]/10000.;
      $dfy = $vals[2]/10000.;
      $angast = $vals[3];
      print "dfx=$dfx dfy=$dfy angast=$angast\n";
    }
  }
  close(FOO);
}
