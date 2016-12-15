#!/usr/bin/perl
while(<>){
  chmod($_);
  open(FOO, "<$_") or die "Cannot open $_: $!\n";
  while($line=<FOO>){
    if( $line =~ /Estimated defocus values(.+)Angstroms/ ){
        @dfxy = split(',', $1);
        $dfxy[0] =~ s/://;
	$dfxy[0] =~ s/\s+//g;
	$dfxy[1] =~ s/\s+//g;
	$dfxy[0] = $dfxy[0]/10000.; # Anstroms to microns
	$dfxy[1] = $dfxy[1]/10000.; # Anstroms to microns 
    }
    if( $line =~ /Estimated azimuth of astigmatism(.+)degrees/){
        $angast = $1;
	$angast =~ s/://;
	$angast =~ s/\s+//g;
	$angast_arr[$cnt] = $angast;
	print "dfx=$dfxy[0] dfy=$dfxy[1] angast=$angast\n";
    }
     
  }
  close(FOO);
}
	
