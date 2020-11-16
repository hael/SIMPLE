#!/usr/bin/perl
$cntx = 0;
$cnty = 0;
$cnta = 0;
@dfx;
@dfy;
@angast;
while(<>){
  chomp($_);
  if( $_ =~ /<defocusU unit="nm">(\d+\.\d+)<\/defocusU>/ ){
    $dfx[$cntx] = $1;
    $cntx++;
  }
  if( $_ =~ /<defocusV unit="nm">(\d+\.\d+)<\/defocusV>/ ){
    $dfy[$cnty] = $1;
    $cnty++;
  }
  if( $_ =~ /<defocusUAngle unit="deg">([-]*\d+\.\d+)<\/defocusUAngle>/ ){
    $angast[$cnta] = $1;
    $cnta++;
  }
	
}
if( $cntx == $cnty ){
   if( $cntx == $cnta ){
      foreach $i (0 .. $#dfx) {
         $df_x = $dfx[$i]/1000.;
         $df_y = $dfy[$i]/1000.;
         print "dfx=$df_x dfy=$df_y angast=$angast[$i]\n";
      }
   }else{
      die "Nonconforming number of dfx/angasts";
   }
}else{
   die "Nonconforming number of dfx/dfys";
}
