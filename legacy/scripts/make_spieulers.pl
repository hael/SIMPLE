#!/usr/bin/perl
$cnt = 0;
while(<>){
  $cnt++;
  chomp($_);
  if( $_ =~ /e1=([-]*\d+\.\d+)/ ){
    $e1 = $1;
  }
  if( $_ =~ /e2=([-]*\d+\.\d+)/ ){
    $e2 = $1;
  }
  if( $_ =~ /e3=([-]*\d+\.\d+)/ ){
    $e3 = $1;
  }
  print "$cnt\t3\t$e3\t$e2\t$e1\n";
}
