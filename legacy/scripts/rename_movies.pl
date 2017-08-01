#!/usr/bin/perl
while(<>){
  chomp($_);
  if( $_ =~ /(\d+)/ ){
    $new_fname = 'movie'.$1.'.mrc';
    system("mv $_ $new_fname");
  }
}
