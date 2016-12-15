#!/usr/bin/perl
$best_res = 50.;
while(<>){
	if( $_ =~ /ctfres\=(.+)/ ){
		if( $1 < $best_res ){
			$best_res = $1;
		}
	}
}
print "Best resolution (observable thon rings): $best_res\n";
