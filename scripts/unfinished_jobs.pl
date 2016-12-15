#!/usr/bin/perl
@jobs = glob("JOB_FINISHED_*");
chomp(@jobs);
if( scalar(@ARGV) != 1 ){
  die "Give nr of partitions to script!\n"; 
}
@indicator;
foreach $i (1 .. $ARGV[0]){
  $indicator[$i] = 0;
}
foreach $job (@jobs){
  if( $job =~ /(\d+)/ ){
    $indicator[$1] = 1;
  }
}
foreach $i (1 .. $ARGV[0]){
  if( $indicator[$i] == 0 ){
    print "This job is unfinished: $i\n";
  }
}
