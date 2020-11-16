#!/usr/bin/perl
$ncsyms = $ARGV[0];
@pgrplabels;
$cnt = 0;
foreach $i (0 .. $ncsyms-1){
  $rot = $i + 1;
  $pgrplabels[$cnt] = 'c'.$rot;
  $cnt++;
}
foreach $i (0 .. $ncsyms-1){
  $rot = $i + 1;
  $pgrplabels[$cnt] = 'd'.$rot;
  $cnt++;
}
foreach $i (1 .. 3){
  if( $i == 1 ){
    $pgrplabels[$cnt] = 't';
  }
  if( $i == 2 ){
    $pgrplabels[$cnt] = 'o';
  }
  if( $i == 3 ){
    $pgrplabels[$cnt] = 'i';
  }
  $cnt++;
}
@output = <STDIN>;
print "sym,score,Z\n";
foreach $pgrp (@pgrplabels){
  foreach $line (@output){
    if( $line =~ /\s$pgrp\s/ ){
      @line_arr = split(/\s+/, $line);
      print $line_arr[3].','.$line_arr[5].','.$line_arr[9]."\n";
    }
  }
}






