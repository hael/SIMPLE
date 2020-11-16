#!/usr/bin/perl
while(<>){
  if( $_ =~ /^distr_simple_script_/ or $_ =~ /^Submitted/ or $_ =~ /^\>\>\> RESOLUTION\:/ or $_ =~ /^\>\>\> DONE/ or $_ =~ /^\*\*\*\*/ or $_ =~ /^DISTRIBUTED/ or $_ =~ /\>\>\> CALCULATING/){
    # don't print anything
  }else{
    print $_;
  }
}
