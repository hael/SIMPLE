#!/usr/bin/perl
while(<>){
    chomp($_);
    if($_ =~ /<pixelSize><x><numericValue>(.+)<\/numericValue>/ ){
        @arr = split('</', $1);
        print $arr[0], "\n";
    }
}