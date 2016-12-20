#!/usr/bin/perl
use warnings;
use strict;
my %all_args;
#-------------------------------------------
# Extract the program names
#-------------------------------------------
my @simple_files = glob("simple_*");
chomp(@simple_files);
my $i = 0;
while( $i <= $#simple_files ){
    if( $simple_files[$i] =~ /\.f90$/ or ($simple_files[$i] =~ /\.pl$/ or ($simple_files[$i] =~ /\.txt$/ or ($simple_files[$i] =~ /\.log$/ or ($simple_files[$i] =~ /simple_test/ or ($simple_files[$i] =~ /simple_unit/ or ($simple_files[$i] =~ /simple_h5/ or $simple_files[$i] =~ /simple_lib/))))))){      
        splice @simple_files, $i, 1;
    } else {
        $i++;
    }
}

my $duplicates = 0;
foreach my$file (@simple_files){
    my %args = parse_args($ENV{SIMPLEBIN}.$file);
    foreach my $key(keys %args){
        if( defined($all_args{$key}) ){
            $duplicates++;
        }
        $all_args{$key} = $args{$key};
    }
}
#print "Nr of duplicates: $duplicates\n";
foreach my $key(sort {lc $a cmp lc $b} keys %all_args) {
    print "$key $all_args{$key}\n";
}

sub parse_args{
    my$prg = shift; # get program name with absolute path
    # make command line instructions
    my$instructions = `$prg`;
    # parse the instructions
    my@splitted = split(/\s/,$instructions);
    my%args;
    foreach (@splitted) {
        if( $_ =~ /(\w+\d*)\=((.+))/){
            if( $1 !~ /vol\d+/ ){
                $args{$1} = $2;
            }  
        }
    }
    return %args;
}