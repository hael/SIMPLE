#!/usr/bin/perl

parse_args($ARGV[0]);

sub parse_args{
    my$prg = shift; # get program name with absolute path
    # make command line instructions
    my$instructions = `/Users/hael/src/fortran/simple_wfrederic/Simple_Restruct.projet/bin/simple_exec prg=$prg`;
    # parse the instructions
    my@splitted = split(/\n/,$instructions);
    my@args;
    foreach (@splitted) {
        if( $_ =~ /^(\w+\d*)\s+\=/ ){
            if( $1 !~ /vol\d+/ ){
                push(@args,$1);
            }
        }
    }
    return @args;
}