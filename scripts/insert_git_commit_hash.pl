#!/usr/bin/perl
use warnings;
use strict;
use Tie::File;
use Cwd;
my $dir = getcwd;
chomp($dir);
my $simple_exec    = `readlink -f ../../production/simple_exec/simple_exec.f90`;
my $single_exec    = `readlink -f ../../production/single_exec/single_exec.f90`;
my $stream_exec    = `readlink -f ../../production/simple_stream/simple_stream.f90`;
my $git_commit_tag = `git rev-parse --short HEAD`;
chomp($simple_exec);
chomp($single_exec);
chomp($stream_exec);
chomp($git_commit_tag);

if( decide_to_substitute($simple_exec) == 1 ){
    substitute($simple_exec);
}

if( decide_to_substitute($single_exec) == 1 ){
    substitute($single_exec);
}

if( decide_to_substitute($stream_exec) == 1 ){
    substitute($stream_exec);
}

sub substitute{
    my $exec = $_[0];
    my @lines;
    tie @lines, 'Tie::File', $exec or die "Cannot tie to file: $exec\n";
    foreach my $line (@lines){
        if( $line =~ /call simple_print_git_version\(\'(.+)\'\)/ ){
            chomp($1);
            $line =~ s/$1/$git_commit_tag/;
        }
    }
}

sub decide_to_substitute{
    my $exec     = $_[0];
    my @exec_doc = read_file_into_array($exec);
    # decide whether to substitute git commit tag without changing the timestamp
    foreach my $line (@exec_doc){
        if( $line =~ /call simple_print_git_version\(\'(.+)\'\)/ ){
            chomp($1);
            if( $1 eq $git_commit_tag ){
                return 0;
            } else {
                return 1;
            }
        }
    }
}


sub read_file_into_array {
    my $filename = $_[0];
    # Use a lexical file handle and the three-argument open for safety
    open my $fh, '<', $filename or die "Could not open '$filename': $!";
    # Read the entire file content into an array in list context
    chomp(my @lines = <$fh>); # 'chomp' removes the trailing newline from each line
    close $fh or warn "Error closing file '$filename': $!";
    # Return the array
    return @lines;
}

