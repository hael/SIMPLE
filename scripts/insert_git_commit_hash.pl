#!/usr/bin/perl
use warnings;
use strict;
use Tie::File;
use Cwd;
my $dir = getcwd;
chomp($dir);
my $git_version_file = "$ENV{SIMPLE_PATH}/lib/simple/SimpleGitVersion.h";
my $git_commit_tag   = `git describe --abbrev=8 --always`;
chomp($git_commit_tag);
if( decide_to_substitute($git_version_file) == 1 ){
    substitute($git_version_file);
}

sub substitute{
    my $exec = $_[0];
    my @lines;
    tie @lines, 'Tie::File', $exec or die "Cannot tie to file: $exec\n";
    foreach my $line (@lines){
        if( $line =~ /SIMPLE_GIT_VERSION\s*=\s*"(.+)"/ ){
            chomp($1);
            $line =~ s/\Q$1\E/$git_commit_tag/;
        }
    }
}

sub decide_to_substitute{
    my $exec     = $_[0];
    my @exec_doc = read_file_into_array($exec);
    # decide whether to substitute git commit tag without changing the timestamp
    foreach my $line (@exec_doc){
        if( $line =~ /SIMPLE_GIT_VERSION\s*=\s*"(.+)"/ ){
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
