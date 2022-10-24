#!/usr/bin/perl
use warnings;
use strict;
use Tie::File;
use Cwd;
my $dir = getcwd;
chomp($dir);
my $simple_exec     = `readlink -f ../../production/simple_exec/simple_exec.f90`;
my $single_exec     = `readlink -f ../../production/single_exec/single_exec.f90`;
my $git_commit_hash = `git rev-parse --short HEAD`;
chomp($simple_exec);
chomp($single_exec);
chomp($git_commit_hash);
my @simple_exec_doc;
tie @simple_exec_doc, 'Tie::File', $simple_exec or die "Cannot tie to file: $simple_exec\n";
foreach my $line (@simple_exec_doc){
    if( $line =~ /call simple_print_git_version\(\'(.+)\'\)/ ){
        chomp($1);
        $line =~ s/$1/$git_commit_hash/;
    }
}
my @single_exec_doc;
tie @single_exec_doc, 'Tie::File', $single_exec or die "Cannot tie to file: $single_exec\n";
foreach my $line (@single_exec_doc){
    if( $line =~ /call simple_print_git_version\(\'(.+)\'\)/ ){
        chomp($1);
        $line =~ s/$1/$git_commit_hash/;
    }
}
