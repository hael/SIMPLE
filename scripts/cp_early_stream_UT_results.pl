#!/usr/bin/perl

use strict;
use warnings;
use POSIX qw(strftime);

my$cavgs_name           = 'cavgs_iter024_ranked.mrc';
my$classdoc_ranked_name = 'classdoc_ranked.txt';
my$target_folder        = '/home/elmlundho/data/early_stream_UT/';
my$test_folder          = '/home/elmlundho/data/picker_testing/';

# Create date-time string (example format: 20260216_143522)
my $datetime = strftime("%Y%m%d_%H%M%S", localtime);

# Define read-only list of folders
my @folders = (
    'Apoferritin/7_check_refpick/',
    'CLC7/6_check_refpick/',
    'FlipQR/9_check_refpick/',
    'HCN/6_check_refpick/',
    'Membralin/3_check_refpick/',
    'MotAB/5_check_refpick/',
    'Not/4_check_refpick/',
    'PepT2/3_check_refpick/',
    'RNApol/2_check_refpick/',
    'RYPER/3_check_refpick/',
    'SLC/3_check_refpick/',
    'TatBC/2_check_refpick/',
    'TRPM4/3_check_refpick/',
);

my $target_subfolder;
my $copy_from;
foreach my $folder (@folders) {
    $target_subfolder = $folder;
    $target_subfolder =~ s{^([^/]+)/.+}{$1/$datetime/};
    $target_subfolder = $target_folder.$target_subfolder; 
    system("mkdir -p $target_subfolder");
    $copy_from = $test_folder.$folder.$cavgs_name;
    system("cp $copy_from $target_subfolder");
    $copy_from = $test_folder.$folder.$classdoc_ranked_name;
    system("cp $copy_from $target_subfolder");
}

