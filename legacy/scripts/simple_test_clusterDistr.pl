#!/usr/bin/perl
use simple_user_input;
use strict;
use warnings;
our (@ISA, @EXPORT);
use Exporter;
@ISA = ("Exporter");
use simple_clusterDistr;
print generate_distr_script(10, 10, 'here', 'hejhopp', 1, 5, 1);
print "###############\n";
print generate_shmem_distr_script(10, 10, 'here', 'hejhopp');