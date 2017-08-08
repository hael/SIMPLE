#!/usr/bin/perl

use warnings;
use Tie::File;

if( scalar(@ARGV) < 1 ){
  die "Need *star file input on command line\n";
}
my $file = $ARGV[0];
my @reliondoc;
tie @reliondoc, 'Tie::File', $file or die "Cannot tie to file: $file\n";
my @relion_tags;
my @simple_tags;
push @relion_tags, '_rlnVoltage';             push @simple_tags, 'kv';
push @relion_tags, '_rlnSphericalAberration'; push @simple_tags, 'cs';
push @relion_tags, '_rlnAmplitudeContrast';   push @simple_tags, 'fraca';
push @relion_tags, '_rlnDefocusU';            push @simple_tags, 'dfx';
push @relion_tags, '_rlnDefocusV';            push @simple_tags, 'dfy';
push @relion_tags, '_rlnDefocusAngle';        push @simple_tags, 'angast';
push @relion_tags, '_rlnAngleRot';            push @simple_tags, 'e1';
push @relion_tags, '_rlnAngleTilt';           push @simple_tags, 'e2';
push @relion_tags, '_rlnAnglePsi';            push @simple_tags, 'e3';
push @relion_tags, '_rlnOriginX';             push @simple_tags, 'x';
push @relion_tags, '_rlnOriginY';             push @simple_tags, 'y';
my @which_columns;
my @which_tags;
my $linecnt = 0;
my $last_header_line;
foreach my $line (@reliondoc){
    if( $line =~ /^_r/ ){
        my($tag,$col) = look4relion_column($line);
        if( $col != 0 ){
            push @which_columns, $col;
            push @which_tags,    $tag;

        }
        $last_header_line = $linecnt;
    }
    $linecnt++;
    if( $linecnt > 200 ){
        last;
    }
}

foreach my $i ($last_header_line + 1 .. $#reliondoc){
    my $line = $reliondoc[$i];
    my $val;
    chomp ($line);
    $line =~ s/^\s+//;
    my @line_split = split(/\s+/, $line);
    foreach my $i (0 .. $#which_columns){
        if( $simple_tags[$which_tags[$i]] =~ /df/ ){
            $val = $line_split[$which_columns[$i]-1]/10000.;
        }else{
            $val = $line_split[$which_columns[$i]-1];
        }
        print "$simple_tags[$which_tags[$i]]=$val ";
    }
    print "\n";
}

sub look4relion_column{ 
    my $line = shift;
    my $tag  = 0;
    my $col  = 0;
    foreach my $i (0 .. $#relion_tags){
        if( $line =~ /$relion_tags[$i]/ ){
            if( $line =~ /#(\d*)/ ){
                $tag = $i;
                $col = $1;
            }
        }
    }
    return ($tag,$col);
}
