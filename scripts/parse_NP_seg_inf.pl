#!/usr/bin/perl
use strict;
use warnings;
my $file = $ARGV[0] or die "Need input csv file on the command line\n";
my $t_per_frame = 0.0025;
my @AoH;
open(my $data, '<', $file) or die "Could not open '$file' $!\n";
while (my $line = <$data>) {
    chomp $line; 
    my @fields = split "," , $line;
    push @AoH, {SEGMENT => $fields[0], FRAME_START => $fields[1], FRAME_END => $fields[2]};
}

for my $i (0 .. $#AoH){
    # parse
    $AoH[$i]{NP_STATS_CSV}                               = $AoH[$i]{SEGMENT}.'/1_atoms_stats/nanoparticle_stats.csv';
    $AoH[$i]{CN_STATS_CSV}                               = $AoH[$i]{SEGMENT}.'/1_atoms_stats/cn_dependent_stats.csv';
    $AoH[$i]{DET_ATMS_LOG}                               = $AoH[$i]{SEGMENT}.'/DET_ATMS';
    ($AoH[$i]{NATOMS},$AoH[$i]{VC_AVG},$AoH[$i]{VC_SIG}) = parse_det_atms($AoH[$i]{DET_ATMS_LOG});
    ($AoH[$i]{DIAM},$AoH[$i]{RADIAL_STRAIN})             = parse_np_stats($AoH[$i]{NP_STATS_CSV});
    $AoH[$i]{PERCEN_CORE}                                = parse_cn_stats($AoH[$i]{CN_STATS_CSV});
    # calculate
    $AoH[$i]{LIFETIME}  = $t_per_frame *  ($AoH[$i]{FRAME_END} - $AoH[$i]{FRAME_START} + 1);
    $AoH[$i]{HALFLIFE}  = $AoH[$i]{LIFETIME}/2;
    $AoH[$i]{TIMESTAMP} = $t_per_frame * (($AoH[$i]{FRAME_END} - $AoH[$i]{FRAME_START} + 1)/2 + $AoH[$i]{FRAME_START});
    if( $i == 0 ){
        $AoH[$i]{ATMS_PER_SEC} = 0;
    }else{
        $AoH[$i]{ATMS_PER_SEC} = ($AoH[$i]{NATOMS} - $AoH[$i-1]{NATOMS}) / ($AoH[$i]{TIMESTAMP} - $AoH[$i-1]{TIMESTAMP});
    }
}

print "SEGMENT,FRAME_START,FRAME_END,LIFETIME,TIMESTAMP,HALFLIFE,NATOMS,DIAM,NATOMS/s,RADIAL_STRAIN,\%CN>8,VC_AVG,VC_SIG\n";
for my $i (0 .. $#AoH){
    print "$AoH[$i]{SEGMENT},$AoH[$i]{FRAME_START},$AoH[$i]{FRAME_END},$AoH[$i]{LIFETIME},$AoH[$i]{TIMESTAMP},";
    print "$AoH[$i]{HALFLIFE},$AoH[$i]{NATOMS},$AoH[$i]{DIAM},$AoH[$i]{ATMS_PER_SEC},$AoH[$i]{RADIAL_STRAIN},";
    print "$AoH[$i]{PERCEN_CORE},$AoH[$i]{VC_AVG},$AoH[$i]{VC_SIG}\n";
}

#for my $i (0 .. $#AoH){
#    print "SEGMENT:       ", $AoH[$i]{SEGMENT},       "\n";
#    print "FRAME_START:   ", $AoH[$i]{FRAME_START},   "\n";
#    print "FRAME_END:     ", $AoH[$i]{FRAME_END},     "\n";
#    print "NP_STATS_CSV:  ", $AoH[$i]{NP_STATS_CSV},  "\n";
#    print "CN_STATS_CSV:  ", $AoH[$i]{CN_STATS_CSV},  "\n";
#    print "DET_ATMS_LOG:  ", $AoH[$i]{DET_ATMS_LOG},  "\n";
#    print "NATOMS:        ", $AoH[$i]{NATOMS},        "\n";
#    print "VC_AVG:        ", $AoH[$i]{VC_AVG},        "\n";
#    print "VC_SIG:        ", $AoH[$i]{VC_SIG},        "\n";
#    print "DIAM:          ", $AoH[$i]{DIAM},          "\n";
#    print "RADIAL_STRAIN: ", $AoH[$i]{RADIAL_STRAIN}, "\n";
#    print "PERCEN_CORE:   ", $AoH[$i]{PERCEN_CORE},   "\n";
#    print "TIMESTAMP:     ", $AoH[$i]{TIMESTAMP},     "\n";
#    print "LIFETIME:      ", $AoH[$i]{LIFETIME},      "\n";
#    print "HALFLIFE:      ", $AoH[$i]{HALFLIFE},      "\n";
#    print "ATMS_PER_SEC   ", $AoH[$i]{ATMS_PER_SEC},  "\n";
#    print "\n";
#}

sub parse_det_atms{
    my $file = shift; 
    my $natoms;
    my $vc_avg;
    my $vc_sig;
    open(my $data, '<', $file) or die "Could not open '$file' $!\n";
    while(my $line = <$data>){
        chomp $line; 
        if( $line =~ /# atoms, final\s+(\d+)/ ){
            $natoms = $1;
        }
        if( $line =~ /VALID_CORR Average:\s+(\d\.\d+)/ ){
            $vc_avg = $1;
        }
        if( $line =~ /VALID_CORR Sigma  :\s+(\d\.\d+)/ ){
            $vc_sig = $1;
        }
    }
    close($file);
    return($natoms,$vc_avg,$vc_sig);
}

sub parse_np_stats{
    my $file = shift;
    my $diam;
    my $radial_strain;
    my $cnt = 0;
    open(my $data, '<', $file) or die "Could not open '$file' $!\n";
    while (my $line = <$data>) {
        $cnt++;
        if( $cnt == 2 ){
            chomp $line;
            my @fields = split "," , $line;
            $diam          = $fields[2];
            $diam          =~ s/^\s+|\s+$//g;
            $radial_strain = $fields[57];
            $radial_strain =~ s/^\s+|\s+$//g;
            close($file);
            return($diam,$radial_strain);
        }
    }
}

sub parse_cn_stats{
    my $file = shift;
    open(my $data, '<', $file) or die "Could not open '$file' $!\n";
    my $cnt = 0;
    my @cn_natms;
    while (my $line = <$data>) {
        $cnt++;
        if( $cnt > 1 ){
            chomp $line; 
            my @fields = split "," , $line;
            $cn_natms[$fields[0]] = $fields[1];
        }
   }
   my $natms_core;
   my $natms_all;
   foreach my $i (0 .. $#cn_natms){
       if( exists $cn_natms[$i] ){
           $natms_all += $cn_natms[$i];
	   if( $i > 8 ){
	       $natms_core += $cn_natms[$i];
	   }
       }
   }
   return 100*($natms_core/$natms_all);
}
