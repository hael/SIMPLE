#!/usr/bin/perl
use lib '/Users/hael/src/fortran/simple3.0/scripts';
use lib '/Users/hael/src/fortran/simple3.0';
use strict;
use warnings;
our (@ISA, @EXPORT);
use Exporter;
@ISA = qw(Exporter);
use simple_user_input;
use simple_clusterSpecs;
use simple_clusterDistr;
use Cwd qw(getcwd);
use List::Util qw(max min);
# check so that at least one argument
if(scalar(@ARGV) == 0){die "need at least one argument: prg=<simple program name>\n"};
################################################################################
# Declare global variables                                                     #
################################################################################
my%name_value;
my$folder;
my$primedoc;
my$lpstart;
my$lpstop;
my$dynlp;
my$command;
my$exec_mode;
my$ext;
my$round_cnt;
my$frac_srch_space=0.0;
my$SIMPLEBIN=$SIMPLE_PATH.'/bin';
my@supported_prgs;
################################################################################
# set programs supported by distr script                                       #
################################################################################
push(@supported_prgs, 'comlin_smat');
push(@supported_prgs, 'eo_recvol');
push(@supported_prgs, 'integrate_movies');
push(@supported_prgs, 'prime2D');
push(@supported_prgs, 'prime3D');
push(@supported_prgs, 'recvol');
push(@supported_prgs, 'simimgs');
push(@supported_prgs, 'unblur_movies');
push(@supported_prgs, 'stackops');
push(@supported_prgs, 'multiptcl_init');
################################################################################
# set dependent programs                                                       #
################################################################################
$name_value{'merge_algndocs'}     = $SIMPLEBIN.'/simple_exec prg=merge_algndocs';
$name_value{'merge_similarities'} = $SIMPLEBIN.'/simple_exec prg=merge_similarities';
$name_value{'volassemble'}        = $SIMPLEBIN.'/simple_exec prg=volassemble';
$name_value{'eo_volassemble'}     = $SIMPLEBIN.'/simple_exec prg=eo_volassemble';
$name_value{'cavgassemble'}       = $SIMPLEBIN.'/simple_exec prg=cavgassemble';
$name_value{'check3D_conv'}       = $SIMPLEBIN.'/simple_exec prg=check3D_conv';
$name_value{'check2D_conv'}       = $SIMPLEBIN.'/simple_exec prg=check2D_conv';
$name_value{'check_nptcls'}       = $SIMPLEBIN.'/simple_exec prg=check_nptcls';
$name_value{'check_box'}          = $SIMPLEBIN.'/simple_exec prg=check_box';
$name_value{'resrange'}           = $SIMPLEBIN.'/simple_exec prg=resrange';
$name_value{'split_pairs'}        = $SIMPLEBIN.'/simple_exec prg=split_pairs';
################################################################################
# set other dependencies                                                       #
################################################################################
if( $name_value{'prg'} =~ /unblur/ ){
    $name_value{'fbody'} = 'intgmovie';
}else{
    $name_value{'fbody'} = $ALGNDOC_FBODY;
}
$name_value{'execdir'} = getcwd();
################################################################################
# Start of the execution commands                                              #
################################################################################
# parse the command line
# %name_value is the hash with input data
foreach my $cmd (@ARGV){
    chomp($cmd); # remove return characters
    my @tmp = split('=', $cmd);
    $tmp[0] =~ s/\s//g; # remove whitespaces
    $name_value{$tmp[0]} = $tmp[1];
}
# check so that prg is set (program name)
if(defined $name_value{'prg'}){
}else{
    die "prg=<simple program name> need to be set\n";
}
# check so that program is supported by distr_script
my$supported = 0;
foreach my$i (0 .. $#supported_prgs){
    if( $name_value{'prg'} =~ /$supported_prgs[$i]/ ){
        $supported = 1;
    }
}
if( $supported == 0 ){
    die "The $name_value{'prg'} program is not supported by distr_simple.pl\n";
}
# check so that nparts is set (number of partitions)
if( !defined $name_value{'nparts'} ){
    die "nparts need to be set\n";
}
# check wether we are executing prime or not
$exec_mode = 0;
if( $name_value{'prg'} =~ /prime3D/ ){
    $exec_mode = 1;
}
if( $name_value{'prg'} =~ /oasis/ ){
    $exec_mode = 2;
}
if( $name_value{'prg'} =~ /prime2D/ ){
    $exec_mode = 3;
}
if( $name_value{'prg'} =~ /smat/ ){
    $exec_mode = -1;
}
# make absolute path to the exec program
my $abs_exec = $SIMPLEBIN.'/simple_exec';
delete $name_value{'exec'};
$name_value{'exec'} = $abs_exec;
# make command line instructions
my$instructions = `$name_value{'exec'} prg=$name_value{'prg'}`;
$instructions =~ s/SIMPLE_/DISTR_SIMPLE_/;
$instructions =~ s/\[nthr\=\<nr of OpenMP threads\{1\}\>\]/nparts\=\<nr of partitions\>/;
$instructions = $instructions."\n";
# if too few argument given, print instructions
if(scalar(@ARGV) < 3){
    print $instructions;
}
# check so that a stack is inputted
if( !defined($name_value{'filetab'}) ){
    if( !defined($name_value{'stk'}) ){
        if( $name_value{'prg'} !~ /simimgs/ and  $name_value{'prg'} !~ /volume_smat/ ){
            die "Need a stack of particle images 4 distr_simple exec!\n";
        }
    }
}
# copy vollist (if defined) to filetab key
if( defined($name_value{'vollist'}) ){
    if( defined($name_value{'filetab'}) ){
        delete $name_value{'filetab'};
    }
    $name_value{'filetab'} = $name_value{'vollist'};
}
my @pieces;
if( defined($name_value{'filetab'}) ){
    open(FHANDLE, "<$name_value{'filetab'}") or die "Cannot open $name_value{'filetab'} for reading: $!\n";
    my$first_movie_file;
    my$nlines = 0;
    while(<FHANDLE>){
        chomp($_);
        $nlines++;
        if( $nlines == 1 ){
            $first_movie_file = $_;
        }
    }
    # set the number of particles = $nlines (number of movies)
    if( $exec_mode < 0 ){
    }else{
        $name_value{'nptcls'} = $nlines;
    }
    # extract file-suffix
    @pieces = split(/\./,$first_movie_file);
    $ext = '.'.$pieces[-1];
    # remove the vollist from the filetab key
    if( defined($name_value{'vollist'}) ){
        delete $name_value{'filetab'};
    }
}else{
    # extract file-suffix
    if( defined($name_value{'stk'}) ){
        @pieces = split(/\./,$name_value{'stk'});
    }elsif( defined($name_value{'vol1'}) ){
        @pieces = split(/\./,$name_value{'vol1'});
    }else{
        die "Need stack (stk) or volume (vol1) to be defined on the command line to figure out the file extension\n";
    }
    $ext = '.'.$pieces[-1];
    if( $ext eq '.mrcs'){
        $ext = '.mrc;'
    }
    print "file extension: $ext\n";
}
if( !defined($name_value{'filetab'}) ){
    # check so that box is set
    if( !defined($name_value{'box'}) ){
        $name_value{'box'} = exec_check_box();
    }
    # set the number of particles
    if( $name_value{'prg'} =~ /simimgs/ ){
        if( !defined($name_value{'nptcls'}) ){
            die "Need nptcls to be set at the command line\n";
        }
    }else{
        if( defined($name_value{'nptcls'}) ){
            delete $name_value{'nptcls'};
        }
        my $nptcls_tmp = exec_check_nptcls();
        if( $exec_mode == -1 ){
            $name_value{'nptcls'} = ($nptcls_tmp*($nptcls_tmp-1))/2;
        }else{
            $name_value{'nptcls'} = $nptcls_tmp;
        }
    }
    if( $exec_mode == 3 ){
        if( !defined($name_value{'ncls'}) ){
            # set the number of particles
            $name_value{'ncls'} = exec_check_ncls();
        }
    }
}
# check so that nthr is NOT set (number of threads)
if( defined $name_value{'nthr'} ){
    print "Setting nthr (number of threads) on the command line is no longer allowed!\n";
    print "Please modify the /simple/scripts/simple_clusterDistr.pm module to configure your distributed execution.\n";
    print "The cluster specifications are defined in /simple/scripts/simple_clusterSpecs.pm\n";
    die "If you need assistance with creating a specification for a new environment, please contact: hans.elmlund\@monash.edu\n";
}
# check so that ndocs is set (number of documents=number of partitions)
delete $name_value{'ndocs'};
$name_value{'ndocs'} = $name_value{'nparts'};
# check so that pgrp is set
if( !defined($name_value{'pgrp'}) ){
    $name_value{'pgrp'} = 'c1';
}
# check so that outvol is not defined
if( defined $name_value{'outvol'} ){
    die "Param: outvol cannot be defined! Used internally...";
}
if( !defined($name_value{'maxits'}) ){
    $name_value{'maxits'} = 500;
}
# check so that refine is set
if( !defined($name_value{'refine'})){
    if( $exec_mode == 1 or $exec_mode == 3 ){
        $name_value{'refine'} = 'no';
    }elsif($exec_mode == 2){
        $name_value{'refine'} = 'soft';
    }
}
if( $exec_mode > 0 ){
    if( $exec_mode < 2 ){
        # check so that a volume is inputted
        if( !defined($name_value{'vol1'}) ){
            die "Need at least one starting volume 4 exec!\n";
        }
    }elsif( $name_value{'prg'} =~ /cluster/ ){
        # check so that references are inputted
        if( !defined($name_value{'refs'}) ){
            die "Need starting references 4 prime2D exec!\n";
        }
    }
    # check so that fstep is set
    if( !defined $name_value{'fstep'} ){
        $name_value{'fstep'} = 1;
    }
    # check so that eo is set
    if( !defined($name_value{'eo'}) ){
        if( $exec_mode == 1 or $exec_mode == 3 ){
            $name_value{'eo'} = 'no';
        }else{
            $name_value{'eo'} = 'yes';
        }
    }
    # set dynlp
    if( $exec_mode == 1 ){
        $dynlp = 'yes';
    }else{
        $dynlp = 'no';
    }
    if( defined($name_value{'dynlp'}) ){
        $dynlp = $name_value{'dynlp'};
    }elsif($name_value{'eo'} eq 'yes' or defined($name_value{'lp'}) ){
        $dynlp = 'no';
    }
    # check so that nspace is set (number of projection directions)
    if( !defined $name_value{'nspace'} ){
        $name_value{'nspace'} = 1000;
    }
    if( $exec_mode == 1 and  $dynlp eq 'yes' ){
        # determine initial and final Fourier index
        if( $dynlp eq 'yes' ){
            ($lpstart, $lpstop) = exec_resrange();
            if( defined($name_value{'lp'}) ){
                if( defined($name_value{'lpstart'}) ){
                    $name_value{'lp'} = $name_value{'lpstart'};
                }else{
                    $name_value{'lp'} = $lpstart;
                }
            }
            if( not defined($name_value{'lpstop'}) ){
                $name_value{'lpstop'} = $lpstop;
            }
            # determine the initial low-pass limit
            delete $name_value{'find'};
            $name_value{'find'} = get_find($name_value{'lp'});
            # determine the final low-pass limit
            delete $name_value{'fmax'};
            if( defined($name_value{'lpstop'}) ){
                $name_value{'fmax'} = get_find($name_value{'lpstop'});
            }else{
                $name_value{'fmax'} = get_find($lpstop);
            }
        }
    }
    # set number of states (by counting input volumes)
    my @keys    = keys %name_value;
    my $nstates = 0;
    foreach my $key(@keys){
        if($key =~ /vol\d/ ){$nstates++};
    }
    $name_value{'nstates'} = $nstates;
}else{
    if( $name_value{'prg'} =~ /recvol/ ){
        if( defined($name_value{'eo'}) ){
        }else{
            if( $name_value{'prg'} =~ /eo/ ){
                $name_value{'eo'} = 'yes';
            }else{
                $name_value{'eo'} = 'no';
            }
        }
    }
}
# do the refinement
my$converged = 0;
my$round;
my$round_stop;
my$lp;
my$refs;
my$round_str;
if( defined($name_value{'startit'}) ){
    $round = $name_value{'startit'}-1;
} else {
    $round = 0;
}
my$updated_res = 0;
$round_cnt = 0;
if( $exec_mode > 0 ){
    while( $round < $name_value{'maxits'} ){
        # make refinement folder
        $round++;
        $round_cnt++;
        $round_str = zero_pad_intg($round,500);
        if( $exec_mode == 1 ){
            $folder   = 'prime3Dround_'.$round_str;
            $primedoc = 'prime3Ddoc_'.$round_str.'.txt';
        }elsif( $exec_mode == 2 ){
            $folder   = 'oasis3Dround_'.$round_str;
            $primedoc = 'oasis3Ddoc_'.$round_str.'.txt';
        }elsif( $exec_mode == 3 ){
            $folder   = 'prime2Dround_'.$round_str;
            $primedoc = 'prime2Ddoc_'.$round_str.'.txt';
            $refs     = 'cavgs_iter'.$round_str.$ext;
            delete $name_value{'which_iter'};
            $name_value{'which_iter'} = $round;
        }
        mkdir $folder;
        # update shift range
        if( $frac_srch_space >= 90.0 ){
            if( not defined($name_value{'trs'}) ){
                $name_value{'trs'} =    0.025*$name_value{'box'};
                $name_value{'trs'} = max 2.0, $name_value{'trs'};
                $name_value{'trs'} = min 6.0, $name_value{'trs'};
            } 
        }
        # parallell execution
        system("rm -f JOB_FINISHED_*");
        exec_para($name_value{'prg'});
        # determine when the jobs have finished
        sleep(10) while( njobs_finished() < $name_value{'nparts'} );
        system("rm -f JOB_FINISHED_* OUT* distr_script*");
        # merge alignment docs
        exec_merge_algndocs($primedoc);
        if( $exec_mode <= 2 ){
            # reconstruct/assemble volumes
            if( $name_value{'eo'} eq 'no' ){
                exec_assemble_volumes();
            }else{
                $lp = exec_assemble_eo_volumes();
            }
        }elsif( $exec_mode == 3 ){
            while( defined($name_value{'refs'}) ){
                delete $name_value{'refs'};
            }
            $name_value{'refs'} = $refs;
            exec_assemble_cavgs();
        }else{
            die "unknown exec_mode mode; distr_simple.pl\n";
        }
        # check convergence
        my $conv;
        if( $exec_mode == 3 ){
            $conv = exec_check2D_conv();
        }else{
            $conv = exec_check3D_conv();
        }
        if( $round > 1 ){
            if( $conv == 2 ){
                cleanup();
                die "**** DISTR_SIMPLE_PRIME NORMAL STOP ****\n";
            }elsif ( $conv == 1 ){
                if( !defined($name_value{'find'}) ){
                    die "Want to do dynamic low-pass update, but find is not defined, weird!\n";
                }
                $name_value{'find'} = $name_value{'find'}+$name_value{'fstep'};
                $updated_res = 1;
                if( $name_value{'find'} > $name_value{'fmax'} ){
                    $name_value{'find'} = $name_value{'fmax'};
                }
            }
        }
    }
}else{ # not PRIME execution
    if( $exec_mode < 0 ){
        # we are operating on pairs, so we need to create mapping
        exec_split_pairs();
    }
    exec_para($name_value{'prg'});
    if( $name_value{'prg'} =~ /recvol/ ){
        # determine when the jobs have finished
        sleep(10) while( njobs_finished() < $name_value{'nparts'} );
        system("rm -f JOB_FINISHED_* OUT* distr_script*");
        # reconstruct/assemble volumes
        if( $name_value{'eo'} eq 'no' ){
            exec_assemble_volumes();
        }else{
            $lp = exec_assemble_eo_volumes();
        }
    }elsif( $exec_mode < 0 ){
        # determine when the jobs have finished
        sleep(10) while( njobs_finished() < $name_value{'nparts'} );
        system("rm -f JOB_FINISHED_* OUT* distr_script*");
        # we generated pairwise similarities, so we need to merge them
        exec_merge_similarities();
        system("rm -f pairs_part* similarities_part*");
    }elsif( $name_value{'prg'} =~ /unblur/ ){
        # determine when the jobs have finished
        sleep(10) while( njobs_finished() < $name_value{'nparts'} );
        system("rm -f JOB_FINISHED_* OUT* distr_script*");
        exec_merge_docs('unblur_movies_params_part', 'unblur_movies_params_merged.txt');
    }
}
################################################################################
# Subroutines                                                                  #
################################################################################
#
sub exec_para{
    my$prog  = shift;
    my$i;
    my$vol;
    my@submitted;
    my$subout;
    my@args = parse_args($prog);
    my$cmd_string = $SIMPLEBIN.'/simple_exec prg='.$prog.' ';
    $i = 1;
    $vol = 'vol'.$i;
    while( defined($name_value{$vol}) ){
        $cmd_string = add2string($vol, $cmd_string);
        $i++;
        $vol = 'vol'.$i;
    }
    foreach (@args) {
        $cmd_string = add2string($_, $cmd_string);
    }
    # these need to be manually added in (not necessarily part of command line)
    $cmd_string = add2string('box',   $cmd_string);
    $cmd_string = add2string('find',  $cmd_string);
    $cmd_string = add2string('nparts', $cmd_string);
    $cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $cmd_string, "\n";
    my $ptcls_per_part = int(int($name_value{'nptcls'})/int($name_value{'nparts'}));
    my$leftover = int($name_value{'nptcls'})-$ptcls_per_part*int($name_value{'nparts'});
    my $stop  = 0;
    my $start = 0;
    for(my $i=1; $i<=$name_value{'nparts'}; $i++){
        my$i_str = zero_pad_intg($i,$name_value{'nparts'});
        if( $i == $name_value{'nparts'} ){
            $start = $stop+1;
            $stop  = $name_value{'nptcls'};
        }else{
            if( $leftover == 0 ){
                $start = $stop+1;
                $stop  = $start+$ptcls_per_part-1;
            }else{
                $stop  = $i*($ptcls_per_part+1);
                $start = $stop-($ptcls_per_part+1)+1;
                $leftover--;
            }
        }
        my $time  = $TIME_PER_IMAGE*($stop-$start+1);
        my $hours = min 23,int($time/3600);
        my $min   = (($time/60)%60);
        open(FHANDLE, ">distr_script_$i_str") or die "Cannot open distr_script_$i_str for writing: $!\n";
        print FHANDLE generate_distr_script($hours, $min, $name_value{'execdir'}, $cmd_string, $start, $stop, $i_str);
        close(FHANDLE);
        chmod 0777, 'distr_script_'.$i_str;
        if( $SIMPLESYS eq 'LOCAL' ){
            system("nohup ./distr_script_$i_str &\n");
        }else{
            my$subcmd = "$SUBMITCMD ./distr_script_$i_str";
            $subout = `$subcmd`;
            # nothing else 4 now
        }
    }
}

sub njobs_finished{
    my$nfini = 0;
    my@jobs;
    @jobs = <JOB_FINISHED_*>;
    return scalar(@jobs);
}

sub exec_assemble_volumes{
    my$i;
    my$vol;
    my@refvols;
    my@recvols;
    my@args = parse_args('volassemble');
    my$assemble_cmd_string = $name_value{'volassemble'}.' ';
    foreach (@args) {
        $assemble_cmd_string = add2string($_, $assemble_cmd_string);
    }
    $assemble_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $assemble_cmd_string, "\n";
    # check number of kernels and recbins
    my@kernels = <rho*part*>;
    chomp(@kernels);
    my@recbins = <recvol_state*part*>;
    chomp(@recbins);
    # execute assemble command
    system($assemble_cmd_string);
    # delete intermediate files
    foreach my$i (0 .. $#kernels){
        system("rm -f $kernels[$i] $recbins[$i]");
    }
    # take care of the reconstructed volumes
    @recvols = <recvol_state*$ext>;
    chomp(@recvols);
    @recvols = grep {$_ =~ /^recvol_state\d+$ext$/} @recvols;
    my $nrecvols = scalar(@recvols);
    @refvols = @recvols;
    foreach $i ( 1 .. $nrecvols ){
        $vol = 'vol'.$i;
        if( defined($name_value{$vol}) ){
            delete $name_value{$vol};
        }
        $name_value{$vol} = $refvols[$i-1];
        if( defined($folder) ){
            if( -e $folder ){
                system("cp $name_value{$vol} ./$folder");
            }
        }
    }
}

sub exec_assemble_eo_volumes{
    my$i;
    my$vol;
    my@refvols;
    my@recvols;
    my$lp;
    my@args = parse_args('eo_volassemble');
    my$assemble_cmd_string = $name_value{'eo_volassemble'}.' ';
    foreach (@args) {
        $assemble_cmd_string = add2string($_, $assemble_cmd_string);
    }
    $assemble_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $assemble_cmd_string, "\n";
    # check number of kernels and recbins
    my@even_kernels = <rho*part*_even*>;
    chomp(@even_kernels);
    my@even_recbins = <recvol_state*part*_even*>;
    chomp(@even_recbins);
    my@odd_kernels = <rho*part*_odd*>;
    chomp(@odd_kernels);
    my@odd_recbins = <recvol_state*part*_odd*>;
    chomp(@odd_recbins);
    my$assemble_eo_info = `$assemble_cmd_string`;
    my@assemble_eo_info_lines = split(/\n/,$assemble_eo_info);
    foreach my$line (@assemble_eo_info_lines){
        print $line, "\n";
        if( $line =~ /\>\>\> LOW-PASS LIMIT\:\s+(\d+.\d+)/ ) { $lp = $1 };
    }
    # delete intermediate files
    foreach my$i (0 .. $#even_kernels){
        system("rm -f $even_kernels[$i] $even_recbins[$i]");
        system("rm -f $odd_kernels[$i] $odd_recbins[$i]");
        system("rm -f rho_recvol_state*_even*");
        system("rm -f recvol_state*_even*");
        system("rm -f rho_recvol_state*_odd*");
        system("rm -f recvol_state*_odd*");
    }
    # take care of the reconstructed volumes
    @recvols = <recvol_state*$ext>;
    chomp(@recvols);
    @recvols = grep {$_ !~ /^recvol_state\d+_even.+$/} @recvols;
    @recvols = grep {$_ !~ /^recvol_state\d+_odd.+$/} @recvols;
    @recvols = grep {$_ =~ /^recvol_state\d+$ext$/} @recvols;
    my $nrecvols = scalar(@recvols);
    @refvols = @recvols;
    foreach $i ( 1 .. $nrecvols ){
        $vol = 'vol'.$i;
        if( defined($name_value{$vol}) ){
            delete $name_value{$vol};
        }
        $name_value{$vol} = $refvols[$i-1];
        if( defined($folder) ){
            if( -e $folder ){
                system("cp $name_value{$vol} ./$folder");
                system("cp fsc_state*bin ./$folder");
            }
        }
    }
    return $lp;
}

sub exec_assemble_cavgs{
    my@args = parse_args('cavgassemble');
    my$assemble_cmd_string = $name_value{'cavgassemble'}.' ';
    foreach (@args) {
        $assemble_cmd_string = add2string($_, $assemble_cmd_string);
    }
    $assemble_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $assemble_cmd_string, "\n";
    # execute assemble command
    system($assemble_cmd_string);
    # delete intermediate files
    system("rm -f cavgs_part*");
    # stash the final cavgs
    my $fbody = 'cavgs_iter'.$round_str.'*';
    system("cp $fbody ./$folder");
}

sub exec_check3D_conv{
    my$status;
    my@args = parse_args('check3D_conv');
    my$check_conv_cmd_string = $name_value{'check3D_conv'}.' ';
    foreach (@args) {
        $check_conv_cmd_string = add2string($_, $check_conv_cmd_string);
    }
    $check_conv_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $check_conv_cmd_string, "\n";
    my$conv_info = `$check_conv_cmd_string`;
    my@conv_info_lines = split(/\n/,$conv_info);
    $status = 0;
    foreach my$line (@conv_info_lines){
        print $line, "\n";
        if( $line =~ /\>\>\> PERCENTAGE OF SEARCH SPACE SCANNED\:\s+(\d+\.\d)/ ){
            $frac_srch_space = $1;
        }
        if( $line =~ /\>\>\> UPDATE LOW-PASS LIMIT\: \.YES\./ ){ $status = 1 };
        if( $line =~ /\>\>\> CONVERGED\: \.YES\./ )            { $status = 2 };
    }
    if( $status == 1 and $updated_res == 1 ){
        $status      = 0;
        $updated_res = 0; # to prevent resolution limit update in consequtive iterations
    }
    return $status;
}

sub exec_check2D_conv{
    my$status;
    $status = 0;
    my@args = parse_args('check2D_conv');
    my$check_conv_cmd_string = $name_value{'check2D_conv'}.' ';
    foreach (@args) {
        $check_conv_cmd_string = add2string($_, $check_conv_cmd_string);
    }
    $check_conv_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $check_conv_cmd_string, "\n";
    if( defined($name_value{'trs'}) ){
        if( $round_cnt <= 2 ){
            return $status;
        }
    }
    my$conv_info = `$check_conv_cmd_string`;
    my@conv_info_lines = split(/\n/,$conv_info);
    foreach my$line (@conv_info_lines){
        print $line, "\n";
        if( $line =~ /\>\>\> PERCENTAGE OF SEARCH SPACE SCANNED\:\s+(\d+\.\d)/ ){
            $frac_srch_space = $1;
        }
        if( $line =~ /\>\>\> CONVERGED\: \.YES\./ ) { $status = 2 };
    }
    return $status;
}

sub exec_merge_algndocs{
    my $new_oritab = shift;
    my $old_oritab;
    if( defined($name_value{'oritab'}) ){
        $old_oritab = $name_value{'oritab'};
        delete $name_value{'oritab'};
    }
    $name_value{'oritab'} = $new_oritab;
    my@args = parse_args('merge_algndocs');
    my$merge_algndocs_cmd_string = $name_value{'merge_algndocs'}.' ';
    foreach my$i (0.. $#args) {
        if( $args[$i] !~ 'oritab' ){
            $merge_algndocs_cmd_string = add2string($args[$i], $merge_algndocs_cmd_string);
        }
    }
    $merge_algndocs_cmd_string =~ s/\s+$//; # removes trailing whitespace
    if(defined($old_oritab)){
        if( -e $old_oritab ){
            $merge_algndocs_cmd_string = $merge_algndocs_cmd_string." outfile=".$new_oritab." oritab=".$old_oritab;
        }else{
            $merge_algndocs_cmd_string = $merge_algndocs_cmd_string." outfile=".$new_oritab;
        }
    } else {
        $merge_algndocs_cmd_string = $merge_algndocs_cmd_string." outfile=".$new_oritab;
    }
    print $merge_algndocs_cmd_string, "\n";
    system($merge_algndocs_cmd_string);
    system("rm $ALGNDOC_FBODY*");
    system("cp $primedoc ./$folder");
}

sub exec_merge_docs{
    my $fbody   = shift;
    my $outfile = shift;
    my $fbody_old;
    my $outfile_old;
    # store & remove old hash values
    if( defined($name_value{'fbody'}) ){
        $fbody_old = $name_value{'fbody'};
        delete $name_value{'fbody'};
    }
    if( defined($name_value{'outfile'}) ){
        $outfile_old = $name_value{'outfile'};
        delete $name_value{'outfile'};
    }
    # replace with the new args
    $name_value{'fbody'}   = $fbody;
    $name_value{'outfile'} = $outfile;
    # parse command line
    my@args = parse_args('merge_algndocs');
    my$merge_docs_cmd_string = $name_value{'merge_algndocs'}.' ';
    foreach my$i (0.. $#args) {
        $merge_docs_cmd_string = add2string($args[$i], $merge_docs_cmd_string);
    }
    $merge_docs_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $merge_docs_cmd_string, "\n";
    # execute
    system($merge_docs_cmd_string);
    # put back old values
    if( defined($fbody_old) ){
        delete $name_value{'fbody'};
        $name_value{'fbody'} = $fbody_old;
    }
    if( defined($outfile_old) ){
        delete $name_value{'outfile'};
        $name_value{'outfile'} = $outfile_old;
    }
}

sub exec_merge_similarities{
    my@args = parse_args('merge_similarities');
    my$merge_similarities_cmd_string = $name_value{'merge_similarities'}.' ';
    my $nptcls_tmp = exec_check_nptcls();
    my $nptcls_old = $name_value{'nptcls'};
    delete $name_value{'nptcls'};
    $name_value{'nptcls'} = $nptcls_tmp;
    foreach my$i (0.. $#args) {
        $merge_similarities_cmd_string = add2string($args[$i], $merge_similarities_cmd_string);
    }
    $merge_similarities_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $merge_similarities_cmd_string, "\n";
    delete $name_value{'nptcls'};
    $name_value{'nptcls'} = $nptcls_old;
    system($merge_similarities_cmd_string);
}

sub exec_split_pairs{
    my@args = parse_args('split_pairs');
    my$split_pairs_cmd_string = $name_value{'split_pairs'}.' ';
    my $nptcls_tmp = exec_check_nptcls();
    my $nptcls_old = $name_value{'nptcls'};
    delete $name_value{'nptcls'};
    $name_value{'nptcls'} = $nptcls_tmp;
    foreach my$i (0 .. $#args) {
        $split_pairs_cmd_string = add2string($args[$i], $split_pairs_cmd_string);
    }
    $split_pairs_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $split_pairs_cmd_string, "\n";
    delete $name_value{'nptcls'};
    $name_value{'nptcls'} = $nptcls_old;
    system($split_pairs_cmd_string);    
}

sub exec_check_nptcls{
    my$nptcls;
    if( defined($name_value{vollist}) ){
        open(VOLL, "<$name_value{vollist}") or die "Cannot open file: $name_value{vollist}, $!\n";
        while (<VOLL>) {}
        return $.;
        close(VOLL);
    }
    my$check_nptcls_cmd_string = $name_value{'check_nptcls'}.' ';
    $check_nptcls_cmd_string = add2string('stk', $check_nptcls_cmd_string);
    $check_nptcls_cmd_string = add2string('box', $check_nptcls_cmd_string);
    $check_nptcls_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $check_nptcls_cmd_string, "\n";
    my$nptcls_info = `$check_nptcls_cmd_string`;
    if( $nptcls_info =~ />>>\sNPTCLS\:\s+(\d+)/ ){
        $nptcls = $1;
    }
    return $nptcls;
}

sub exec_check_ncls{
    my$ncls;
    my$check_ncls_cmd_string = $name_value{'check_nptcls'}.' ';
    $check_ncls_cmd_string = add2string('refs', $check_ncls_cmd_string);
    $check_ncls_cmd_string = add2string('box', $check_ncls_cmd_string);
    $check_ncls_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $check_ncls_cmd_string, "\n";
    my$ncls_info = `$check_ncls_cmd_string`;
    if( $ncls_info =~ />>>\sNPTCLS\:\s+(\d+)/ ){
        $ncls = $1;
    }
    return $ncls;
}

sub exec_check_box{
    my$box;
    my$check_box_cmd_string = $name_value{'check_box'}.' ';
    $check_box_cmd_string = add2string('stk', $check_box_cmd_string);
    $check_box_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $check_box_cmd_string, "\n";
    my$box_info = `$check_box_cmd_string`;
    if( $box_info =~ />>>\sBOX\:\s+(\d+)/ ){
        $box = $1;
    }
    return $box;
}

sub exec_resrange{
    my$lpstart;
    my$lpstop;
    my$resrange_cmd_string = $name_value{'resrange'}.' ';
    $resrange_cmd_string = add2string('box',    $resrange_cmd_string);
    $resrange_cmd_string = add2string('smpd',   $resrange_cmd_string);
    $resrange_cmd_string = add2string('nspace', $resrange_cmd_string);
    $resrange_cmd_string = add2string('pgrp',   $resrange_cmd_string);
    $resrange_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $resrange_cmd_string, "\n";
    my$resrange_info = `$resrange_cmd_string`;
    if( $resrange_info =~ />>>\sLP\sSTART\:\s+(\d+)/ ){
        $lpstart = $1;
    }
    if( $resrange_info =~ />>>\sLP\sSTOP\:\s+(\d+)/ ){
        $lpstop = $1;
    }
    return ($lpstart,$lpstop);
}

sub add2string{
    my$var = shift;
    my$str = shift;
    if( defined($name_value{$var}) ){
        $str = $str.$var.'='.$name_value{$var}.' ';
    }
    return $str;
}

sub parse_args{
    my$prg = shift; # get program name with absolute path
    # make command line instructions
    my$instructions = `$name_value{'exec'} prg=$prg`;
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

sub cleanup{
    system("rm -rf FOO OUT* algndoc_* distr_script_* recvol_state*_part* rho* fort.0 primedoc_* recvol* JOB_FINISHED_* errfile.* outfile.* shdistr_script SHMEMJOBOUT shmemerrfile.* shmemoutfile.* prime2Ddoc_* ctfsqsums_part* noisespecs_part*");
}

sub get_find{
    my$res = shift;
    return int((($name_value{'box'}-1)*$name_value{'smpd'})/$res)
}

# sub get_lp{
#     my$find = shift;
#     return (($name_value{'box'}-1)*$name_value{'smpd'})/$find
# }

sub zero_pad_intg{
    my$intg   = shift;
    my$numlen = shift;
    while( length($intg) < length($numlen) ){
        $intg = '0'.$intg;
    }
    return $intg;
}
