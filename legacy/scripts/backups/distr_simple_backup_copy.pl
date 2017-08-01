#!/usr/bin/perl
use lib './';
use lib '../';
use simple_ClusterDistr;
use warnings;
use strict;
use Cwd qw(getcwd);
use List::Util qw(max min);
# check so that at least one argument
if(scalar(@ARGV) == 0){die "need at least one argument: prg=<simple program name>\n"};
# global variables
my%name_value;
my$folder;
my$primedoc;
my$lpstart;
my$lpstop;
my$dynlp;
my$command;
my$prime_exec;
my$ext;
$name_value{'prime2'}         = $SIMPLEBIN.'/simple_prime2';
$name_value{'prime2_cluster'} = $SIMPLEBIN.'/simple_prime2_cluster';
$name_value{'merge_algndocs'} = $SIMPLEBIN.'/simple_merge_algndocs';
$name_value{'volassemble'}    = $SIMPLEBIN.'/simple_volassemble';
$name_value{'eo_volassemble'} = $SIMPLEBIN.'/simple_eo_volassemble';
$name_value{'cavgassemble'}   = $SIMPLEBIN.'/simple_cavgassemble';
$name_value{'automask'}       = $SIMPLEBIN.'/simple_automask';
$name_value{'check_conv'}     = $SIMPLEBIN.'/simple_check_conv';
$name_value{'check2d_conv'}   = $SIMPLEBIN.'/simple_check2d_conv';
$name_value{'check_nptcls'}   = $SIMPLEBIN.'/simple_check_nptcls';
$name_value{'check_box'}      = $SIMPLEBIN.'/simple_check_box';
$name_value{'resrange'}       = $SIMPLEBIN.'/simple_resrange';
$name_value{'recvol'}         = $SIMPLEBIN.'/simple_recvol';
$name_value{'eo_recvol'}      = $SIMPLEBIN.'/simple_eo_recvol';
$name_value{'fsc_filt'}       = $SIMPLEBIN.'/simple_fsc_filt';
$name_value{'extr_ptcls'}     = $SIMPLEBIN.'/simple_extr_ptcls';
$name_value{'execdir'}        = getcwd();
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
# check wether we are executing prime or not
$prime_exec = 0;
if( $name_value{'prg'} =~ /prime/ ){
    $prime_exec = 1;
}
if( $name_value{'prg'} =~ /oasis/ ){
    $prime_exec = 2;
}
if( $name_value{'prg'} =~ /cluster/ ){
    $prime_exec = 3;
}
# check so that prg is given with the simple_prefix
if( $name_value{'prg'} =~ /^simple_/ ){
    # all good!
}else{
    my $tmp = 'simple_'.$name_value{'prg'};
    delete $name_value{'prg'};
    $name_value{'prg'} = $tmp;
}
# make command line instructions
my$instructions = `$name_value{'prg'}`;
$instructions =~ s/SIMPLE_/DISTR_SIMPLE_/;
$instructions = $instructions."npart=<nr of partitions>\n"; #[nthr_shmem=<nr shared-mem threads>]\n";
# if too few argument given, print instructions
if(scalar(@ARGV) < 3){
    print $instructions;
}
# make absolute path to the program (needed on the node)
my $abs_prg = $SIMPLEBIN.'/'.$name_value{'prg'};
delete $name_value{'prg'};
$name_value{'prg'} = $abs_prg;
# check so that a stack is inputted
if( !defined($name_value{'filetab'}) ){
    if( !defined($name_value{'stk'}) ){
        die "Need a stack of particle images 4 distr_simple exec!\n";
    }
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
    $name_value{'nptcls'} = $nlines;
    # extract file-suffix
    @pieces = split(/\./,$first_movie_file);
    $ext = '.'.$pieces[-1];
}else{
    # extract file-suffix
    @pieces = split(/\./,$name_value{'stk'});
    $ext = '.'.$pieces[-1];
    if( $ext eq '.mrcs'){
        $ext = '.mrc;'
    }
    print "file extension of stack: $ext\n";
}
if( !defined($name_value{'filetab'}) ){
    # check so that box is set
    if( !defined($name_value{'box'}) ){
        # set the number of particles
        $name_value{'box'} = exec_check_box(); 
    }    
    # set memory requirements
    # if( $name_value{'box'} <= 200 ){
    #    $name_value{'mem'} = '10000MB';
    # }elsif( $name_value{'box'} <= 240 ){
    #    $name_value{'mem'} = '15000MB';
    # }elsif($name_value{'box'} <= 300){
    #    $name_value{'mem'} = '20000MB';
    # }elsif($name_value{'box'} <= 400){
    #    $name_value{'mem'} = '35000MB';
    # }else{
    #    $name_value{'mem'} = '40000MB';
    # }
    # set the number of particles
    $name_value{'nptcls'} = exec_check_nptcls();
    # check so that the mask parameter is set
    if( !defined($name_value{'msk'}) ){
        $name_value{'msk'} = $name_value{'box'}/2
    }
    if( $prime_exec == 3 ){
        if( !defined($name_value{'ncls'}) ){
            # set the number of particles
            $name_value{'ncls'} = exec_check_ncls();
        }
    }
# }else{
#     $name_value{'mem'} = '40000MB';
}
# check so that npart is set (number of partitions)
if( !defined $name_value{'npart'} ){
    die "npart need to be set\n";
}
# check so that ndocs is set (number of documents=number of partitions)
if( defined $name_value{'ndocs'} ){
    delete $name_value{'ndocs'};
}
$name_value{'ndocs'} = $name_value{'npart'};
# check so that pgrp is set
if( !defined($name_value{'pgrp'}) ){
    $name_value{'pgrp'} = 'c1'; 
}
# check so that nthr is set (number of threads)
# if( defined $name_value{'nthr'} ){
# }else{
#     if( $SIMPLESYS eq 'NEWMASSIVE' or $SIMPLESYS eq 'OXFORD' ){
#         $name_value{'nthr'} = 8;
#     }else{
#         $name_value{'nthr'} = 1;
#     }
# }
# check so that nthr_shmem is set (number of threads)
# if( defined($name_value{'nthr_shmem'}) ){
# }else{
#     $name_value{'nthr_shmem'} = $name_value{'nthr'};
# }
# check so that outvol is not defined
if( defined $name_value{'outvol'} ){
    die "Param: outvol cannot be defined! Used internally...";
}
# check so that file body is set
if( defined $name_value{'fbody'} ){
    delete $name_value{'fbody'}; 
}
$name_value{'fbody'} = 'algndoc_';
# check so that maxits is set 
if( !defined($name_value{'maxits'}) ){
    $name_value{'maxits'} = 100;
}
# check so that refine is set
if( !defined($name_value{'refine'})){
    if( $prime_exec == 1 or $prime_exec == 3 ){
        $name_value{'refine'} = 'no';
    }elsif($prime_exec == 2){
        $name_value{'refine'} = 'soft';   
    }
}
if( $prime_exec > 0 ){
    my@should_not_be_there = glob("recvol_state*");
    if( scalar(@should_not_be_there) >= 1 ){
        die "ERROR: files match the pattern recvol_state*, remove these from cwd\n";
    }
    if( $prime_exec < 2 ){
        # check so that a volume is inputted
        if( !defined($name_value{'vol1'}) ){
            die "Need at least one starting volume 4 exec!\n";
        }
    }elsif( $name_value{'prg'} =~ /cluster/ ){
        # check so that references are inputted
        if( !defined($name_value{'refs'}) ){
            die "Need starting references 4 prime2_cluster exec!\n";
        }
    }
    # check so that fstep is set
    if( !defined $name_value{'fstep'} ){
        $name_value{'fstep'} = 1;
    }
    # indicate whether diversify is defined in the command line args
    if( defined($name_value{'diversify'}) ){
        $name_value{'diversify_defined'} ='yes';
    }else{
        $name_value{'diversify_defined'} ='no';
    }
    if( $prime_exec == 1 ){
        # set diversify
        if( $name_value{'refine'} eq 'shc' ){
            if( defined($name_value{'diversify'}) ){
                delete $name_value{'diversify'};
            }
            $name_value{'diversify'} = 'no';
        }else{
            if( !defined($name_value{'oritab'}) ){
                if( defined($name_value{'diversify'}) ){
                    delete $name_value{'diversify'};
                }
                $name_value{'diversify'} = 'no';
            }
        }
    }
    if( $prime_exec == 3 ){
        if( !defined($name_value{'oritab'}) ){
            die "Need input clustering solution (oritab) for prime2_cluster\n"; 
        }
    }
    # check so that eo is set
    if( !defined($name_value{'eo'}) ){
        if( $prime_exec == 1 ){
            $name_value{'eo'} = 'no';
        }else{
            $name_value{'eo'} = 'yes';
        }
    }
    # if eo==yes, then let check_conv determine lp from fsc-file
    if( $name_value{'eo'} eq 'yes' ){
        $name_value{'fsc'} = 'fsc_state1.bin';
    }
    # set dynlp
    if( $prime_exec == 1 ){
        $dynlp = 'yes';
    }else{
        $dynlp = 'no';
    }
    if( defined($name_value{'dynlp'}) ){
        $dynlp = $name_value{'dynlp'};
    }elsif($name_value{'eo'} eq 'yes' or defined($name_value{'lp'}) ){
        $dynlp = 'no';
    }
    # set norec
    # $norec = 'no';
    # if( defined($name_value{'norec'}) ){
    #     $norec = $name_value{'norec'};
    # }
    # check so that nspace is set (number of projection directions)
    if( !defined $name_value{'nspace'} ){
        $name_value{'nspace'} = 1000;
    }
    # my $update_amsklp = 'yes';
    # # check if amsklp is defined
    # if( defined($name_value{'amsklp'}) ){
    #     $update_amsklp = 'no';
    # }
    if( $prime_exec == 1 and  $dynlp eq 'yes' ){
        # determine initial and final Fourier index
        if( $dynlp eq 'yes' ){
            ($lpstart, $lpstop) = exec_resrange();
            if( not defined($name_value{'lp'}) ){
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
            if( defined($name_value{'find'}) ){
                delete $name_value{'find'};
            }
            $name_value{'find'} = get_find($name_value{'lp'});
            # determine the final low-pass limit
            if( defined($name_value{'fmax'}) ){
                delete $name_value{'fmax'};
            }
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
    $name_value{'nstates'} = max 1,$nstates;
}
# check so that nstates is set 
if( !defined($name_value{'nstates'}) ){
    $name_value{'nstates'} = 1;
}
# check so that nspace is set 
# if( !defined($name_value{'nspace'}) ){
#     $name_value{'nspace'} = 1000;
# }
# set time per image
# if( !defined($name_value{'time_per_image'}) ){
#     $name_value{'time_per_image'} = 100*$name_value{'nstates'}*($name_value{'nspace'}/1000.);
# }
# set time per shmemjob
# if( !defined($name_value{'time_shmemjob'}) ){
#     $name_value{'time_shmemjob'} = 5000;
# }
# do the refinement
my$converged = 0;
my$round;
my$round_stop;
my$lp;
#my$amsklp;
my$refs;
if( defined($name_value{'startit'}) ){
    $round = $name_value{'startit'}-1; 
} else {
    $round = 0;
}
my$updated_res = 0;
if( $prime_exec > 0 ){
    while( $round < $name_value{'maxits'} ){
        # make refinement folder
        $round++;
        if( $prime_exec == 1 ){
            $folder   = 'prime_round_'.$round;
            $primedoc = 'primedoc_'.$round.'.txt';
        }elsif( $prime_exec == 2 ){
            $folder   = 'oasis_round_'.$round;
            $primedoc = 'oasisdoc_'.$round.'.txt';
        }elsif( $prime_exec == 3 ){
            $folder   = 'prime2d_round_'.$round;
            $primedoc = 'prime2ddoc_'.$round.'.txt';
            $refs     = 'cavgs_iter'.$round.'msk'.$ext;
            if( defined($name_value{'which_iter'}) ){
                delete $name_value{'which_iter'};
            }
            $name_value{'which_iter'} = $round;
        }
        mkdir $folder;
        # parallell execution
        if( $prime_exec == 1 and $name_value{'diversify_defined'} eq 'no' and $name_value{'refine'} ne 'shc' ){
            if( defined($name_value{'diversify'}) ){
                delete $name_value{'diversify'};
                $name_value{'diversify'} = 'yes';
            }
        }
        system("rm -f JOB_FINISHED_*");
        exec_para($name_value{'prg'});
        # determine when the jobs have finished
        sleep(10) while( njobs_finished() < $name_value{'npart'} );
        system("rm -f JOB_FINISHED_* OUT* distr_script*");
        # merge alignment docs
        exec_merge_algndocs($primedoc);
        if( $prime_exec <= 2 ){
            # reconstruct/assemble volumes
            if( $name_value{'eo'} eq 'no' ){
                # if( $norec eq 'yes' ){
                #     exec_para($name_value{'recvol'});
                #     # determine when the jobs have finished
                #     sleep(10) while( njobs_finished() < $name_value{'npart'} );
                #     system("rm -f JOB_FINISHED_*");
                # }
                exec_assemble_volumes();
            # }elsif( ($norec eq 'yes') and ($name_value{'eo'} eq 'yes') ){
            #     if( $name_value{'nstates'} > 1 ){
            #         die "Lowmem rec not available for nstates > 1";
            #     }
            #     exec_eorec_lowmem('even');
            #     exec_eorec_lowmem('odd');
            #     exec_eoassemble_lowmem();
            #     ($lp, $amsklp) = exec_fsc_filt();
            }else{
                # if( $norec eq 'yes' ){
                #     exec_para($name_value{'eo_recvol'});
                #     # determine when the jobs have finished
                #     sleep(10) while( njobs_finished() < $name_value{'npart'} );
                #     system("rm -f JOB_FINISHED_*");
                # }
                $lp = exec_assemble_eo_volumes();
            }
            # automask
            if( defined($name_value{'mw'}) or defined($name_value{'nvox'}) ){
                exec_automask();
            }
        }elsif( $prime_exec == 3 ){
            if( defined($name_value{'refs'}) ){
                delete $name_value{'refs'};
            }
            $name_value{'refs'} = $refs;            
            exec_assemble_cavgs();
        }else{
            die "unknown prime_exec mode; distr_simple.pl\n";
        }
        # check convergence
        my $conv;
        if( $prime_exec == 3 ){
            $conv = exec_check2d_conv();
        }else{
            $conv = exec_check_conv();
        }
        if( $round > 1 ){
            if( $conv == 2 ){
                die "**** DISTR_SIMPLE_PRIME NORMAL STOP ****\n";
            }elsif ( $conv == 1 ){
                if( !defined($name_value{'find'}) ){
                    die "Want 2 do dynamic low-pass update, but find is not defined, weird!\n";
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
    exec_para($name_value{'prg'});
}

sub exec_para{
    my$prog = shift;
    my$i;
    my$vol;
    my@submitted;
    my$subout;
    my@args = parse_args($prog);
    my$cmd_string = $prog.' ';
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
    $cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $cmd_string, "\n";
    my $stop=0;
    my $start;
    my $ptcls_per_part = int($name_value{'nptcls'}/$name_value{'npart'});
    my$leftover = $name_value{'nptcls'}-$ptcls_per_part*$name_value{'npart'};
    for(my $i=1; $i<=$name_value{'npart'}; $i++){
        if( $i == $name_value{'npart'} ){
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
        my $time  = $name_value{'time_per_image'}*($stop-$start+1);
        my $hours = int($time/3600);
        my $min   = (($time/60)%60);
        open(FHANDLE, ">distr_script_$i") or die "Cannot open distr_script_$i for writing: $!\n";
        print FHANDLE generate_distr_script($hours, $min, $name_value{'execdir'}, $cmd_string, $start, $stop, $i);
#         if ($SIMPLESYS eq 'BIOX3' or $SIMPLESYS eq 'LOCAL'){
#             print FHANDLE "#!/bin/bash
# #PBS -N distr_simple
# #PBS -l nodes=1:ppn=1
# #PBS -l walltime=$time
# #PBS -o outfile.\$PBS_JOBID
# #PBS -e errfile.\$PBS_JOBID
# #PBS -q SP
# cd $name_value{'execdir'}
# $cmd_string fromp=$start top=$stop part=$i outfile=$name_value{'fbody'}$i.txt > OUT$i\nexit\n";
#         }elsif ($SIMPLESYS eq 'DELPHI'){
#             print FHANDLE "#!/bin/bash
# #\$ \-N d_simple
# #\$ \-cwd
# #\$ \-pe multi 1
# #\$ \-V
# #\$ \-o outfile.\$PBS_JOBID
# #\$ \-e errfile.\$PBS_JOBID
# cd $name_value{'execdir'}
# $cmd_string fromp=$start top=$stop part=$i outfile=$name_value{'fbody'}$i.txt > OUT$i\nexit\n";
#         }elsif($SIMPLESYS eq 'MASSIVE'){
#             print FHANDLE "#!/bin/bash
# #PBS -N distr_simple
# #PBS -l nodes=1:ppn=1,mem=$name_value{'mem'}
# #PBS -l walltime=$hours:$min:0
# #PBS -o outfile.\$PBS_JOBID
# #PBS -e errfile.\$PBS_JOBID
# cd $name_value{'execdir'}
# $cmd_string fromp=$start top=$stop part=$i outfile=$name_value{'fbody'}$i.txt > OUT$i\nexit\n";
#         }elsif ($SIMPLESYS eq 'OXFORD' ){
#             print FHANDLE "#!/bin/bash
# #PBS -N distr_simple
# #PBS -l nodes=1:ppn=$name_value{'nthr'},mem=9gb
# #PBS -l walltime=$hours:$min:0
# #PBS -o outfile.\$PBS_JOBID
# #PBS -e errfile.\$PBS_JOBID
# #PBS -V
# #PBS -l naccesspolicy=UNIQUEUSER
# cd $name_value{'execdir'}
# mpirun -np 1 --bind-to-socket --cpus-per-proc 8 $cmd_string fromp=$start top=$stop part=$i outfile=$name_value{'fbody'}$i.txt > OUT$i\nexit\n";
#         }elsif ($SIMPLESYS eq 'NEWMASSIVE' ){
#             print FHANDLE "#!/bin/bash
# #SBATCH --mail-user=<hans\.elmlund\@monash\.edu>
# #SBATCH --mail-type=FAIL
# #SBATCH --job-name=distr_simple
# #SBATCH --ntasks=1
# #SBATCH --ntasks-per-socket=1
# #SBATCH --cpus-per-task=$name_value{'nthr'}
# #SBATCH --mem=32000
# #SBATCH --time=0-$hours:$min:0
# #SBATCH --output=outfile.%j
# #SBATCH --error=errfile.%j
# #SBATCH --partition=cryoem
# #SBATCH --qos=vip_m2
# cd $name_value{'execdir'}
# $cmd_string fromp=$start top=$stop part=$i outfile=$name_value{'fbody'}$i.txt > OUT$i\nexit\n";
#         }elsif ($SIMPLESYS eq 'CVL' ){
#             print FHANDLE "#!/bin/bash
# #SBATCH --mail-user=<hans\.elmlund\@monash\.edu>
# #SBATCH --mail-type=FAIL
# #SBATCH --job-name=distr_simple
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=16
# #SBATCH --mem=32000
# #SBATCH --time=0-$hours:$min:0
# #SBATCH --output=outfile.%j
# #SBATCH --error=errfile.%j
# #SBATCH --partition=compute
# #SBATCH --account=cvl
# cd $name_value{'execdir'}
# $cmd_string fromp=$start top=$stop part=$i outfile=$name_value{'fbody'}$i.txt > OUT$i\nexit\n";
#         } else {
#             die "Unrecognized system!\n";
#         }
        close(FHANDLE);
        chmod 0777, 'distr_script_'.$i;
        if( $SIMPLESYS eq 'LOCAL' ){
            system("nohup ./distr_script_$i &\n");
        # }elsif( $SIMPLESYS eq 'BIOX3' ){
        #     my$subcmd = "qsub ./distr_script_$i";
        #     $subout = `$subcmd`;
        #     # nothing else 4 now
        # }elsif( $SIMPLESYS eq 'DELPHI' ){
        #     my$subcmd = "qsub ./distr_script_$i";
        #     $subout = `$subcmd`;
        #     # nothing else 4 now
        # }elsif($SIMPLESYS eq 'MASSIVE' or $SIMPLESYS eq 'OXFORD'){
        #     my$subcmd = "qsub ./distr_script_$i";
        #     $subout = `$subcmd`;
        #     # nothing else 4 now
        # }elsif($SIMPLESYS eq 'NEWMASSIVE' or $SIMPLESYS eq 'CVL'){
        #     my$subcmd = "sbatch ./distr_script_$i";
        #     $subout = `$subcmd`;
        #     # nothing else 4 now
        }else{
            # die "Unrecognized system!\n";
            my$subcmd = "$SUBMITCMD ./distr_script_$i";
            $subout = `$subcmd`;
            # nothing else 4 now
        }
    }
    # if( $SIMPLESYS eq 'BIOX3' ){
#         # resubmission if failed
#         for(my $i=1; $i<=$name_value{'npart'}; $i++){
#             if( $submitted[$i] == 0 ){
#                 my$subcmd = "qsub ./distr_script_$i";
#                 $subout = `$subcmd`;
#                 if( $subout !~ /\d+\.biox3\-frontend\-1\.stanford\.edu/ ){
#                     chomp($subout);
#                     print "ERROR, $subout\n";
#                 }
#             }
#         }
#     } elsif( $SIMPLESYS eq 'DELPHI' ){
#         # resubmission if failed
#         for(my $i=1; $i<=$name_value{'npart'}; $i++){
#             if( $submitted[$i] == 0 ){
#                 my$subcmd = "qsub ./distr_script_$i";
#                 $subout = `$subcmd`;
#                 if( $subout !~ /has been submitted/ ){
#                     chomp($subout);
#                     print "ERROR, $subout\n";
#                 }
#             }
#         }
#     } elsif($SIMPLESYS eq 'MASSIVE' or $SIMPLESYS  eq 'OXFORD'){
#         # nothing 4 now
#     }
}

sub exec_shmem_para{
    my$prog = shift;
    my$i;
    my$vol;
    my$subout;
    my@args = parse_args($prog);
    my$cmd_string = $prog.' ';
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
    $cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $cmd_string, "\n";
    my $time  = $name_value{'time_shmemjob'};
    my $hours = int($time/3600);
    my $min   = (($time/60)%60);
    my $numeric_mem;
    # if( $name_value{'mem'} =~ /(\d+)/ ){
    #     $numeric_mem = $1;
    # }
    # $numeric_mem = int($numeric_mem/$name_value{'nthr_shmem'});
    # my $mem = $numeric_mem.'MB';
    # my $mem_per_cpu = $mem;
    open(FHANDLE, ">shmemdistr_script") or die "Cannot open shmemdistr_script for writing: $!\n";
    print FHANDLE generate_shmem_distr_script($hours, $min, $name_value{'execdir'}, $cmd_string);
#     if( $SIMPLESYS eq 'MASSIVE' or $SIMPLESYS eq 'LOCAL' ){
#         print FHANDLE "#!/bin/bash
# #PBS -N shmemdistr_simple
# #PBS -l nodes=1:ppn=$name_value{'nthr_shmem'},mem=$name_value{'mem'}
# #PBS -l walltime=$hours:$min:0
# #PBS -o outfile.\$PBS_JOBID
# #PBS -e errfile.\$PBS_JOBID
# cd $name_value{'execdir'}
# $cmd_string > SHMEMJOBOUT\nexit\n";
#     } elsif( $SIMPLESYS eq 'OXFORD'){
#         print FHANDLE "#!/bin/bash
# #PBS -N shmemdistr_simple
# #PBS -l nodes=1:ppn=$name_value{'nthr_shmem'},,mem=9gb
# #PBS -l walltime=$hours:$min:0
# #PBS -o outfile.\$PBS_JOBID
# #PBS -e errfile.\$PBS_JOBID
# #PBS -V
# #PBS -l naccesspolicy=UNIQUEUSER
# cd $name_value{'execdir'}
# mpirun -np 1 --bind-to-socket --cpus-per-proc 8 $cmd_string > SHMEMJOBOUT\nexit\n";
#     } elsif( $SIMPLESYS eq 'NEWMASSIVE' ){
#         print FHANDLE "#!/bin/bash
# #SBATCH --mail-user=<hans\.elmlund\@monash\.edu>
# #SBATCH --mail-type=FAIL
# #SBATCH --job-name=shmemdistr_simple
# #SBATCH --ntasks=1
# #SBATCH --ntasks-per-socket=1
# #SBATCH --cpus-per-task=$name_value{'nthr_shmem'}
# #SBATCH --mem=32000
# #SBATCH --time=0-$hours:$min:0
# #SBATCH --output=shmemoutfile.%j
# #SBATCH --error=shmemerrfile.%j
# #SBATCH --partition=cryoem
# #SBATCH --qos=vip_m2
# cd $name_value{'execdir'}
# $cmd_string > SHMEMJOBOUT\nexit\n";
#     } elsif( $SIMPLESYS eq 'CVL' ){
#         print FHANDLE "#!/bin/bash
# #SBATCH --mail-user=<hans\.elmlund\@monash\.edu>
# #SBATCH --mail-type=FAIL
# #SBATCH --job-name=shmemdistr_simple
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=16
# #SBATCH --mem=32000
# #SBATCH --time=0-$hours:$min:0
# #SBATCH --output=shmemoutfile.%j
# #SBATCH --error=shmemerrfile.%j
# #SBATCH --partition=compute
# #SBATCH --account=cvl
# cd $name_value{'execdir'}
# $cmd_string > SHMEMJOBOUT\nexit\n";
#     } else {
#         die "Unrecognized system!\n";
#     }
    close(FHANDLE);
    chmod 0777, 'shmemdistr_script';
    if( $SIMPLESYS eq 'LOCAL' ){
        system("nohup ./shmemdistr_script &\n");
    #     }elsif($SIMPLESYS eq 'MASSIVE' or $SIMPLESYS eq 'OXFORD'){
    #         my$subcmd = "qsub ./shmemdistr_script";
    #         $subout = `$subcmd`;
    #         # nothing else 4 now
    # }elsif($SIMPLESYS eq 'NEWMASSIVE' or $SIMPLESYS eq 'CVL'){
    #         my$subcmd = "sbatch ./shmemdistr_script";
    #         $subout = `$subcmd`;
    #         # nothing else 4 now
    }else{
        my$subcmd = "$SUBMITCMD ./shmemdistr_script";
        $subout = `$subcmd`;
        # nothing else 4 now 
    }
}

sub njobs_finished{
    my$nfini = 0;
    my@jobs;
    @jobs = <JOB_FINISHED_*>;
    return scalar(@jobs);     
}

sub shmemjob_finished{
    my $state = 0;
    if( -e 'SHMEMJOBOUT' ){
        open(CHECKHANDLE, 'SHMEMJOBOUT') or die "Cannot open SHMEMJOBOUT: $!\n";
        while(my$line=<CHECKHANDLE>){
            if( $line =~ /NORMAL\sSTOP\s\*\*\*\*/ ){ $state = 1 };
        }
        close(CHECKHANDLE);
    }
    return $state;
}

# sub exec_eorec_lowmem{
#     my$which = shift;
#     if( ($which eq 'even') or ($which eq 'odd') ){
#         if( defined($name_value{$which}) ){
#             delete $name_value{$which};
#         }
#         $name_value{$which} = 'yes';
#         exec_para($name_value{'recvol'});
#         # determine when the jobs have finished
#         sleep(10) while( njobs_finished() < $name_value{'npart'} );
#         system("rm -f JOB_FINISHED_*");
#         delete $name_value{$which};
#     }else{
#         die "Unknown which parameter; eorec_lowmem";
#     }
# }

# sub exec_eoassemble_lowmem{
#     if( defined($name_value{'even'}) ){
#         delete $name_value{'even'};
#     }
#     if( defined($name_value{'odd'}) ){
#         delete $name_value{'odd'};
#     }
#     $name_value{'even'} = 'yes';
#     exec_assemble_volumes();
#     delete $name_value{'even'};
#     $name_value{'odd'} = 'yes';
#     exec_assemble_volumes();
#     delete $name_value{'odd'};
#     exec_assemble_volumes();
# }

# sub exec_fsc_filt{
#     my$i;
#     my$vol;
#     my$lp;
#     my$amsklp;
#     my$fsc_filt_cmd_string = $name_value{'fsc_filt'}.' ';
#     $fsc_filt_cmd_string = $fsc_filt_cmd_string.' vol1=recvol_state1'.$ext;
#     $fsc_filt_cmd_string = $fsc_filt_cmd_string.' vol2=recvol_state1_even'.$ext;
#     $fsc_filt_cmd_string = $fsc_filt_cmd_string.' vol3=recvol_state1_odd'.$ext.' ';
#     $fsc_filt_cmd_string = $fsc_filt_cmd_string.' outvol=recvol_state1'.$ext.' ';
#     my@args = parse_args($name_value{'fsc_filt'});
#     foreach (@args) {
#         $fsc_filt_cmd_string = add2string($_, $fsc_filt_cmd_string);
#     }
#     $fsc_filt_cmd_string =~ s/\s+$//; # removes trailing whitespace
#     print $fsc_filt_cmd_string, "\n";
#     # execute assemble command and parse low-pass & mask low-pass
#     my$fsc_filt_info = `$fsc_filt_cmd_string`;
#     my@fsc_filt_info_lines = split(/\n/,$fsc_filt_info);
#     foreach my$line (@fsc_filt_info_lines){
#         print $line, "\n";
#         if( $line =~ /\>\>\> LOW-PASS LIMIT\:\s+(\d+.\d+)/ )     { $lp     = $1 };
#         if( $line =~ /\>\>\> MASK LOW-PASS LIMIT\:\s+(\d+.\d+)/ ){ $amsklp = $1 };
#     }
#     $vol = 'vol1';
#     if( defined($name_value{$vol}) ){
#         delete $name_value{$vol};
#     }
#     $name_value{$vol} = 'recvol_state1'.$ext;
#     system("cp $name_value{$vol} ./$folder");
#     system("cp fsc_state*bin ./$folder");
#     return ($lp, $amsklp);
# }

sub exec_assemble_volumes{
    my$i;
    my$vol;
    my@refvols;
    my@recvols;
    my@args = parse_args($name_value{'volassemble'});
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
    if( ($name_value{'odd'} eq 'yes') or ($name_value{'even'} eq 'yes') ){
        # do nothing
    }else{
        # delete intermediate files
        foreach my$i (0 .. $#kernels){
            system("rm -f $kernels[$i] $recbins[$i]");
        }
        # take care of the reconstructed volumes
        @recvols = <recvol_state*$ext>;
        chomp(@recvols);
        @recvols = grep {$_ =~ /^recvol_state\d+$ext$/} @recvols;
        my $nrecvols = scalar(@recvols);
        if( $nrecvols != $name_value{'nstates'} ){
            print "nstates: $name_value{'nstates'}\n";
            print "nrecvols: $nrecvols\n";
            print "ext: $ext\n";
            print "CONTENT IN RECVOLS ARRAY:";
            foreach(@recvols){
                print $_, "\n";
            }
            die "Number of reference volumes not equal to number of states";   
        } else {
            @refvols = @recvols;
        }
        foreach $i ( 1 .. $name_value{'nstates'} ){
            $vol = 'vol'.$i;
            if( defined($name_value{$vol}) ){
                delete $name_value{$vol};
            }
            $name_value{$vol} = $refvols[$i-1];
            system("cp $name_value{$vol} ./$folder");
        }
    }
}

sub exec_assemble_eo_volumes{
    my$i;
    my$vol;
    my@refvols;
    my@recvols;
    my$lp;
    my@args = parse_args($name_value{'eo_volassemble'});
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
    # execute assemble command and parse low-pass & mask low-pass
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
    }
    # take care of the reconstructed volumes
    @recvols = <recvol_state*$ext>;
    chomp(@recvols);
    @recvols = grep {$_ !~ /^recvol_state\d+_even.+$/} @recvols;
    @recvols = grep {$_ !~ /^recvol_state\d+_odd.+$/} @recvols;
    @recvols = grep {$_ =~ /^recvol_state\d+$ext$/} @recvols;    
    my $nrecvols = scalar(@recvols);
    if( $nrecvols != $name_value{'nstates'} ){
        print "nstates: $name_value{'nstates'}\n";
        print "nrecvols: $nrecvols\n";
        die "Number of reference volumes not equal to number of states";   
    } else {
        @refvols = @recvols;
    }
    foreach $i (1 .. $name_value{'nstates'}){
        $vol = 'vol'.$i;
        if( defined($name_value{$vol}) ){
            delete $name_value{$vol};
        }
        $name_value{$vol} = $refvols[$i-1];
        system("cp $name_value{$vol} ./$folder");
        system("cp fsc_state*bin ./$folder");
    }
    return $lp;
}

sub exec_assemble_cavgs{
    my@args = parse_args($name_value{'cavgassemble'});
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
    my $fbody = 'cavgs_iter'.$round.'*';
    system("cp $fbody ./$folder");
}

sub exec_automask{
    exec_shmem_para($name_value{'automask'});
    # determine when the job has finished
    sleep(10) while( shmemjob_finished() == 0 );
    system("rm SHMEMJOBOUT");
    my$fbody      = $name_value{'vol1'};
    $fbody        =~ s/\d+$ext//;
    my@maskedvols = <$fbody*msk$ext>;
    foreach my $i ( 1 .. $name_value{'nstates'} ){
        my $vol = 'vol'.$i;
        system("cp $maskedvols[$i-1] ./$folder");
        system("mv $maskedvols[$i-1] $name_value{$vol}");
    }
}

sub exec_check_conv{
    my$status;
    my@args = parse_args($name_value{'check_conv'});
    my$check_conv_cmd_string = $name_value{'check_conv'}.' ';
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
        if( $line =~ /\>\>\> UPDATE LOW-PASS LIMIT\: \.YES\./ ){ $status = 1 };
        if( $line =~ /\>\>\> CONVERGED\: \.YES\./ )            { $status = 2 };
    }
    if( $status == 1 and $updated_res == 1 ){ 
        $status = 0;
        $updated_res = 0; # to prevent resolution limit update in consequtive iterations
    }
    return $status;
}

sub exec_check2d_conv{
    my$status;
    my@args = parse_args($name_value{'check2d_conv'});
    my$check_conv_cmd_string = $name_value{'check2d_conv'}.' ';
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
    my@args = parse_args($name_value{'merge_algndocs'});
    my$merge_algndocs_cmd_string = $name_value{'merge_algndocs'}.' ';
    foreach my$i (0.. $#args-1) {
        $merge_algndocs_cmd_string = add2string($args[$i], $merge_algndocs_cmd_string);
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
    system("rm $name_value{'fbody'}*");
    system("cp $primedoc ./$folder");
}

sub exec_check_nptcls{
    my$nptcls;
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
    my$instructions = `$prg`;
    # parse the instructions
    my@splitted = split(/\s/,$instructions);
    my@args;
    foreach (@splitted) {
        if( $_ =~ /(\w+\d*)\=/ ){
            if( $1 !~ /vol\d+/ ){
                push(@args,$1);
            }  
        }
    }
    return @args;
}

sub get_find{
    my$res = shift;
    return int((($name_value{'box'}-1)*$name_value{'smpd'})/$res)
}

sub get_lp{
    my$find = shift;
    return (($name_value{'box'}-1)*$name_value{'smpd'})/$find
}