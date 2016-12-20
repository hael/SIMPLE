use Term::ANSIColor;
use strict;
use Cwd qw(getcwd);

my $finishdir = getcwd();
finishing_up_compilation($finishdir);
################################################################################
# Subroutines                                                                  #
################################################################################
#subroutine to finish up script.
sub finishing_up_compilation{
    my ($finishdir_in) = @_;
    print color('bold blue');
    print"SIMPLE library has finished compilation in dir: ";print color('bold green');
    print"$finishdir_in\n"; print color('bold blue');
    print "*******************************************************************************\n";
    print "* Compilation and linkage is complete for Simple-v2.1                         *\n";
    print "* You may run all simple checks  --- List check options                       *\n";
    print "* ";print color('bold yellow');print"> make [OPTIONS] [VAR=VALUE]";
    print color('bold blue');print"   --- ";
    print color('bold yellow');print"> make check_help";
    print color('bold blue');print"                        *\n";
    print "* ";print color('bold yellow');print"            ";
    print color('bold blue');print"                   --- ";
    print color('bold yellow');print"> make check_news";
    print color('bold blue');print"                        *\n";
    print "* ";print color('bold yellow');print"[OPTIONS]:";
    print color('bold blue');print" --- ";
    print color('bold yellow');print"> {check, check_cpu, check_gpu, "; 
    print color('bold blue');print"                             *\n";
    print color('bold yellow'); print"                    ";
    print"bench, check_news, wc, tar}";
    print color('bold blue');print"                               *\n";
    print "* ";print color('bold yellow');print"[VAR]:";
    print color('bold blue');print"     --- ";
    print color('bold yellow');print"> {use_gpu,bench_gpu,fix_gpu,set_gpu,help}";
    print color('bold blue');print"                   *\n";
    print "* ";print color('bold yellow');print"[VALUE]:";
    print color('bold blue');print"   --- ";
    print color('bold yellow');print"> {yes,no},{0,..,MAX_N_GPU}";
    print color('bold blue');print"                                  *\n";
    print "* ";print color('bold yellow');print"[Example]:";
    print color('bold blue');print" --- ";
    print color('bold yellow');print"> make check use_gpu=yes help=no";
    print color('bold blue');print"                             *\n";
    print "* ";print color('bold yellow');print"          ";
    print color('bold blue');print" --- ";
    print color('bold yellow');print"> make bench use_gpu=yes fix_gpu=yes set_gpu=0 bench_gpu=no";
    print color('bold blue');print"  *\n";
    print "* ";print color('bold blue');print"         ";
    print color('bold blue');print"     ";
    print color('bold yellow');print"                                  ";
    print color('bold blue');print"                            *\n";
    print "* ";print color('bold blue');print"Cleaners:";
    print color('bold blue');print"  --- ";
    print color('bold yellow');print"> make {clean,cleanall,clean_check_cpu,"; 
    print color('bold blue');print"                      *\n";
    print color('bold yellow'); print"                         ";
    print"clean_check_gpu}";
    print color('bold blue');print"                                     *\n";
    print "* ";print color('bold blue');print"New Rel.:";
    print color('bold blue');print"  --- ";
    print color('bold yellow');print"> make check_news             ";
    print color('bold blue');print"                               *\n";
    print "* ";print color('bold blue');print"Lne Cntr:";
    print color('bold blue');print"  --- ";
    print color('bold yellow');print"> make wc             ";
    print color('bold blue');print"                                       *\n";
    print "*******************************************************************************\n";
    print color('reset');
}
