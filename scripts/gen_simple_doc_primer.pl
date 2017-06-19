#!/usr/bin/perl
use lib './';
use warnings;
use strict;
use Env;
use simple_user_input;
my $maxlen = 80;           # maximum string length of shell commands
my $maxlc  = 53;           # maximum line count (LaTEX table)
my @cmdlinkeys;            # command line keys
my %cmdlindict;            # command line dictonary
my %prginstr;              # program instructions
my %prginstr_distr;        # distributed program instructions
my %prgcmdlineinstr;       # program command line instructions
my %prgcmdlineinstr_distr; # distributed program command line instructions
my @prgkeys;               # program keys
my @prgkeys_distr;         # distributed program keys
if( scalar(@ARGV) < 1 ){
    die "Need to know which kind of documentation to generate (tex|web|f90)\n";
}
my $doc = $ARGV[0];
my $simple_exec_abspath = "$SIMPLE_PATH/production/simple/simple_exec/simple_exec.f90";
my $simple_distr_exec_abspath = "$SIMPLE_PATH/production/simple/simple_distr_exec/simple_distr_exec.f90";

# PARSERS

# extract program names
my @prgnames;
open( EXEC, "<$simple_exec_abspath", ) or die "Cannot open: $!\n";
while(my$line=<EXEC>){
    if( $line =~ /case\(\s*\'(\w+)\'\s*\)/ ){
        push(@prgnames,$1);
    }
}
close(EXEC);
my @prgnames_distr;
open( EXEC, "<$simple_distr_exec_abspath", ) or die "Cannot open: $!\n";
while(my$line=<EXEC>){
    if( $line =~ /case\(\s*\'(\w+)\'\s*\)/ ){
        push(@prgnames_distr,$1);     
    }
}
close(EXEC);

# extract comments
my $comments = extract_comments($simple_exec_abspath);
chomp($comments);
my $comments_distr = extract_comments($simple_distr_exec_abspath);
chomp($comments_distr);

# extract program instructions
my @comments_split = split('==', $comments);
my $cnt = 0;
foreach my $comment (@comments_split){
    if( $comment =~ /\s*\<(.+)\/begin\>((.+))\<.+\/end\>/ ){
        $cnt++;
        $prginstr{$1} = $2;
        push(@prgkeys,$1);
    }
}
@comments_split = split('==', $comments_distr);
$cnt = 0;
foreach my $comment (@comments_split){
    if( $comment =~ /\s*\<(.+)\/begin\>((.+))\<.+\/end\>/ ){
        $cnt++;
        $prginstr_distr{$1} = $2;
        push(@prgkeys_distr,$1);
    }
}

# generate command-line dictonary
my $foo     = `$SIMPLE_PATH/bin/simple_exec prg=print_cmd_dict outfile=cmd_dict.txt`;
my $nparams = 0;
open(CMDDICT, "<cmd_dict.txt") or die "Cannot open file: cmd_dict.txt for reading, $!\n";
while(my $keyval=<CMDDICT>){
    chomp($keyval);
    if( $keyval =~ /^(.+)\=((.+))$/ ){
        push(@cmdlinkeys,$1);
        my $key = str2latex($1);
        my $val = str2latex($2);
        $cmdlindict{$key} = $val;
        $nparams++; # number of command line parameters
    }
}

# generate per-program command line dictonaries
foreach my $prg (@prgnames){
    $prgcmdlineinstr{$prg} = `$SIMPLE_PATH/bin/simple_exec prg=$prg`;
    $prgcmdlineinstr{$prg} =~ s/\n$//; # removes trailing newline character
}
foreach my $prg (@prgnames_distr){
    $prgcmdlineinstr_distr{$prg} = `$SIMPLE_PATH/bin/simple_distr_exec prg=$prg`;
    $prgcmdlineinstr_distr{$prg} =~ s/\n$//; # removes trailing newline character
}

# CODE GENERATORS
my @prgnames_sorted       = sort @prgnames;
my @prgnames_distr_sorted = sort @prgnames_distr;
if( $doc eq 'tex' ){
    print_full_latex_cmdlindict();
    print '\section{Distributed SIMPLE Workflows}'."\n";
    foreach my $prg (@prgnames_distr_sorted){
        print_latex_instr(1, $prg, %prginstr_distr);
        print_latex_usage($prg, %prgcmdlineinstr_distr);
    }
    print '\section{SIMPLE Programs}'."\n";
    foreach my $prg (@prgnames_sorted){
        print_latex_instr(0, $prg, %prginstr);
        print_latex_usage($prg, %prgcmdlineinstr);
    }  
}elsif( $doc eq 'web' ){
    foreach my $prg (@prgnames_sorted){
        print_html_instr($prg, %prginstr);
    }
    foreach my $prg (@prgnames_distr_sorted){
        my$isthere = 0;
        foreach my $prg_tmp (@prgnames_sorted){
            if( $prg =~ /$prg_tmp$/ ){
                $isthere = 1;
            }
        }
        if( $isthere == 0 ){
            print_html_instr($prg, %prginstr_distr);
        }        
    }
}elsif( $doc eq 'f90' ){
    print "module simple_gen_doc\n";
    print "implicit none\n\n";
    print "contains\n\n";
    foreach my $prg (@prgnames_sorted){
        prginstr2f90($prg, %prginstr);
    }
    foreach my $prg (@prgnames_distr_sorted){
        my$isthere = 0;
        foreach my $prg_tmp (@prgnames_sorted){
            if( $prg =~ /$prg_tmp$/ ){
                $isthere = 1;
            }
        }
        if( $isthere == 0 ){
            prginstr2f90($prg, %prginstr_distr);
        }
    }
    print "    subroutine list_all_simple_programs\n";
    foreach my $prg (@prgnames_sorted){
        print "        write(*,'(A)') '$prg'\n";
    }
    print "        stop\n";
    print "    end subroutine list_all_simple_programs\n\n";
    print "    subroutine list_all_simple_distr_programs\n";
    foreach my $prg (@prgnames_distr_sorted){
        print "        write(*,'(A)') '$prg'\n";
    }
    print "        stop\n";
    print "    end subroutine list_all_simple_distr_programs\n\n";

    print "end module simple_gen_doc\n";
}else{
    die "Unknown doc type inputted!\n";
}

# DRIVERS

sub print_full_latex_cmdlindict{
    print '\section{Command Line Dictionary}', "\n";
    print '\begin{tabular}{ll}', "\n";
    my $linecnt = 0;
    foreach my $key (sort keys %cmdlindict) {
        $linecnt++;
        print '\texttt{'.$key.'}&{'.$cmdlindict{$key}.'}\\';
        print '\\', "\n";     
        if( $linecnt == $maxlc ){
            print '\end{tabular}', "\n\n";
            print '\begin{tabular}{ll}', "\n";
            $linecnt = 0;
        }
    }
    print '\end{tabular}', "\n\n";
}

sub print_latex_instr{
    # input is label and prginstr
    my $distr = shift;
    my $label = shift;
    my %instr = @_;
    my $label_latex = $label;
    $label_latex =~ s/\_/\\_/g;
    if( $distr == 0 ){
        print '\subsection{Program: \prgname{'."$label_latex}}\n";
    }else{
        print '\subsection{Distributed Workflow: \prgname{'."$label_latex}}\n";
    }    
    print '\label{'."$label}\n";
    if( defined($instr{$label}) ){
        my $prginstr_latex = str2latex($instr{$label});
        print '\prgname{'."$label_latex} ".$prginstr_latex."\n\n";
    }
}

sub print_latex_usage{
    # input is label and prgcmdlineinstr
    my $label        = shift;
    my %cmdlineinstr = @_;
    print '\begin{verbatim}', "\n";
    print $cmdlineinstr{$label};
    print '\end{verbatim}', "\n\n";
}

sub print_html_instr{
    # input is label and prginstr
    my $label = shift;
    my %instr = @_;
    if( defined($instr{$label}) ){
        my $instr_html = latex2html(str2latex($instr{$label}));
        my $prg2print = '************'.$label.'************';
        print "$prg2print\n\n";
        print qq[<div class="AccordionPanel">\n];
        print qq[  <div class="AccordionPanelTab">$label</div>\n];
        print qq[  <div class="AccordionPanelContent">\n];
        print qq[    <p><strong>$label</strong> $instr_html</p>\n];
        print qq[</div>\n];
        print qq[</div>\n];
        print "\n\n";
    }
}

sub str2latex{
    my $str2convert = shift;
    $str2convert =~ s/\{/\\{/g;
    $str2convert =~ s/\}/\\}/g;
    $str2convert =~ s/\_/\\_/g;
    $str2convert =~ s/%/\\%/g;
    $str2convert =~ s/#/\\#/g;
    $str2convert =~ s/in A/in \\AA\{\}/;
    return $str2convert;
}

sub latex2html{
    my $str = shift;
    if( defined($str) ){
        # remove any maths
        $str =~ s/\$.+\$//g;
        # replace the special LaTex char \&
        $str =~ s/\\&/\&amp;/g;
        # replace the special LaTex char \_
        $str =~ s/\\_/_/g;
        # replace the special LaTex char \%
        $str =~ s/\\%/%/g;
        # replace < (less than)
        $str =~ s/\</\&lt;/g;
        # replace > (greater than)
        $str =~ s/\>/\&gt;/g;
        # replace terminal LaTex font with italic html
        $str =~ s/\\texttt\{(\w+\.*\=*\w*\.*\w*)\}/<em>$1<\/em>/g;
        # remove the LaTex \prgname macro
        $str =~ s/\\prgname\{(\w+\.*\w*)\}/<em>$1<\/em>/g;
        # replace the LaTex citation
        $str =~ s/\\citep\{(\w+\:\d+\w*)\}/\($1\)/g;
        # remove the italic LaTex formatting
        $str =~ s/\\textit\{(\w+)\}/$1/g;
        $str =~ s/\\textit\{(.+)\}/$1/g;
        # fix the LaTex Anstrom
        $str =~ s/\\AA\{\}/A/g;
        # fix the 2d->2D and 3d->3D issue
        $str =~ s/2d/2D/g;
        $str =~ s/3d/3D/g;
        # remove junk
        $str =~ s/\&lt;comment\/end\&gt;//g;
        $str =~ s/\{//g;
        $str =~ s/\}//g;
    }
    return $str;
}

sub prginstr2f90{
    # input is label and prginstr
    my $key   = shift;
    my %instr = @_;
    if( defined($instr{$key}) ){
        my @lines    = unpack("(A80)*", $instr{$key});
        my $sub_name = 'print_doc_'.$key;
        print "    subroutine $sub_name\n";
        if( $#lines > 0 ){
            foreach my $i ( 0 .. $#lines-1 ){
                print "        write(*,'(A)', advance='no') '$lines[$i]'\n";
            }
        }
        print "        write(*,'(A)') '$lines[-1]'\n";
        print "        stop\n";
        print "    end subroutine $sub_name\n\n";
    }
}

sub extract_comments{
    my $fname = shift;
    open(FORFILE, "<$fname") or die "Cannot open $fname for reading: $!\n";
    my $comments;
    while(my $line=<FORFILE>){
        chomp($line);
        if( $line =~ /^\s*!\s*(.+)$/ ){
            my $line_clean = $1;
            $line_clean =~ s/\s+$//; # removes trailing whitespace
            $line_clean =~ s/^\s+//; # removes leading whitespace
            if( defined($comments) ){
                $comments = $comments.' '.$line_clean;
            }else{
                $comments = ' '.$line_clean;
            }
           
        }
    }
    close(FORFILE);
    return $comments;
}
