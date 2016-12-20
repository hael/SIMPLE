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
my %prgcmdlineinstr;       # program command line instructions
my @prgclasses;            # program classes
my @prgkeys;               # program keys
my @event_indices;         # critical events that define program classification
my %when2print_prgclasses; # controls timing of printing 
if( scalar(@ARGV) < 1 ){
    die "Need to know which kind of documentation to generate (tex|web|f90)\n";
}
my $doc = $ARGV[0];
my $simple_exec_abspath = "$SIMPLE_PATH/production/simple/simple_exec/simple_exec.f90";

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

# generate labels
my @labels = @prgnames;
foreach my $label (@labels){
    $label =~ s/simple_//;
}

# extract comments
my $comments = extract_comments($simple_exec_abspath);
chomp($comments);

# extract program instructions and program classes
my @comments_split = split('==', $comments);
my $passed = 0;
my $cnt = 0;
foreach my $comment (@comments_split){
    if( $comment =~ /make all programs have the simple_ prefix/ ){
        $passed = 1
    }
    if( $passed == 1 ){
        if( $comment =~ /([A-Z2-3]+[-]*\s*[A-Z2-3]*\s*[A-Z2-3]+\sPROGRAMS)/ ){
            push(@prgclasses,$1);
            push(@event_indices,$cnt)
        }
    }
    if( $comment =~ /\s*\<(.+)\/begin\>((.+))\<.+\/end\>/ ){
        $cnt++;
        my $key = remove_simple($1);
        my $val = remove_simple($2);
        $prginstr{$key} = $val;
        push(@prgkeys,$key);
    }
}
foreach my $i (0 .. $#event_indices){
    $prgclasses[$i] = lc $prgclasses[$i];
    $prgclasses[$i] =~ s/2d/2D/;
    $prgclasses[$i] =~ s/3d/3D/;
    $prgclasses[$i] = ucfirst $prgclasses[$i];
    $when2print_prgclasses{$prgkeys[$event_indices[$i]]} = $prgclasses[$i];
}

# generate command-line dictonary
my $foo     = `$SIMPLE_PATH/bin/simple_exec prg=print_cmd_dict outfile=cmd_dict.txt`;
my $nparams = 0;
open(CMDDICT, "<cmd_dict.txt") or die "Cannot open file: cmd_dict.txt for reading, $!\n";
while(my $keyval=<CMDDICT>){
    chomp($keyval);
    if( $keyval =~ /^(.+)\=((.+))$/ ){
        push(@cmdlinkeys,$1);
        my $key = str2latex(remove_simple($1));
        my $val = str2latex(remove_simple($2));
        $cmdlindict{$key} = $val;
        $nparams++; # number of command line parameters
    }
}

# generate per-program command line dictonaries
foreach my $label (@labels){
    $prgcmdlineinstr{$label} = `$SIMPLE_PATH/bin/simple_exec prg=$label`;
    $prgcmdlineinstr{$label} =~ s/\n$//; # removes trailing newline character
}

# CODE GENERATORS

if( $doc eq 'tex' ){
    # print the command-line dictonary precursor
    print_full_latex_cmdlindict();
    # generate instruction precursors for insertion into the latex manual
    print_latex_subsection($prgclasses[0]);
    foreach my $label (@labels){
        print_latex_instr($label);
        print_latex_usage($label);
        if( defined($when2print_prgclasses{$label}) ){
            print_latex_subsection($when2print_prgclasses{$label});
        }
    }
}elsif( $doc eq 'web' ){
    # generate accordion precursors for insertion into the webpage
    foreach my $label (@labels){
        print_html_instr($label);
    }
}elsif( $doc eq 'f90' ){
    print "module simple_gen_doc\n";
    print "implicit none\n\n";
    print "contains\n\n";
    foreach my $label (@labels){
        prginstr2f90($label);
    }
    print "    subroutine list_all_simple_programs\n";
    my @labels_sorted = sort @labels;
    foreach my $label (@labels_sorted){
        my $prg = 'simple_'.$label;
        print "        write(*,'(A)') '$prg'\n";
    }
    print "        stop\n";
    print "    end subroutine list_all_simple_programs\n\n";
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
    my $label = shift;
    my $label_latex = $label;
    $label_latex =~ s/\_/\\_/g;
    print '\subsubsection{Program: \prgname{'."$label_latex}}\n";
    print '\label{'."$label}\n";
    if( defined($prginstr{$label}) ){
        my $prginstr_latex = str2latex($prginstr{$label});
        print '\prgname{'."$label_latex} ".$prginstr_latex."\n\n";
    }
}

sub print_latex_usage{
    my $label = shift;
    print '\begin{verbatim}', "\n";
    print $prgcmdlineinstr{$label};
    print '\end{verbatim}', "\n\n";
}

sub print_latex_subsection{
    my $heading = shift;
    print '\subsection{'."$heading}\n\n";
}

sub print_html_instr{
    my $label = shift;
    if( defined($prginstr{$label}) ){
        my $instr_html = latex2html(str2latex($prginstr{$label}));
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
    my $key = shift;
    if( defined($prginstr{$key}) ){
        my @lines    = unpack("(A80)*", $prginstr{$key});
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
            # $line_clean =~ s/^\s+//; # removes leading whitespace
            if( defined($comments) ){
                $comments = $comments.' '.$line_clean;
            }else{
                $comments = $line_clean;
            }
        }
    }
    close(FORFILE);
    return $comments;
}

sub remove_simple{
    my $str = shift;
    $str    =~ s/simple_//g; # removes all simple_ prefixes
    return $str;
}
