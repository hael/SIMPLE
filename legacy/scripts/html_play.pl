#!/usr/bin/perl
use warnings;
use strict;

my $instr = 'is a great program';

my $cmdlin = 'SIMPLE_AUTOMASK vol1=<invol.ext> [vol2=<invol2.ext> etc.] smpd=<sampling distance(in A)> [mw=<molecular weight(in kD)>] [amsklp=<low-pass limit(in A){20}>] [nthr=<nr of OpenMP threads{1}>]**less commonly used**[edge=<edge size for softening molecular envelope(in pixels){3}>] [dens=<density(e.g. 9.368 Da/A3 4 gold clusters){0.}>] [nvox=<nr of voxels in mask{0}>] [bin=<yes|no{no}>]';

my $comments = 'no comments';

instr2html($instr, $cmdlin, $comments);

sub instr2html{
    my $instr    = shift;
    my $cmdlin   = shift;
    my $comments = shift;
    my $prgname;
    print qq[<div class="AccordionPanel">\n];
    if( $cmdlin =~ /^(SIMPLE_\w+)\s/ ){
        $prgname = lc $1;
        print qq[  <div class="AccordionPanelTab">$prgname</div>\n];
    }
    print qq[  <div class="AccordionPanelContent">\n];
    print qq[    <p><strong>$prgname</strong> $instr</p>\n];
    print qq[    <p align="left"><span class="h2_manual">Usage:<br>\n];
    # replace < (less than)
    $cmdlin =~ s/\</\&lt;/g;
    # replace > (greater than)
    $cmdlin =~ s/\>/\&gt;/g;
    my $replacement = '    </span>  &gt;&gt;&gt;<strong>'.$prgname.' </strong>';
    $cmdlin =~ s/^SIMPLE_\w+\s/$replacement/;
    while( $cmdlin =~ /(\w+\d*)=/ ){
        my $new_val = '<em>'.$1.'</em>';
        $cmdlin =~ s/$1/$new_val/;
    }
    if( $cmdlin =~ /\*\*less commonly used\*\*/ ){
        my @parts = split(/\*\*less commonly used\*\*/, $cmdlin);
        print qq[$parts[0]<br>\n];
        print qq[    <strong>** less commonly used**</strong> $parts[1]</p>\n];
    }
    
    if( defined($comments) ){
        print qq[<p><span class="h2_manual">Comments:<br>\n];
        print qq[</span> $comments</p>\n];
    }
    print qq[</div>\n];
    print qq[</div>\n];
}

sub texdescr2html{
    my $str = shift;
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
    $str =~ s/\\texttt\{(.+)\}/<em>$1<\/em>/g;
    # remove the LaTex \prgname macro
    $str =~ s/\\prgname\{(\w+\.*\w*)\}/$1/g;
    # replace the LaTex citation
    $str =~ s/\\citep\{(\w+\:\d+\w*)\}/\($1\)/g;
    # remove the italic LaTex formatting
    $str =~ s/\\textit\{(\w+)\}/$1/g;
    $str =~ s/\\textit\{(.+)\}/$1/g;
    # fix the LaTex Anstrom
    $str =~ s/\\AA\{\}//g;
    # remove junk
    $str =~ s/\&lt;comment\/end\&gt;//g;
    $str =~ s/\{//g;
    $str =~ s/\}//g;
    return $str;
}