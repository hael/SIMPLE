#!/usr/bin/perl
use lib './';
use warnings;
use strict;
use Cwd qw(getcwd);
use Env;
use Config;
use Tie::File;
use Term::ANSIColor;
use XML::Parser;
use Data::Dumper;
################################################################################
# Perl script to Parse the ctfpm files after they have been CTFfind            #
#   **** Do not modify this script, this sript automates the compilation ****  #
#   ****                         procedure                               ****  #
#                                                                              #
# Input variables should be modifed in the file:                               #
#                                                                              #
#                          simple_ParseCTFfind_input.pm                        #
#                                                                              #
################################################################################
use simple_ParseCTFfind_input;
################################################################################
# Declare global variables                                                     #
################################################################################
my $sptr    = "/";            # set the name for the separator.
my $all     = "*";            # unix commands for all files eg: file*
my $ext_txt = ".txt";         # txt file name extension: file*.txt
my $ext_xml = ".xml";         # xml file name extension: file*.xml
my $graph   = "/usr/bin/graph -T X --bg-color black --frame-color white";
my $execdir = getcwd();
my $f1;
my $f2;
my $file;
my $file_box;
my $file_ctf;
my $file_def;
my $xml_file;                 #XML file to be parsed
my $xml_dir;                  #XML dir to be parsed
my $nrows = 0;                #number of lines in the file target.txt
my $irows;
my $message;
my @ctf_values;
my @file_array;
my @xml_file_array;
my @xml_dir_array;
my @targetCTFfind_pspec_file; #array for the CTFfind_pspec_files
my @nparticles;               #array for the nparticles per boxed files
my @boxed_file;               #array of files names of the bnoxed files
my @lines;                    #Array for the XML parsing files
#XML parsed values
my $x_pixSze;                 #The parsed x PixelSize
my $y_pixSze;                 #The parsed y PixelSize
my $TemMag;                   #Nominal Magnification
my $AccVlt;                   #Acceleration Voltage (KeV)
#The file names
my $target_CTFfind_pspec_file = "target_CTFFind_pspec.txt";
my $target_Boxfile_file = "target_Boxfile.txt";
my $dftb_stk_file = "Deftabs_stack.asc";
my $stat_all_file = "dfxdfy_all_stat.asc";
my $stat_box_file = "dfxdfy_box_stat.asc";
#File hadnlers
my $ffoo;
my $fboo;
my $fhoo;
my $fqoo;
my $fdoo;
################################################################################
# Start of the execution commands                                              #
################################################################################
print qq[-------------------------------------------------------------------\n];
print qq[   Perl script to strip and analyse mrc files from Krios-Titan     \n];
print qq[           By Frederic D.R. Bonnet date: 11th May. 2015.           \n];
print qq[-------------------------------------------------------------------\n];
print qq[Our working path is       : $SIMPLE_DATA_PATH\n];
print qq[The execution dir         : $execdir\n];
print qq[mrc files to be stripped  : $NCONFIG\n];
if ( $SIMPLE_DATA_PATH ne $execdir ) {
    print qq[\n];
    print qq[-----Our working path is not the same as The execution dir-----\n];
    print qq[---------Check that the paths are the same in input module-----\n];
    print qq[-----------------The script is terminated----------------------\n];
    print qq[\n];
    die;
}
################################################################################
# Setting up the target file                                                   #
################################################################################
$f1 = $target_CTFfind_pspec_file;  # setting up the target file name

if ($FILE_HANDLER =~ /distribute/ || $FILE_HANDLER =~ /donothing/ ) {
    generate_Target_CTFfind_pspec_list($f1);
}
################################################################################
# Extract the files names from the targetfile                                  #
################################################################################
get_target_CTFfind_pspec_file($f1);
################################################################################
# Extract the files names from the targetfile                                  #
################################################################################
open($fhoo,'>', $stat_all_file) or die "cannot open input file $stat_all_file";

#creatign the target_Boxfile.txt
create_boxfile_target();
#getting the number of particles from the target_Boxfile.txt file
get_nboxed_2_defTab_boxfile($target_Boxfile_file);

#parsing the entire data set and printing the stats from all files
for ($irows = 0 ; $irows < ($#targetCTFfind_pspec_file) ; $irows++ ) {
    $file = "$execdir/CTFfind_txt/$targetCTFfind_pspec_file[$irows]";
    @file_array = split('.txt', $file);
    #print qq[0:$file_array[0] 1:$file_array[1]\n];
    $file = "$execdir/CTFfind_txt/$file_array[1]$ext_txt";
    #parsing the pspecfiles for the entire data set
    parsing_target_CTFfind_pspec_file($file);

    #parsing the xml files.
    if ($FILE_HANDLER =~ /regroup/) {
	@xml_file_array = split('_frames',$file_array[1]);
	@xml_file_array = split('/',$xml_file_array[0]);
	$xml_file = $xml_file_array[1].$ext_xml;
	parsing_target_xml_file($xml_file);
    } elsif ($FILE_HANDLER =~ /distribute/) {
	#getting the directory from files 
	@xml_dir_array = split('_frames',$file_array[1]);
	@xml_dir_array = split('/',$xml_dir_array[0]);
	@xml_dir_array = split('_',$xml_dir_array[1]);
	$xml_dir = $xml_dir_array[1];
	#constructing files with the directory structure
	@xml_file_array = split('_frames',$file_array[1]);
	@xml_file_array = split('/',$xml_file_array[0]);
	$xml_file = "$execdir$sptr$xml_dir$sptr$xml_file_array[1]$ext_xml";
	parsing_target_xml_file($xml_file);
    }
}
close $fhoo;

#parsing the entire data set and printing the stats from all files
# and creating the defteb files
open($fdoo,'>', $dftb_stk_file) or die "cannot open input file $dftb_stk_file";
open($fhoo,'>', $stat_box_file) or die "cannot open input file $stat_box_file";
system("mkdir ./Deftabs");
for ($irows = 0 ; $irows < ($#nparticles) ; $irows++ ) {
    $file_box = "$boxed_file[$irows]";
    print qq[$file_box\n];

    @file_array = split('.box', $file_box);
    @file_array = split('Boxfile', $file_array[0]);
    @file_array = split('/', $file_array[1]);
    
    print qq[0:$file_array[0] 1:$file_array[1]\n];

    $file_ctf = "$execdir/CTFfind_txt/$file_array[1]"."_pspec.txt";
    #system("cat $file_ctf");
    parsing_target_CTFfind_pspec_file($file_ctf);

    $file_def = "$execdir/Deftabs/$file_array[1]"."_def.txt";

    write_nboxed_2_Deftabs_boxfile($nparticles[$irows],$file_def);
}

close $fhoo;
close $fdoo;

################################################################################
#                           Subroutines                                        #
################################################################################
################################################################################
#subroutine to extract the number of particles in each box files               #
#Frederic D.R. Bonnet, Date: 3rd of Aug. 2015.                                 #
################################################################################
sub write_nboxed_2_Deftabs_boxfile
{
    (my $nptc,$f1)=@_;
    my $iptc;
    my $dfx;
    my $dfy;
    my $angast;
    my $ctfres;
    print qq[~~~~~~~~~~~~~~~~~~Writting in the Deftabs dir~~~~~~~~~~~~~~~~~~\n];
    print qq[wrtting the parsed line $nptc from $f1                         \n];

    open($fqoo,'>', $f1) or die "cannot open input file $f1";

    #print qq[$nptc $ctf_values[1] $ctf_values[2] $ctf_values[3]\n];
    for ($iptc = 0 ; $iptc < $nptc ; $iptc++ ) {
	#converting from Angstrom to micron meters
	$dfx = $ctf_values[1] * (0.0001);
	$dfy = $ctf_values[2] * (0.0001);
	$angast=$ctf_values[3];
	$ctfres=$ctf_values[6];
	print $fqoo qq[dfx=$dfx dfy=$dfy angast=$angast ctfres=$ctfres\n];
	print $fdoo qq[dfx=$dfx dfy=$dfy angast=$angast ctfres=$ctfres\n];
    }
    close $fqoo;

    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
}
################################################################################
#subroutine to extract the number of particles in each box files               #
#Frederic D.R. Bonnet, Date: 3rd of Aug. 2015.                                 #
################################################################################
sub get_nboxed_2_defTab_boxfile
{
    ($f1)=@_;
    my @Fld;
    print color('bold magenta');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    print qq[Getting the various entries from $f1.                          \n];
    print color ('reset');
    open($ffoo,'<', $f1) or die "cannot open input file $f1";
    my $end = 0;
    my $irows = 0;
    while ( $end == 0 )
    {
	$_= <$ffoo>;
	chop $_[0];
	@Fld = split(' ', $_);
	$nparticles[$irows] = $Fld[0];
	$boxed_file[$irows] = $Fld[3];
#	print qq[$irows nptcl: $nparticles[$irows] file: $boxed_file[$irows]];
	if ( $Fld[3] eq 'total' ){$end = 1 ;}
	if ( $end == 0 ) {$irows = $irows + 1;}
    }

    $nrows = $irows;

    close $ffoo;
    print color('bold magenta');       print color('bold green');
    print qq[$f1 ];                    print color('bold magenta');
    print qq[reports nrows: ];         print color('bold green');
    print qq[$nrows ];                 print color('bold magenta');
    print qq[with ];                   print color('bold green');
    print qq[$nrows ];                 print color('bold magenta');
    print qq[files to be analysed.\n];
    print qq[And a ];  \               print color('bold green');
    print qq[$boxed_file[$nrows] ];    print color('bold magenta');
    print qq[of ];                     print color('bold green');
    print qq[$nparticles[$nrows] ];    print color('bold magenta');
    print qq[boxed particles\n];
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    print color ('reset');
    return($nrows,@nparticles,@boxed_file);

}
################################################################################
#subroutine to extract the number of particles in each box files               #
#Frederic D.R. Bonnet, Date: 3rd of Aug. 2015.                                 #
################################################################################
sub create_boxfile_target
{
    system("wc Boxfile/*.box > $target_Boxfile_file");
}
################################################################################
#subroutine to parse the lines from the old pspec files.                       #
#Frederic D.R. Bonnet, Date: 31th of Jul. 2015.                                #
################################################################################
sub parsing_target_CTFfind_pspec_file 
{
    ($f2)=@_;
    my $ival;

    print color('bold yellow');
    print qq[* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n];
    print color('bold cyan');
    print qq[Processing file:\n];
    print qq[$f2\n];

    open($fboo,'<', $f2) or die "cannot open input file $f2";
    my @lines = <$fboo>;
    close $fboo;
    print color('bold green');
    print qq[$lines[$#lines]];
    print color ('reset');
    @ctf_values = split(' ',$lines[$#lines]);

    for ($ival = 0 ; $ival < ($#ctf_values-1) ; $ival++ ) {
	print $fhoo qq[$ctf_values[$ival]  ];
    }
    print $fhoo qq[$ctf_values[$#ctf_values]\n];

    print color('bold yellow');
    print qq[* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n];
    print color ('reset');
    return (@ctf_values);
}
################################################################################
#subroutine to parse the lines from the old pspec files.                       #
#Frederic D.R. Bonnet, Date: 3rd of Aug. 2015.                                 #
################################################################################
sub parsing_target_xml_file
{
    ($f2)=@_;

    print color('bold yellow');
    print qq[* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n];
    print color('bold cyan');
    print qq[Processing file:\n];
    print qq[$f2\n];

    #TODO: method to parse the xml files

    open($fboo,'<', $f2) or die "cannot open input file $f2";
    @lines = <$fboo>;
    close $fboo;
    print color('magenta');
    #print qq[$lines[$#lines]\n];
    print color ('reset');

    print qq[~~~~~~~~~~~~~~~~~~~~~~Parsing XML file~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    parse_PixelSize(@lines);
    parse_TemMagnification(@lines);
    parse_AccelerationVoltage(@lines);
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];

    print color('bold yellow');
    print qq[* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n];
    print color ('reset');
}
#The parsers
sub parse_AccelerationVoltage {
    (@lines)=@_;
    #doing the parsing, creating our parser object
    my $psr = new XML::Parser(
	Handlers => {
	    Start => \&hdl_start_AccVlt, End => \&hdl_end_AccVlt,
	    Char  => \&hdl_char,     Default => \&hdl_def, });
    $psr->parse($lines[$#lines]);
}
sub parse_TemMagnification {
    (@lines)=@_;
    #doing the parsing, creating our parser object
    my $psr = new XML::Parser(
	Handlers => {
	    Start => \&hdl_start_TemMag, End => \&hdl_end_TemMag,
	    Char  => \&hdl_char,     Default => \&hdl_def, });
    $psr->parse($lines[$#lines]);
}
sub parse_PixelSize {
    (@lines)=@_;
    #doing the parsing, creating our parser object
    my $psr = new XML::Parser(
	Handlers => {
	    Start => \&hdl_start_pixSze, End => \&hdl_end_pixSze,
	    Char  => \&hdl_char,     Default => \&hdl_def, });
    $psr->parse($lines[$#lines]);
}
#The Handlers
sub hdl_start_AccVlt {
    my ($p, $elt, %atts) = @_;
    return unless $elt eq 'AccelerationVoltage';
    $atts{'_str'} = '';
    $message = \%atts;
}
sub hdl_end_AccVlt {
    my ($p, $elt) = @_;
    format_message_AccVlt($message) if $elt eq 'AccelerationVoltage' && $message && $message->{'_str'} =~ /\S/;
}
sub format_message_AccVlt {
    my $atts = shift;
    $atts->{'_str'} =~ s/\n//g;
    $AccVlt = $atts->{'_str'} / 1000;
    print "Acceleration Voltage (KeV): $AccVlt\n";
}
sub hdl_start_TemMag {
    my ($p, $elt, %atts) = @_;
    return unless $elt eq 'TemMagnification';
    $atts{'_str'} = '';
    $message = \%atts;
}
sub hdl_end_TemMag {
    my ($p, $elt) = @_;
    format_message_TemMag($message) if $elt eq 'TemMagnification' && $message && $message->{'_str'} =~ /\S/;
}
sub format_message_TemMag {
    my $atts = shift;
    $atts->{'_str'} =~ s/\n//g;
    $TemMag = $atts->{'_str'};
    print "Nominal Magnification: $TemMag\n";
}
sub hdl_start_pixSze {
    my ($p, $elt, %atts) = @_;
    return unless $elt eq 'pixelSize';
    $atts{'_str'} = '';
    $message = \%atts;
}
sub hdl_end_pixSze {
    my ($p, $elt) = @_;
    format_message_pixSze($message) if $elt eq 'pixelSize' && $message && $message->{'_str'} =~ /\S/;
}
sub format_message_pixSze {
    my $atts = shift;
    $atts->{'_str'} =~ s/\n//g;

    my @Fld = split('1m',$atts->{'_str'});
    $x_pixSze = $Fld[0];
    $y_pixSze = $Fld[1];
    print "PixelSize(Angstroms): x=$x_pixSze y=$y_pixSze\n";
}
#common handlers
sub hdl_char {
    my ($p, $str) = @_;
    $message->{'_str'} .= $str;
}
sub hdl_def{}
################################################################################
#subroutine to generate the target file for *CTFFind_pspec.txt                 #
#Frederic D.R. Bonnet, Date: 31st of Jul. 2015.                                #
################################################################################
sub generate_Target_CTFfind_pspec_list
{
    ($f1)=@_;
    print qq[Generating the target file for *CTFFind_pspec.txt: $f1\n];
    system("ls CTFfind_txt/FoilHole_1*_pspec.txt > $f1");
}
################################################################################
#subroutine to extract the prices from the old TickData files.                 #
#Frederic D.R. Bonnet, Date: 31st of Jul. 2015.                                #
################################################################################
sub get_target_CTFfind_pspec_file
{
    ($f1)=@_;
    my @Fld;
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    print qq[getting the various entries from the file name.                \n];
    print qq[We are now printing out the details from  the file: $f1        \n];
    open($ffoo,'<', $f1) or die "cannot open input file $f1";
    my $end = 0;
    my $irows = 0;
    while ( $end == 0 )
    {
	$_= <$ffoo>;
	chop $_[0];
	@Fld = split('/', $_);
	$targetCTFfind_pspec_file[$irows] = $Fld[1];
	#print qq[$irows $targetCTFfind_pspec_file[$irows]];
	if ( $Fld[0] eq '' ){$end = 1 ;}
	if ( $end == 0 ) {$irows = $irows + 1;}
    }

    $nrows = $irows;

    close $ffoo;
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];
    return($nrows,@targetCTFfind_pspec_file);
}
