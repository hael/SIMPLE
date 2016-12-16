#!/usr/bin/perl
use File::Temp qw/ tempfile /;
use Tie::File;
if( scalar(@ARGV) < 1 ){
    die "Need filetable as input (integrated movies or micrograps) and remember to set parameters in script\n";
}
# INPUT PARAMETERS
$Pix_Sze = 1.06;  # sampling distance (a/pix)
$Acc_Vlt = 300.0; # acceleration voltage (in thousands of V)
$Sph_Abe = 2.7;   # spherical aberration constant (in mm)
$Amp_Con = 0.07;  # fraction of amplitude contrast (assuming protein, as measured by Unwin/Henderson)
$deftab  = 'ctffind_output.txt'; # file with output parameters
unlink $deftab;
unlink "Diag.mrc";
# CONSTANTS
$Sze_Pwr_Spc = 1024;      # 512x512 2D FT
$Min_Res     = 30.0;     # Angstroms
$Max_Res     = 5.0;      # Angstroms
$Min_Dfc     = 5000.0;   # Angstroms
$Max_Dfc     = 70000.0;  # Angstroms
$Dfc_Srh     = 500.0;    # Angstroms
$Exp_Atg     = 100.0;    # Angstroms
$Add_Phs     = 'no';     # Volta phase plate
# PARSE FILES IN PARTITION
@files;
tie @files, 'Tie::File', $ARGV[0];
@dfx_arr;
@dfy_arr;
@angast_arr;
@ctfres_arr;
$cnt = 0;
open(DEFTAB, ">>$deftab") or die 'Cannot open file for writing: $!';
foreach $file (@files){
	$cnt++;
	$diagnostic_file = $file;
	$diagnostic_file = s/\.mrc//;
	$diagnostic_file = $diagnostic_file.'Diag.mrc';
	# MAKE A TEMPORARY FILE FOR THE CTFFIND PARAMETERS
	($ctffind_handle, $ctffind_param_file) = tempfile();
	# WRITE THE NECESSARY PARAMETERS
	print $ctffind_handle "$file\n";
	print $ctffind_handle "$diagnostic_file\n";
	print $ctffind_handle "$Pix_Sze\n";
	print $ctffind_handle "$Acc_Vlt\n";
	print $ctffind_handle "$Sph_Abe\n";
	print $ctffind_handle "$Amp_Con\n";
	print $ctffind_handle "$Sze_Pwr_Spc\n";
	print $ctffind_handle "$Min_Res\n";
	print $ctffind_handle "$Max_Res\n";
	print $ctffind_handle "$Min_Dfc\n";
	print $ctffind_handle "$Max_Dfc\n";
	print $ctffind_handle "$Dfc_Srh\n";
	print $ctffind_handle "no\n";
        print $ctffind_handle "no\n";
        print $ctffind_handle "yes\n";
	print $ctffind_handle "$Exp_Atg\n";
	print $ctffind_handle "$Add_Phs\n";
	print $ctffind_handle "no\n";
	close($ctffind_handle);
	$ctffind_output = `cat $ctffind_param_file | ctffind`;
	print $ctffind_output;
	if( $ctffind_output =~ /Estimated defocus values(.+)Angstroms/){
		@dfxy = split(',', $1);
		$dfxy[0] =~ s/://;
		$dfxy[0] =~ s/\s+//g;
		$dfxy[1] =~ s/\s+//g;
		$dfx_arr[$cnt] = $dfxy[0]/10000.; # Anstroms to microns
		$dfy_arr[$cnt] = $dfxy[1]/10000.; # Anstroms to microns
	}
	if( $ctffind_output =~ /Estimated azimuth of astigmatism(.+)degrees/){
		$angast = $1;
		$angast =~ s/://;
		$angast =~ s/\s+//g;
		$angast_arr[$cnt] = $angast;
	}
	if( $ctffind_output =~ /Thon rings with good fit up to(.+)Angstroms/){
		$ctfres = $1;
		$ctfres =~ s/://;
		$ctfres =~ s/\s+//g;
		$ctfres_arr[$cnt] = $ctfres;
	}
	print DEFTAB "dfx=$dfx_arr[$cnt] dfy=$dfy_arr[$cnt] angast=$angast_arr[$cnt] ctfres=$ctfres_arr[$cnt]\n";
	$frac = sprintf("%.1f", 100.0*($cnt/scalar(@files)));
	print "$frac percent of the data processed\n";
}
close(DEFTAB);
