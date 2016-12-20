#!/usr/bin/perl
use warnings;
use strict;
#Parsing the xml files for the target_files only need the first 1 or 2 
parse_XML_file();


sub parse_XML_file {
    # only need the kfirst one or two files in the list as the parsed param are
    # are fixed for the entire measurment

    for ($irows = 0; $irows < ($#targetFrames-($#targetFrames-2)); $irows++ ) {
	#constructing files with the directory structure
	@xml_file_array = split("$ext_mrc",$targetFrames[$irows]);
	@xml_file_array = split("$ext_frames",$xml_file_array[0]);

	#parsing the xml files.
	if ($has_been_organised == 0) {
	    $xml_file = "$xml_file_array[0]$ext_xml";	    
	} elsif ($has_been_organised == 1) {
	    $xml_file = "$dir_xml$sptr$xml_file_array[0]$ext_xml";
	}

	#print qq[$xml_file\n];
	parsing_target_xml_file($xml_file);

    }
    print color('bold cyan');
    print qq[The parsed PixelSize x: $x_pixSze, y: $y_pixSze, ];
    print qq[Nominal Magnification: $TemMag and ];
    print qq[Acceleration Voltage (KeV): $AccVlt\n];
    print color('reset');
    print qq[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n];

}