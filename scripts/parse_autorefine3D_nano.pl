#!/usr/bin/perl


my @overlap;
my @fracsrch;

while(<>){
    chmod($_);
    open(FOO, "<$_") or die "Cannot open $_: $!\n";
    while($line=<FOO>){
        # ORIENTATION OVERLAP
        if( $line =~ /ORIENTATION OVERLAP/ ){
            my@row  = split(/\s+/,$line);
            push (@overlap, $row[2]);
            print "@overlap\n";
        }
        # % SEARCH SPACE SCANNED
        if( $line =~ /SEARCH SPACE SCANNED/ ){
            my@row  = split(/\s+/,$line);
            push (@fracsrch, $row[2]);
            print "@fracsrch \n";
        }
        #
        # if( $line =~ // ){
        #     my@row  = split(/\s+/,$line);
        #
        # }
        # #
        # if( $line =~ // ){
        #     my@row  = split(/\s+/,$line);
        #
        # }
        # #
        # if( $line =~ /Estimated defocus values(.+)Angstroms/ ){
        #     my@row  = split(/\s+/,$line);
        #
        # }
    close(FOO);
  }
}
