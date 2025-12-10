#!/usr/bin/perl
# Perl script to delete unused variables in the source code flagged in debug compilation mode
use strict;
use warnings;
sub trim { $_[0] =~ s/^\s+|\s+$//gr }

# compile code in debug mode
print(">>> CLONE REPO & COMPILE SIMPLE IN DEBUG MODE\n");
system("rm -f build_log\* unused_log unused_after_remove_log \;git clone git\@github.com:hael/SIMPLE.git\; cd SIMPLE\; ./compile_debug.sh >../build_log 2>&1\; cd ..\;grep \"Unused variable\" -B 4 build_log > unused_log"); 
system("pwd"); 

print(">>> REMOVE UNUSED VARIABLES\n");
open( my $build_filehandle, "<", "unused_log" ) || die; # opening log compilation output file
my $file;               
my $line_to_edit;       
my $col_to_edit;        
my $word_to_remove;
my @files;
my @fout;
my @line_nums;
my @unused_variables;
my $line;

my $line_number = 0;
while($line = <$build_filehandle>)
{
    $line_number++;
    if($line =~ m/([^:]*):(\d+)\:(\d+)/)
    {
        $file = $1;
        $line_to_edit = $2;
        $col_to_edit = $3;
    }
    elsif($line =~ /Unused\s+variable\s+\‘([a-zA-Z0-9_]+)\’/) {
        $word_to_remove = $1;
        push(@files, $file);
        push(@line_nums, $line_to_edit);
        push(@unused_variables, $word_to_remove);
    }
}
close($build_filehandle);

# sort [file, line_num, word_to_remove] by files then line_nums then word_to_remove
my @merged_columns;
for(my $i = 0; $i < scalar(@files); $i++) # loop over all the replacements to be made
{
    push(@merged_columns, [$files[$i],$line_nums[$i],$unused_variables[$i]]);
}
# if sort is stable, sort by line_nums, then files
sub by_file_then_line
{
    $a->[0] cmp $b->[0] || # by file
    $a->[1] <=> $b->[1];
}

my @sorted_by_line_nums = sort by_file_then_line @merged_columns;
undef @merged_columns;

for(my $i = 0; $i < scalar(@files); $i++) # loop over all the replacements to be made
{
    $files[$i] = $sorted_by_line_nums[$i][0];
    $fout[$i] = $files[$i] =~ s/\.[^.]+$/_replace.f90"/r; 
    $line_nums[$i] = $sorted_by_line_nums[$i][1];
    $unused_variables[$i] = $sorted_by_line_nums[$i][2];
}

for(my $i = 0; $i < scalar(@files); $i++) # loop over all the replacements to be made
{
    open my $IN,  "+<", $files[$i] or die "Cannot open $files[$i]: $!";
    open my $OUT,  ">", $fout[$i] or die "Cannot open $fout[$i]: $!";
    print("file                     : $files[$i]\n");
    print("line number              : $line_nums[$i]\n");
    print("remove unused variable   : ",trim($unused_variables[$i]), "\n");
    my $line_num = 0;
    while (my $line = <$IN>) {
          $line_num++;
          if ($line_num == $line_nums[$i]) {
              print "before deletion          : ",trim($line), "\n";
              # replace 'old' with 'new' on these lines only
              $line =~ s/,(\S)/, $1/g; # insert space after comma 
              # case followed by stuff in parens
              $line =~ s/ ${unused_variables[$i]}\([^\)]*\) [^,!]*//i;
              # case followed by stuff in parens
              $line =~ s/ ${unused_variables[$i]}\([^\)]*\),?//i;
              # case followed by equal and variable 
              $line =~ s/ ${unused_variables[$i]}\s*\=.*?(?:,|$|!)//i;
              # case followed by comma
              $line =~ s/ ${unused_variables[$i]},//i;
              # case end of line
              $line =~ s/ ${unused_variables[$i]}\s*$//i;
              # remove trailing commas
              $line =~ s/,\s*$/\n/;
              # collapse double commas
              $line =~ s/,\s*,/,/;
              # the only variable defined then remove that line
              $line =~ s/^.*::\s*$/\n/;
              if ($line =~ /::\s*!/) {$line="\n"}; # special case for line with a comment
              print "after deletion           : ",trim($line), "\n\n";
          }
          print{$OUT} $line;
    }
    close $IN;
    close $OUT; 
    unlink($files[$i]) or die "Could not delete file: $!";
    rename($fout[$i], $files[$i]) or die "rename failed: $!";
}

my @merged_columns1;
for(my $i = 0; $i < scalar(@files); $i++) # loop over all the replacements to be made
{
    push(@merged_columns1, [$files[$i],$line_nums[$i],$unused_variables[$i]]);
}
my %lines_by_file;
for my $row (@merged_columns1) {
    my ($file, $line_nums, $var) = @$row;
    push @{ $lines_by_file{$file} }, $line_nums;
}
for my $file (keys %lines_by_file) {
    my %seen;
    my @unique = grep { !$seen{$_}++ } @{ $lines_by_file{$file} };
    # store unique values back into the hash
    $lines_by_file{$file} = \@unique;
}
# Print
print(">>> DELETE CREATED BLANK LINES\n");
for my $file (sort keys %lines_by_file) {
    my $outf = $file =~ s/\.[^.]+$/_replace.f90"/r; 
    open my $IN,  "+<", $file or die "Cannot open $file: $!";
    open my $OUT,  ">", $outf or die "Cannot open $outf: $!";

    my @line_nums =();
    for my $line (@{ $lines_by_file{$file} }) {
        push @line_nums, $line;
    }

    my $line_num = 0;
    while (my $line = <$IN>) {
         $line_num++;
         if (grep { $_ == $line_num }  @line_nums) {
         #if ($line_num ~~ @line_nums) { ! introduced in Perl 5.10
             if ($line eq "\n") { 
                 print("file       : $file\n");
                 print("line number: $line_num\n");
                 print "line       : ",trim($line), "\n\n";
             }
             else {
                 print{$OUT} $line;
             }
         }
         else {
             print{$OUT} $line;
         }
    }
    close $IN;
    close $OUT; 
    unlink($file) or die "Could not delete file: $!";
    rename($outf, $file) or die "rename failed: $!";
}

# compile in debug mode to check that everything is fine
print(">>> COMPILE NEW VERSION OF SIMPLE IN DEBUG MODE\n");
system("cd SIMPLE\; ./compile_debug.sh > ../build_log_1 2>&1") == 0 
    or die ">>> FAILED DELETION OF UNUSED VARIABLES - SIMPLE COMPILATION ERROR: $?";;
system("grep \"Unused variable\" build_log_1 > unused_after_removal_log");
if (-z "unused_after_removal_log") {
    print ">>> SUCCESFUL DELETION OF UNUSED VARIABLES\n";
}
else {
    print ">>> FAILED DELETION OF UNUSED VARIABLES\n";
}
