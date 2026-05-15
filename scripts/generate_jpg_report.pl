#!/usr/bin/env perl
use strict;
use warnings;

use File::Find qw(find);
use File::Basename qw(basename);
use File::Spec;
use Cwd qw(abs_path getcwd);
use POSIX qw(strftime);

my $input_dir = $ARGV[0];
my $output_html = 'artifact.html';
my $title = 'SIMPLE Test Report';
my $thumb_width = 220;

die usage() unless defined $input_dir && @ARGV == 1;

$input_dir = abs_path($input_dir)
    or die "ERROR: Cannot resolve input directory: $input_dir\n";

die "ERROR: Input path is not a directory: $input_dir\n"
    unless -d $input_dir;

my $output_abs;
if (File::Spec->file_name_is_absolute($output_html)) {
    $output_abs = $output_html;
} else {
    my $cwd = getcwd();
    $output_abs = File::Spec->catfile($cwd, $output_html);
}

my ($vol, $out_dirs, $out_file) = File::Spec->splitpath($output_abs);
my $output_dir = File::Spec->catpath($vol, $out_dirs, '');

if (length $output_dir && !-d $output_dir) {
    my @parts = File::Spec->splitdir($output_dir);
    my $path = File::Spec->rootdir();
    for my $part (@parts) {
        next if $part eq q{};
        $path = File::Spec->catdir($path, $part);
        mkdir $path unless -d $path;
    }
}

my @images;

find(
    {
        wanted => sub {
            return unless -f $_;
            return unless /cavgs/i;
            return unless /\.(?:jpe?g)$/i;
            push @images, abs_path($File::Find::name);
        },
        no_chdir => 1,
    },
    $input_dir
);

@images = sort @images;

my $last_iter = -1;
my @last_iter_images;

for my $img_abs (@images) {
    my $name = basename($img_abs);
    my $iter;

    if ($name =~ /(?:iter|iteration)[_\-\s]*([0-9]+)/i) {
        $iter = $1;
    } elsif ($name =~ /cavgs[^0-9]*([0-9]+)/i) {
        $iter = $1;
    }

    next unless defined $iter;

    if ($iter > $last_iter) {
        $last_iter = $iter;
        @last_iter_images = ($img_abs);
    } elsif ($iter == $last_iter) {
        push @last_iter_images, $img_abs;
    }
}

if ($last_iter >= 0) {
    @images = sort @last_iter_images;
}

open(my $out, '>', $output_abs) or die "ERROR: Cannot write $output_abs: $!\n";

my $time_str = strftime('%Y-%m-%d %H:%M:%S', localtime());
my $status = @images ? 'Test Passed' : 'No matching cavgs JPG files found';

print {$out} "<!doctype html>\n";
print {$out} "<html lang=\"en\">\n";
print {$out} "<head>\n";
print {$out} "  <meta charset=\"utf-8\">\n";
print {$out} "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n";
print {$out} "  <title>$title</title>\n";
print {$out} "  <style>\n";
print {$out} "    body { font-family: Arial, sans-serif; margin: 0; background: #f6f8fa; color: #1f2328; }\n";
print {$out} "    .wrap { max-width: 1200px; margin: 0 auto; padding: 24px; }\n";
print {$out} "    .card { background: #fff; border: 1px solid #d0d7de; border-radius: 12px; padding: 16px 20px; margin-bottom: 18px; }\n";
print {$out} "    h1 { margin: 0 0 10px; font-size: 28px; }\n";
print {$out} "    .meta { color: #57606a; font-size: 14px; }\n";
print {$out} "    .status { display: inline-block; margin-top: 10px; padding: 6px 10px; border-radius: 999px; background: #dafbe1; color: #116329; font-weight: 600; }\n";
print {$out} "    .grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(240px, 1fr)); gap: 14px; }\n";
print {$out} "    .tile { background: #fff; border: 1px solid #d0d7de; border-radius: 10px; overflow: hidden; }\n";
print {$out} "    .tile img { width: 100%; height: auto; display: block; }\n";
print {$out} "    .tile .cap { padding: 8px 10px; font-size: 12px; color: #57606a; word-break: break-all; }\n";
print {$out} "    .empty { padding: 20px; text-align: center; color: #57606a; background: #fff; border: 1px dashed #d0d7de; border-radius: 12px; }\n";
print {$out} "  </style>\n";
print {$out} "</head>\n";
print {$out} "<body>\n";
print {$out} "  <div class=\"wrap\">\n";
print {$out} "    <div class=\"card\">\n";
print {$out} "      <h1>$title</h1>\n";
print {$out} "      <div class=\"meta\">Generated: $time_str</div>\n";
print {$out} "      <div class=\"meta\">Input directory: $input_dir</div>\n";
if ($last_iter >= 0) {
    print {$out} "      <div class=\"meta\">Selected iteration: $last_iter</div>\n";
} else {
    print {$out} "      <div class=\"meta\">Selected iteration: not detected (showing all matching cavgs images)</div>\n";
}
print {$out} "      <div class=\"meta\">Images found: " . scalar(@images) . "</div>\n";
print {$out} "      <div class=\"status\">$status</div>\n";
print {$out} "    </div>\n";

if (@images) {
    print {$out} "    <div class=\"grid\">\n";
    for my $img_abs (@images) {
        my $rel = File::Spec->abs2rel($img_abs, $output_dir || abs_path(File::Spec->curdir()));
        my $name = basename($img_abs);
        print {$out} qq{      <div class="tile">\n};
        print {$out} qq{        <a href="$rel" target="_blank" rel="noopener noreferrer">\n};
        print {$out} qq{          <img src="$rel" alt="$name" width="$thumb_width">\n};
        print {$out} qq{        </a>\n};
        print {$out} qq{        <div class="cap">$name</div>\n};
        print {$out} qq{      </div>\n};
    }
    print {$out} "    </div>\n";
} else {
    print {$out} "    <div class=\"empty\">No JPG/JPEG files containing 'cavgs' were found in the provided directory.</div>\n";
}

print {$out} "  </div>\n";
print {$out} "</body>\n";
print {$out} "</html>\n";

close($out);

print "Wrote report: $output_abs\n";
exit 0;

sub usage {
    return <<'USAGE';
Usage:
    generate_jpg_report.pl <directory>

Notes:
    - Only one argument is supported.
    - Output file is generated as ./artifact.html
USAGE
}
