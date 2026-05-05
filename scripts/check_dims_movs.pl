#!/usr/bin/env perl
# Check dimensions of movie files from a filetab or directory (MRC/EER/TIFF)
use strict;
use warnings;
use Fcntl qw(SEEK_SET);
use Cwd 'abs_path';

if (@ARGV > 1) {
    die "Usage: $0 [filetab.txt|directory]\n";
}

my $input = @ARGV ? $ARGV[0] : '.';
my @files;

if (-d $input) {
    my $base_dir = abs_path($input);
    opendir(my $dh, $input) or die "Cannot open directory $input: $!\n";
    my @entries = readdir($dh);
    closedir($dh);

    my @mrc_files  = sort grep { /\.mrc$/i        && -f "$base_dir/$_" } @entries;
    my @eer_files  = sort grep { /\.eer$/i        && -f "$base_dir/$_" } @entries;
    my @tiff_files = sort grep { /\.(?:tif|tiff)$/i && -f "$base_dir/$_" } @entries;

    my @selected;
    my $kind;
    if (@mrc_files) {
        @selected = @mrc_files;
        $kind = 'MRC';
    } elsif (@eer_files) {
        @selected = @eer_files;
        $kind = 'EER';
    } elsif (@tiff_files) {
        @selected = @tiff_files;
        $kind = 'TIFF';
    } else {
        die "No .mrc, .eer, .tif, or .tiff files found in $base_dir\n";
    }

    @files = map { abs_path("$base_dir/$_") } @selected;
    print "Searched directory: $base_dir\n";
    print "Using $kind entries: " . scalar(@files) . " file(s)\n";
} else {
    my $filetab = $input;
    open my $ft, '<', $filetab or die "Cannot open $filetab: $!\n";
    while (my $line = <$ft>) {
        chomp $line;
        $line =~ s/^\s+|\s+$//g;
        next if $line eq '';
        next if $line =~ /^\s*#/;
        my ($f) = split /\s+/, $line;
        push @files, $f;
    }
    close $ft;
    print "Read filetab: $filetab\n";
}

my %dims;        # key -> [ [idx, file], ... ]
my @unreadable;  # [idx, file, reason]

my $idx = 0;
for my $f (@files) {
    $idx++;

    my ($nx, $ny, $nz);
    my $ok = eval {
        ($nx, $ny, $nz) = read_dims($f);
        1;
    };

    if (!$ok) {
        my $err = $@ || 'unknown error';
        chomp $err;
        push @unreadable, [$idx, $f, $err];
        next;
    }

    my $key = join(',', $nx, $ny, $nz);
    push @{ $dims{$key} }, [$idx, $f];
}

print "\n=== Dimension groups ===\n";
for my $k (sort { scalar(@{$dims{$b}}) <=> scalar(@{$dims{$a}}) } keys %dims) {
    print "($k) ", scalar(@{$dims{$k}}), "\n";
}

if (!keys %dims) {
    if (@unreadable) {
        print "\n=== Unreadable entries ===\n";
        for my $u (@unreadable) {
            print "  entry #$u->[0]: $u->[1]\n";
            print "    reason: $u->[2]\n";
        }
    }
    die "No readable entries found.\n";
}

my ($major_key) = sort { scalar(@{$dims{$b}}) <=> scalar(@{$dims{$a}}) } keys %dims;
my $major_count = scalar @{ $dims{$major_key} };
print "\nMajority dimension: ($major_key) (count=$major_count)\n";

my @diff = grep { $_ ne $major_key } keys %dims;
if (!@diff) {
    print "No different files: all entries match the majority.\n";
} else {
    print "\n=== Different from majority ===\n";
    for my $k (sort { scalar(@{$dims{$b}}) <=> scalar(@{$dims{$a}}) } @diff) {
        my $items = $dims{$k};
        print "Dimension ($k) (count=", scalar(@$items), "):\n";
        for my $it (@$items) {
            print "  entry #$it->[0]: $it->[1]\n";
        }
    }
}

if (@unreadable) {
    print "\n=== Unreadable entries ===\n";
    for my $u (@unreadable) {
        print "  entry #$u->[0]: $u->[1]\n";
        print "    reason: $u->[2]\n";
    }
}

exit 0;

sub read_dims {
    my ($path) = @_;
    if ($path =~ /\.(?:eer|tif|tiff)$/i) {
        return read_tiff_dims($path);
    }
    return read_mrc_dims($path);
}

sub read_mrc_dims {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "cannot open MRC '$path': $!";
    my $buf = '';
    my $n = read($fh, $buf, 12);
    close $fh;
    die "cannot read MRC header from '$path'" if !defined $n || $n < 12;

    # MRC header starts with nx, ny, nz as 32-bit integers (little-endian)
    my ($nx, $ny, $nz) = unpack('l<l<l<', $buf);
    die "invalid MRC dimensions in '$path'" if $nx <= 0 || $ny <= 0 || $nz <= 0;
    return ($nx, $ny, $nz);
}

sub read_tiff_dims {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "cannot open TIFF/EER '$path': $!";

    my $hdr = '';
    my $n = read($fh, $hdr, 16);
    die "cannot read TIFF header from '$path'" if !defined $n || $n < 8;

    my $endian = substr($hdr, 0, 2);
    my $le;
    if ($endian eq 'II') { $le = 1; }
    elsif ($endian eq 'MM') { $le = 0; }
    else { die "invalid TIFF byte order in '$path'"; }

    my $magic = u16(substr($hdr, 2, 2), $le);
    my ($bigtiff, $ifd_off, $inline_bytes);

    if ($magic == 42) {
        $bigtiff = 0;
        $ifd_off = u32(substr($hdr, 4, 4), $le);
        $inline_bytes = 4;
    } elsif ($magic == 43) {
        die "truncated BigTIFF header in '$path'" if $n < 16;
        my $offsize = u16(substr($hdr, 4, 2), $le);
        my $zero    = u16(substr($hdr, 6, 2), $le);
        die "unsupported BigTIFF offset size in '$path'" if $offsize != 8 || $zero != 0;
        $bigtiff = 1;
        $ifd_off = u64(substr($hdr, 8, 8), $le);
        $inline_bytes = 8;
    } else {
        die "not a TIFF/EER file (magic=$magic) in '$path'";
    }

    die "invalid first IFD offset in '$path'" if !$ifd_off;

    my ($width, $height, $pages) = (undef, undef, 0);

    while ($ifd_off) {
        $pages++;
        seek($fh, $ifd_off, SEEK_SET) or die "seek IFD failed in '$path'";

        my ($n_entries, $entry_size, $next_ptr_size);
        if ($bigtiff) {
            my $buf8 = read_exact($fh, 8, $path);
            $n_entries = u64($buf8, $le);
            $entry_size = 20;
            $next_ptr_size = 8;
        } else {
            my $buf2 = read_exact($fh, 2, $path);
            $n_entries = u16($buf2, $le);
            $entry_size = 12;
            $next_ptr_size = 4;
        }

        if (!defined $width || !defined $height) {
            for (my $i = 0; $i < $n_entries; $i++) {
                my $e = read_exact($fh, $entry_size, $path);
                my ($tag, $type, $count, $valoff);
                if ($bigtiff) {
                    $tag   = u16(substr($e, 0, 2), $le);
                    $type  = u16(substr($e, 2, 2), $le);
                    $count = u64(substr($e, 4, 8), $le);
                    $valoff = substr($e, 12, 8);
                } else {
                    $tag   = u16(substr($e, 0, 2), $le);
                    $type  = u16(substr($e, 2, 2), $le);
                    $count = u32(substr($e, 4, 4), $le);
                    $valoff = substr($e, 8, 4);
                }

                next unless ($tag == 256 || $tag == 257); # ImageWidth / ImageLength
                next unless $count >= 1;

                my $v = read_tiff_value($fh, $type, $count, $valoff, $inline_bytes, $le, $path);
                if ($tag == 256) { $width = $v; }
                if ($tag == 257) { $height = $v; }
            }

            # If width/height still missing, skip remaining entries and continue.
            if (!defined $width || !defined $height) {
                # already consumed all entries
            }
        } else {
            # skip entries quickly
            seek($fh, $n_entries * $entry_size, 1) or die "seek entries failed in '$path'";
        }

        my $nextbuf = read_exact($fh, $next_ptr_size, $path);
        $ifd_off = $bigtiff ? u64($nextbuf, $le) : u32($nextbuf, $le);
    }

    close $fh;

    die "could not read TIFF width/height from '$path'" if !defined($width) || !defined($height);
    die "invalid TIFF dimensions in '$path'" if $width <= 0 || $height <= 0 || $pages <= 0;

    return ($width, $height, $pages);
}

sub read_tiff_value {
    my ($fh, $type, $count, $valoff_raw, $inline_bytes, $le, $path) = @_;

    my %type_size = (
        1  => 1,  # BYTE
        2  => 1,  # ASCII
        3  => 2,  # SHORT
        4  => 4,  # LONG
        5  => 8,  # RATIONAL
        16 => 8,  # LONG8 (BigTIFF)
        17 => 8,  # SLONG8 (BigTIFF)
    );

    my $sz = $type_size{$type} // die "unsupported TIFF tag type=$type in '$path'";
    my $total = $sz * $count;

    my $data;
    if ($total <= $inline_bytes) {
        $data = substr($valoff_raw, 0, $total);
    } else {
        my $off = ($inline_bytes == 8) ? u64($valoff_raw, $le) : u32($valoff_raw, $le);
        seek($fh, $off, SEEK_SET) or die "seek tag value failed in '$path'";
        $data = read_exact($fh, $total, $path);
    }

    # Return first value (count>=1)
    if ($type == 3) { return u16(substr($data, 0, 2), $le); }
    if ($type == 4) { return u32(substr($data, 0, 4), $le); }
    if ($type == 16 || $type == 17) { return u64(substr($data, 0, 8), $le); }
    if ($type == 1 || $type == 2) { return ord(substr($data, 0, 1)); }

    die "unsupported TIFF value conversion for type=$type in '$path'";
}

sub read_exact {
    my ($fh, $len, $path) = @_;
    my $buf = '';
    my $got = read($fh, $buf, $len);
    die "unexpected EOF while reading '$path'" if !defined($got) || $got != $len;
    return $buf;
}

sub u16 {
    my ($b, $le) = @_;
    return unpack($le ? 'v' : 'n', $b);
}

sub u32 {
    my ($b, $le) = @_;
    return unpack($le ? 'V' : 'N', $b);
}

sub u64 {
    my ($b, $le) = @_;
    my ($lo, $hi) = unpack($le ? 'V2' : 'N2', $b);
    return $le ? ($lo + 4294967296 * $hi) : ($hi + 4294967296 * $lo);
}
