// LibTIFF wrapper edited from Unblur
// Copyright 2014 Howard Hughes Medical Institute
// All rights reserved
// Use is subject to Janelia Farm Research Campus Software Copyright 1.1
// license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
// Custom EER related routines adapted from Relion
#include <tiffio.h>

TIFFErrorHandler warn = NULL;

int TIFFRawStripSizer(TIFF *my_tiff, int strip)
{
#ifdef TIFF430
    return TIFFRawStripSize(my_tiff, (uint16_t)strip);
#else
    return TIFFRawStripSize(my_tiff, (uint16)strip);
#endif
}

void EERDecode_7bit(tdata_t bytes, int bitpos, int *p1, int *s1, int *p2, int *s2)
{
	const unsigned int bit_offset_in_first_byte = ((unsigned int)bitpos) & 7;
    const unsigned long int chunk = (*(unsigned int*)bytes) >> bit_offset_in_first_byte;
	*p1 = (int)chunk & 127;
	*s1 = (int)((chunk >> 7) & 15) ^ 0x0A;
	*p2 = (int)(chunk >> 11) & 127;
	*s2 = (int)((chunk >> 18) & 15) ^ 0x0A;
}

void EERDecode_8bit(unsigned char b0, unsigned char b1, unsigned char b2,
    int *p1, int *s1, int *p2, int *s2)
{
    *p1 = b0;
    *s1 = (int)(b1 & 0x0F) ^ 0x0A;
    *p2 = (int)(b1 >> 4) | (b2 << 4);
    *s2 = (int)(b2 >> 4) ^ 0x0A;
}

void EERdecodePos4K(int p, int *x, int *y)
{
	*x = (p & 4095) + 1;
	*y = (p >> 12)  + 1;
}

void EERdecodePos8K(int p, int s, int *x, int *y)
{
    *x = (((p & 4095) << 1) | ((s & 2) >> 1)) + 1;
	*y = (((p >> 12)  << 1) | ((s & 8) >> 3)) + 1;
}

void EERdecodePos16K(int p, int s, int *x, int *y)
{
	*x = (((p & 4095) << 2) |  (s & 3)) + 1;
    *y = (((p >> 12)  << 2) | ((s & 12) >> 2)) + 1;
}

void TIFFPrintInfo(TIFF *my_tiff)
{
    TIFFPrintDirectory(my_tiff,stdout,TIFFPRINT_NONE);
}

int TIFFGetWidth(TIFF *my_tiff)
{
	int success, width_cast;
#ifdef TIFF430
    uint32_t width;
#else
	uint32 width;
#endif
	success = TIFFGetField(my_tiff,TIFFTAG_IMAGEWIDTH,&width);
	if ( success == 1 ) {
		width_cast = (int) width;
	} else {
		width_cast = 0;
	}
	return width_cast;
}

int TIFFGetLength(TIFF *my_tiff)
{
	int success, length_cast;
#ifdef TIFF430
    uint32_t length;
#else
    uint32 length;
#endif

	success = TIFFGetField(my_tiff,TIFFTAG_IMAGELENGTH,&length);
	if ( success == 1 ) {
		length_cast = (int) length;
	} else {
		length_cast = 0;
	}
	return length_cast;
}

int TIFFGetSampleFormat(TIFF *my_tiff)
{
    int success;
#ifdef TIFF430
    uint16_t format;
#else
    uint16 format;
#endif
    success = TIFFGetField(my_tiff,TIFFTAG_SAMPLEFORMAT,&format);
    return (int) format;
}

int TIFFGetCompression(TIFF *my_tiff)
{
    int success, compression_cast;
#ifdef TIFF430
    uint16_t compression;
#else
    uint16 compression;
#endif

    success = TIFFGetField(my_tiff,TIFFTAG_COMPRESSION,&compression);
	if ( success == 1 ) {
		compression_cast = (int)compression;
	} else {
		compression_cast = 0;
	}
    return compression_cast;
}

int TIFFNumDirectories(TIFF *my_tiff)
{
    int success, n_cast;
#ifdef TIFF430
    uint16_t n;
#else
    uint16 n;
#endif

    n = TIFFNumberOfDirectories(my_tiff);
	if ( n >= 1 ) {
		n_cast = (int)n;
	} else {
		n_cast = 0;
	}
    return n_cast;
}

int TIFFGetSamplesPerPixel(TIFF *my_tiff)
{
    int success;
#ifdef TIFF430
    uint16_t samples_per_pixel;
#else
    uint16 samples_per_pixel;
#endif
    success = TIFFGetField(my_tiff,TIFFTAG_SAMPLESPERPIXEL,&samples_per_pixel);
    return (int) samples_per_pixel;
}

int TIFFGetBitsPerSample(TIFF *my_tiff)
{
    int success;
#ifdef TIFF430
    uint16_t bits_per_sample;
#else
    uint16 bits_per_sample;
#endif
    success = TIFFGetField(my_tiff,TIFFTAG_BITSPERSAMPLE,&bits_per_sample);
    return (int) bits_per_sample;
}

int TIFFGetRowsPerStrip(TIFF *my_tiff)
{
    int success;
#ifdef TIFF430
    uint32_t rows_per_strip;
#else
    uint32 rows_per_strip;
#endif

    success = TIFFGetField(my_tiff,TIFFTAG_ROWSPERSTRIP,&rows_per_strip);
    return (int) rows_per_strip;
}

int TIFFGetMinVal(TIFF *my_tiff)
{
    int success;
#ifdef TIFF430
    uint32_t minval;
#else
    uint32 minval;
#endif
    success = TIFFGetField(my_tiff,TIFFTAG_MINSAMPLEVALUE,&minval);
    return (int) minval;
}

int TIFFGetMaxVal(TIFF *my_tiff)
{
    int success;
#ifdef TIFF430
    uint32_t maxval;
#else
    uint32 maxval;
#endif
    success = TIFFGetField(my_tiff,TIFFTAG_MINSAMPLEVALUE,&maxval);
    return (int) maxval;
}

tdata_t TIFFAllocateStripBuffer(TIFF *my_tiff)
{
    return _TIFFmalloc(TIFFStripSize(my_tiff));
}

tdata_t TIFFmalloc(int size)
{
    return _TIFFmalloc(size);
}

int TIFFfree(tdata_t buf)
{
    _TIFFfree(buf);
    return 0;
}

void TIFFMuteWarnings()
{
       warn = TIFFSetWarningHandler(0);
}

void TIFFUnMuteWarnings()
{
       TIFFSetWarningHandler(warn);
}
