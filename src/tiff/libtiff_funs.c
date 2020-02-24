// LibTIFF wrapper edited from Unblur
// Copyright 2014 Howard Hughes Medical Institute
// All rights reserved
// Use is subject to Janelia Farm Research Campus Software Copyright 1.1
// license terms ( http://license.janelia.org/license/jfrc_copyright_1_1.html )
#include <tiffio.h>

void TIFFPrintInfo(TIFF *my_tiff)
{
    TIFFPrintDirectory(my_tiff,stdout,TIFFPRINT_NONE);
}

int TIFFGetWidth(TIFF *my_tiff)
{
	int success;
	uint32 width;
	int width_cast;

	success = TIFFGetField(my_tiff,TIFFTAG_IMAGEWIDTH,&width);
	if ( success == 1 ) {
		width_cast = (int) width;
	}
	else {
		width_cast = 0;
	}
	return width_cast;
}

int TIFFGetLength(TIFF *my_tiff)
{
	int success;
	uint32 length;
	int length_cast;

	success = TIFFGetField(my_tiff,TIFFTAG_IMAGELENGTH,&length);
	if ( success == 1 ) {
		length_cast = (int) length;
	}
	else {
		length_cast = 0;
	}
	return length_cast;
}

int TIFFGetSampleFormat(TIFF *my_tiff)
{
    int success;
    uint16 format;
    success = TIFFGetField(my_tiff,TIFFTAG_SAMPLEFORMAT,&format);
    return (int) format;
}

int TIFFGetSamplesPerPixel(TIFF *my_tiff)
{
    int success;
    uint16 samples_per_pixel;
    success = TIFFGetField(my_tiff,TIFFTAG_SAMPLESPERPIXEL,&samples_per_pixel);
    return (int) samples_per_pixel;
}

int TIFFGetBitsPerSample(TIFF *my_tiff)
{
    int success;
    uint16 bits_per_sample;
    success = TIFFGetField(my_tiff,TIFFTAG_BITSPERSAMPLE,&bits_per_sample);
    return (int) bits_per_sample;
}

int TIFFGetRowsPerStrip(TIFF *my_tiff)
{
    int success;
    uint32 rows_per_strip;
    success = TIFFGetField(my_tiff,TIFFTAG_ROWSPERSTRIP,&rows_per_strip);
    return (int) rows_per_strip;
}

int TIFFGetMinVal(TIFF *my_tiff)
{
    int success;
    uint32 minval;
    success = TIFFGetField(my_tiff,TIFFTAG_MINSAMPLEVALUE,&minval);
    return (int) minval;
}

int TIFFGetMaxVal(TIFF *my_tiff)
{
    int success;
    uint32 maxval;
    success = TIFFGetField(my_tiff,TIFFTAG_MINSAMPLEVALUE,&maxval);
    return (int) maxval;
}

tdata_t TIFFAllocateStripBuffer(TIFF *my_tiff)
{
    return _TIFFmalloc(TIFFStripSize(my_tiff));
}

int TIFFDeallocateStripBuffer(tdata_t buf)
{
    _TIFFfree(buf);
    return 0;
}
