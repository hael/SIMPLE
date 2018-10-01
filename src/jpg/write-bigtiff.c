
/* compile with
 * 	gcc -g -Wall write-bigtiff.c `pkg-config libtiff-4 --cflags --libs` -lm
 */
#include <sys/types.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <string.h>
#include <math.h>
#include <tiff.h>
#include <tiffio.h>

#define BIGTIFF_WRITE bigtiff_write_

int
bigtiff_write (int*nx, int*ny, int*nz, float *mnval, float *mxval, void*imgptr, char * fname ,int* charStringLen , size_t ivf_CharStringLen)
{
	int width = *nx;
	int height = *ny;
  int slices = *nz;
	const double pixels_per_cm = 100.0;
  float range[2] = { *mnval,  *mxval};
	uint16_t* data; //unsigned short *data;
  float* img= imgptr;
	data = malloc((size_t) width * height * 2);

	int cx = width / 2;
	int cy = height / 2;
	int x, y;

	for (y = 0; y < height; y++)
		for (x = 0; x < width; x++) {
      float d =  img[y * width + x];
      d = (d - range[0])/(range[1]-range[0]);
			data[y * width + x] = (int) d * 65534;
		}

	TIFF *tif;

  tif = TIFFOpen(fname, "w8");

  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH     , width        );
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH    , height       );
  TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE  , 16                    );
  TIFFSetField(tif, TIFFTAG_COMPRESSION    , COMPRESSION_LZW       );
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC    , PHOTOMETRIC_MINISBLACK);
  TIFFSetField(tif, TIFFTAG_FILLORDER      , FILLORDER_MSB2LSB     );
  TIFFSetField(tif, TIFFTAG_ORIENTATION    , ORIENTATION_TOPLEFT   );
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1 );
  TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP , 8   );
  TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT , RESUNIT_CENTIMETER    );
  TIFFSetField(tif, TIFFTAG_XRESOLUTION    , pixels_per_cm);
  TIFFSetField(tif, TIFFTAG_YRESOLUTION    , pixels_per_cm);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG   , PLANARCONFIG_CONTIG   );

  tdata_t       l_buffer;
  unsigned long l_buffer_size = width*2; // LINE buffer for 16bit

  l_buffer = _TIFFmalloc(l_buffer_size);

  for (int row = 0; row < height; row++){
    memcpy(l_buffer, &data[row*width], l_buffer_size);

    int ret=TIFFWriteScanline(tif, l_buffer, row, 0);
    if (ret==-1){
      TIFFClose(tif);
      free(data);
      return -1;
    }
  }

  _TIFFfree(l_buffer);
  free(data);
  TIFFClose(tif);

	return 0;
}

static void
swapBytesInScanline(void *buf, uint32 width, TIFFDataType dtype)
{
	switch (dtype) {
		case TIFF_SHORT:
		case TIFF_SSHORT:
			TIFFSwabArrayOfShort((uint16*)buf,
                                             (unsigned long)width);
			break;
		case TIFF_LONG:
		case TIFF_SLONG:
			TIFFSwabArrayOfLong((uint32*)buf,
                                            (unsigned long)width);
			break;
		/* case TIFF_FLOAT: */	/* FIXME */
		case TIFF_DOUBLE:
			TIFFSwabArrayOfDouble((double*)buf,
                                              (unsigned long)width);
			break;
		default:
			break;
	}
}

// Multi page tiff writer based on
// https://raymondlo84.blogspot.com/2015/09/how-to-write-multipage-tiff-file.html
int tiffwriter(char*fname, int *dim1, int *dim2, int *dim3, float * imgptr)
{
    TIFF *out = TIFFOpen(fname,"w") ;
    uint32 imagelength = *dim1;
    uint32 imagewidth = *dim2;

    if (out)
    {
        int NPAGES = *dim3;
        int page;
        for (page = 0; page < NPAGES; page++){
            uint8 * buf;
            uint32 row, col, n, pixel;
            uint16 config, nsamples = 3;
            config = PLANARCONFIG_CONTIG ;

            TIFFSetField(out, TIFFTAG_IMAGELENGTH, imagelength);
            TIFFSetField(out, TIFFTAG_IMAGEWIDTH, imagewidth);
            TIFFSetField(out, TIFFTAG_PLANARCONFIG, config);
            TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, nsamples);
            TIFFSetField(out, TIFFTAG_COMPRESSION, COMPRESSION_LZW) ;
            TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8) ;
            TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out, imagewidth*nsamples));

            /* We are writing single page of the multipage file */
            TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
            TIFFSetField(out, TIFFTAG_PAGENUMBER, page, NPAGES);

            printf("writing %d x %d, nsamples %d", imagewidth, imagelength, nsamples);

            buf = _TIFFmalloc(imagewidth*nsamples);

            for (row = 0; row < imagelength; row++){

                for(col=0; col < imagewidth; col++){
                  pixel = (uint32) imgptr[page*imagelength*imagewidth + row*imagewidth+ col];
                  for(n = 0 ; n < nsamples ; ++n)
                    {

                      buf[col*nsamples+n] = (uint8) ( ( pixel >> (24 - n*8) ) & 255);
                    }
                }
                if (TIFFWriteScanline(out, buf, row, 0) != 1 )
                {
                    printf("Unable to write a row\n") ;
                    break ;
                }
            }
            TIFFWriteDirectory(out);
            _TIFFfree(buf);
        }
        TIFFClose(out);
    }
}

/* static	uint16 compression = (uint16) -1; */
/* static	int jpegcolormode = JPEGCOLORMODE_RGB; */
/* static	int quality = 75;		/\* JPEG quality *\/ */
/* static	uint16 predictor = 0; */
/* int */
/* raw2tiff(float*in,  int* width_, int*length_, int*slices_, */
/*          char* outfilename, */
/*          int* charStringLen, */
/*          size_t ivf_CharStringLen) */
/* { */
/* 	uint32 slices, linebytes, bufsize,	width = 0, length = 0; */
/* 	uint32	nbands = 1;		    /\* number of bands in input image*\/ */
/* 	_TIFF_off_t hdr_size = 0;	    /\* size of the header to skip *\/ */
/* 	TIFFDataType dtype = TIFF_FLOAT; */
/* 	int16	depth = 3;		    /\* bytes per pixel in input image *\/ */
/* 	int	swab = 0;		    /\* byte swapping flag *\/ */
/* 	InterleavingType interleaving = 0;  /\* interleaving type flag *\/ */
/* 	uint32  rowsperstrip = (uint32) -1; */
/* 	uint16	photometric = PHOTOMETRIC_MINISBLACK; */
/* 	uint16	config = PLANARCONFIG_CONTIG; */
/* 	uint16	fillorder = FILLORDER_LSB2MSB; */
/* 	int	fd; */
/* 	char	*outfilename = NULL; */
/* 	TIFF	*out; */

/* 	uint32 row, col, band; */
/* 	int	c; */
/* 	unsigned char *buf = NULL, *buf1 = NULL; */
/*   width=* width_; */
/*   length=*length_; */
/*   slices=*slices_; */
/* 	out = TIFFOpen(outfilename, "w"); */
/* 	if (out == NULL) { */
/* 		fprintf(stderr, "raw2tiff: %s: Cannot open file for output.\n", outfilename); */
/* 		return (-1); */
/* 	} */

/* 	TIFFSetField(out, TIFFTAG_IMAGEWIDTH, width); */
/* 	TIFFSetField(out, TIFFTAG_IMAGELENGTH, length); */
/* 	TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT); */
/* 	TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, nbands); */
/* 	TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, depth * 8); */
/* 	TIFFSetField(out, TIFFTAG_FILLORDER, fillorder); */
/* 	TIFFSetField(out, TIFFTAG_PLANARCONFIG, config); */
/* 	TIFFSetField(out, TIFFTAG_PHOTOMETRIC, photometric); */
/* 	switch (dtype) { */
/* 	case TIFF_BYTE: */
/* 	case TIFF_SHORT: */
/* 	case TIFF_LONG: */
/* 		TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT); */
/* 		break; */
/* 	case TIFF_SBYTE: */
/* 	case TIFF_SSHORT: */
/* 	case TIFF_SLONG: */
/* 		TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_INT); */
/* 		break; */
/* 	case TIFF_FLOAT: */
/* 	case TIFF_DOUBLE: */
/* 		TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP); */
/* 		break; */
/* 	default: */
/* 		TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_VOID); */
/* 		break; */
/* 	} */
/* 	if (compression == (uint16) -1) */
/* 		compression = COMPRESSION_PACKBITS; */
/* 	TIFFSetField(out, TIFFTAG_COMPRESSION, compression); */
/* 	switch (compression) { */
/* 	case COMPRESSION_JPEG: */
/* 		if (photometric == PHOTOMETRIC_RGB */
/* 		    && jpegcolormode == JPEGCOLORMODE_RGB) */
/* 			photometric = PHOTOMETRIC_YCBCR; */
/* 		TIFFSetField(out, TIFFTAG_JPEGQUALITY, quality); */
/* 		TIFFSetField(out, TIFFTAG_JPEGCOLORMODE, jpegcolormode); */
/* 		break; */
/* 	case COMPRESSION_LZW: */
/* 	case COMPRESSION_DEFLATE: */
/* 		if (predictor != 0) */
/* 			TIFFSetField(out, TIFFTAG_PREDICTOR, predictor); */
/* 		break; */
/* 	} */
/* 	switch(interleaving) { */
/* 	case BAND:				/\* band interleaved data *\/ */
/* 		linebytes = width * depth; */
/* 		buf = (unsigned char *)_TIFFmalloc(linebytes); */
/* 		break; */
/* 	case PIXEL:				/\* pixel interleaved data *\/ */
/* 	default: */
/* 		linebytes = width * nbands * depth; */
/* 		break; */
/* 	} */
/* 	bufsize = width * nbands * depth; */
/* 	buf1 = (unsigned char *)_TIFFmalloc(bufsize); */

/* 	rowsperstrip = TIFFDefaultStripSize(out, rowsperstrip); */
/* 	if (rowsperstrip > length) { */
/* 		rowsperstrip = length; */
/* 	} */
/* 	TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, rowsperstrip ); */

/* 	_TIFF_lseek_f(fd, hdr_size, SEEK_SET);		/\* Skip the file header *\/ */
/* 	for (row = 0; row < length; row++) { */
/* 		switch(interleaving) { */
/* 		case BAND:			/\* band interleaved data *\/ */
/* 			for (band = 0; band < nbands; band++) { */
/* 				if (_TIFF_lseek_f(fd, */
/*                                           hdr_size + (length*band+row)*linebytes, */
/*                                           SEEK_SET) == (_TIFF_off_t)-1) { */
/*                                         fprintf(stderr, */
/*                                                 "%s: %s: scanline %lu: seek error.\n", */
/*                                                 argv[0], argv[optind], */
/*                                                 (unsigned long) row); */
/*                                         break; */
/*                                 } */
/* 				if (read(fd, buf, linebytes) < 0) { */
/* 					fprintf(stderr, */
/*                                                 "%s: %s: scanline %lu: Read error.\n", */
/*                                                 argv[0], argv[optind], */
/*                                                 (unsigned long) row); */
/*                                         break; */
/* 				} */
/* 				if (swab)	/\* Swap bytes if needed *\/ */
/* 					swapBytesInScanline(buf, width, dtype); */
/* 				for (col = 0; col < width; col++) */
/* 					memcpy(buf1 + (col*nbands+band)*depth, */
/* 					       buf + col * depth, depth); */
/* 			} */
/* 			break; */
/* 		case PIXEL:			/\* pixel interleaved data *\/ */
/* 		default: */
/* 			if (read(fd, buf1, bufsize) < 0) { */
/* 				fprintf(stderr, */
/* 					"%s: %s: scanline %lu: Read error.\n", */
/* 					argv[0], argv[optind], */
/* 					(unsigned long) row); */
/* 				break; */
/* 			} */
/* 			if (swab)		/\* Swap bytes if needed *\/ */
/* 				swapBytesInScanline(buf1, width, dtype); */
/* 			break; */
/* 		} */

/* 		if (TIFFWriteScanline(out, buf1, row, 0) < 0) { */
/* 			fprintf(stderr,	"%s: %s: scanline %lu: Write error.\n", */
/* 				argv[0], outfilename, (unsigned long) row); */
/* 			break; */
/* 		} */
/* 	} */
/* 	if (buf) */
/* 		_TIFFfree(buf); */
/* 	if (buf1) */
/* 		_TIFFfree(buf1); */
/* 	TIFFClose(out); */
/* 	return (0); */
/* } */
