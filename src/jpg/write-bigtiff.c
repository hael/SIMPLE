
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

#include <tiffio.h>

#define BIGTIFF_WRITE bigtiff_write_

int
bigtiff_write (char * fname, int*nx, int*ny, int*nz, float *mnval, float *mxval, void*imgptr)
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
