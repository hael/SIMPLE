/*
 * Created by Michael Eager, 2018
 * Elmlund Lab, Monash University
 */
#define _THREAD_SAFE
#define _GNU_SOURCE             /* See feature_test_macros(7) */
#ifdef __APPLE__
#define _DARWIN_C_SOURCE
#include <MacTypes.h>
#define _BSD_SOURCE
#define __DARWIN_C_SOURCE
#endif
#include <sys/types.h>
#include <sys/mman.h>
#include <errno.h>          /*     extern int errno;   */
#include <ctype.h>
#include <string.h> /* memset */
#include <unistd.h> /* close */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "tiffio.h"

#define TIFFLIB_READRGBA tifflib_readrgba_
#define TIFFLIB_WRITERGBA tifflib_writergba_
#define TIFFLIB_READSCANLINE tifflib_readscanline_
#define TIFFLIB_READTILEDIMAGE tifflib_readtiledimage_

char *FtoCstring(char *str, int len)
{
  char *res; // C arrays are from 0:len-1
  size_t n = (size_t) len;
  if(len > 1024 && len < 0) {
    fprintf(stderr, "Warning: simple_posix.c: possible strlen error in converting f90 string to c-string\n");
  }
  if((res = (char*)malloc(n + 1))) {
    strncpy(res, str, n);
    res[len] = '\0';
  }
  return res;
}

void tifflib_readrgba( float*realimg, int*width, int*height,
                       char * fname,
                       int* charStringLen,
                       size_t ivf_CharStringLen)
{
  char* cname = FtoCstring(fname, *charStringLen);
  TIFF* tif = TIFFOpen(cname, "r"); free(cname);
  if (tif) {
    TIFFRGBAImage img;
    char emsg[512];

    if (TIFFRGBAImageBegin(&img, tif, 0, emsg)) {
	    size_t npixels;
	    uint32* raster;
      *width=img.width; *height= img.height;
	    npixels = img.width * img.height;
	    raster = (uint32*) _TIFFmalloc(npixels * sizeof (uint32));
	    if (raster != NULL) {
        if (TIFFRGBAImageGet(&img, raster, img.width, img.height)) {
          for(register int i=0; i<npixels;i=i+1) realimg[i] = (float) raster[i];
        }
        _TIFFfree(raster);
	    }
	    TIFFRGBAImageEnd(&img);
    } else{
      TIFFError(fname, "ctif.c:Readrgba failed: %s", emsg);
    }
    TIFFClose(tif);
  }
  TIFFClose(tif);
}
/*

void tifflib_readscanline(float*realimg, int*width, int*height,
                          char * fname,
                          int* charStringLen,
                          size_t ivf_CharStringLen)
{ char* cname = FtoCstring(fname, *charStringLen);
  TIFF* tif = TIFFOpen(fname, "r");free(cname);
  if (tif) {
    uint32 imagelength, imagewidth;
    tdata_t buf;
    uint32 row;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imagewidth);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imagelength);
    *height = (int) imagelength;
    *width= (int) imagewidth;
    buf = _TIFFmalloc(TIFFScanlineSize(tif));
    for (row = 0; row < imagelength; row++){
      TIFFReadScanline(tif, buf, row,0);
      for(register int i=0; i<height;i=i+1) realimg[(row*imagewidth)+ i] = (float) buf[i];
    }
    _TIFFfree(buf);
    TIFFClose(tif);
  }

}

void tifflib_readtiledimage(char * fname, float*realimg, int*width, int*height, int*xtiles, int*ytiles,
                       int* charStringLen,
                       size_t ivf_CharStringLen)
{ char* cname = FtoCstring(fname, *charStringLen);
  TIFF* tif = TIFFOpen(fname, "r");free(cname);
  if (tif) {
    uint32 imageWidth, imageLength;
    uint32 tileWidth, tileLength;
    uint32 x, y;
    tdata_t buf;
    size_t npixels;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &imageWidth);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &imageLength);
    TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileWidth);
    TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileLength);
    buf = _TIFFmalloc(TIFFTileSize(tif));
    npixels = tileWidth*tileLength;
    for (y = 0; y < imageLength; y += tileLength){
      for (x = 0; x < imageWidth; x += tileWidth){
        TIFFReadTile(tif, buf, x, y, 0);
        for(register int i=0; i<npixels;i=i+1) realimg[(y*imageWidth)+x+i] = (float) buf[i];
      }
    }
    _TIFFfree(buf);
    TIFFClose(tif);

    *width=(int)tileWidth;
    *height=(int)tileLength;
    *xtiles=(int)imageWidth/tileWidth;
    *ytiles=(int)imageLength/tileLength;
  }
}
*/

void tifflib_writergba(unsigned char *image, int*width, int*height, int*crgb,char * fname,
                       int* charStringLen,
                       size_t ivf_CharStringLen)
{

  fprintf(stderr, " ctif.c: tifflib_writergba size %d, %d fname: %s fnamesz: %d\n",*width,*height,fname, * charStringLen);
  char* cname = FtoCstring(fname, (int)*charStringLen);
  fprintf(stderr, " ctif.c: tifflib_writergba Cstring :%s:\n", cname);
  int sampleperpixel = *crgb;    // or 3 if there is no alpha channel
  int w=*width,h=*height;
  //  float*img=img_ptr;
  //unsigned char *image = image_ptr;
  //fprintf(stderr, " ctif.c: tifflib_writergba %f scale %f    %10.8f \n",   img[0] ,*scale,   img[0]*(*scale));
  fprintf(stderr, " ctif.c: Allocating image \n");
  /* image = malloc(w*h*sampleperpixel*sizeof(char) ); */
  /* fprintf(stderr," ctif.c: Allocating image %zd", sizeof(image)); */
  /* for (int y = 0; y < h; y += 1){ */
  /*   for (int x = 0; x < w; x += 1){ */
  /*     unsigned int  pixel = (int) (*scale) * img[y*w+x]; */
  /*     image[y*w + x + 0] =  (0xff & pixel); */
  /*     image[y*w + x + 1] =  (0xff & (pixel >> 8)); */
  /*     image[y*w + x + 2] =  (0xff & (pixel >> 16)); */
  /*     if(sampleperpixel == 4) image[y*w + x + 3] =  (0xff & (pixel >> 24)); */
  /*   }} */
  TIFF* out = TIFFOpen(cname, "w"); free(cname);
  if(out) {
    printf(" ctif.c: TIFFOpen success \n" );
  }else{
    printf(" ctif.c: TIFFOpen failed to open file %s\n", fname );
    return ;
  }
  TIFFSetField(out, TIFFTAG_IMAGEWIDTH,  w);  // set the width of the image
  TIFFSetField(out, TIFFTAG_IMAGELENGTH, h);    // set the height of the image
  TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);   // set number of channels per pixel
  TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);    // set the size of the channels
  TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
  //   Some other essential fields to set that you do not have to understand for now.
  TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  fprintf(stderr," ctif.c: linebytes %d x %d\n", sampleperpixel , w);
  tsize_t linebytes = sampleperpixel * w;     // length in memory of one row of pixel in the image.
  unsigned char *buf = NULL;        // buffer used to store the row of pixel information for writing to file
  //    Allocating memory to store the pixels of current row
  if (TIFFScanlineSize(out) != linebytes)
    buf =(unsigned char *)_TIFFmalloc(linebytes);
  else
    buf = (unsigned char *)_TIFFmalloc(TIFFScanlineSize(out));
  fprintf(stderr," ctif.c: strip size %u buf:%ld tiffsclsz:%ld\n", TIFFDefaultStripSize(out, w*sampleperpixel), (long) sizeof(buf), (long) TIFFScanlineSize(out));
  // We set the strip size of the file to be size of one row of pixels
  TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out, w*sampleperpixel));
  printf("Address of 'image[(row)*linebytes]' %pn %pn  0 %12ld\n",&image[0], &image[h*linebytes-1], h*linebytes-1);
  //Now writing image to the file one strip at a time
  for (uint32 row = 0; row < h; row++)
  {
    fprintf(stderr," ctif.c: row %d ",row);
    printf("Address of 'image[(row)*linebytes]' 0:%pn %12ld:%pn  %12ld:%pn\n",&image[0],row*linebytes, &image[(row)*linebytes],h*linebytes-1, &image[h*linebytes-1]);
    _TIFFmemcpy(buf, &image[(row)*linebytes], linebytes);    // check the index here
    if (TIFFWriteScanline(out, buf, row, 0) < 0)
      break;
    fprintf(stderr," ctif.c: \n");
  }

  TIFFClose(out);
}
