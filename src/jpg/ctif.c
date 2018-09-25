/*
 * Created by Michael Eager, 2018
 * Elmlund Lab, Monash University
 */

#include <sys/types.h>
#include <stdio.h>
#include "tiffio.h"


void Tifflib_readRGBA( float*realimg, int*width, int*height,
                       char * fname,
                       int* charStringLen,
                       size_t ivf_CharStringLen)
{
  TIFF* tif = TIFFOpen(fname, "r");
  if (tif) {
    TIFFRGBAImage img;
    char emsg[1024];

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
      TIFFError('Tifflib_read ', emsg);
    }
    TIFFClose(tif);
  }
  TIFFClose(tif);
}

void Tifflib_readScanline(float*realimg, int*width, int*height,
                          char * fname,
                          int* charStringLen,
                          size_t ivf_CharStringLen)
{
  TIFF* tif = TIFFOpen(fname, "r");
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

void Tifflib_readTiledImage(char * fname, float*realimg, int*width, int*height, int*xtiles, int*ytiles)
{
  TIFF* tif = TIFFOpen(fname, "r");
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


void Tifflib_writeRGBA(char * fname, float*img, int*width, int*height)
{
  int sampleperpixel = 4;    // or 3 if there is no alpha channel
  int w=*width,h=*height;
  unsigned char *image;

  image = malloc(w*h*sampleperpixel*sizeof(char) );

  for (int y = 0; y < h; y += 1){
    for (int x = 0; x < w; x += 1){
      unsigned int  pixel = img[y*w+x];
      image[y*w + x + 0] =  (0xff & pixel);
      image[y*w + x + 1] =  (0xff & (pixel >> 8));
      image[y*w + x + 2] =  (0xff & (pixel >> 16));
      image[y*w + x + 3] =  (0xff & (pixel >> 24));
    }}
  TIFF* out = TIFFOpen(fname, "w");
  TIFFSetField (out, TIFFTAG_IMAGEWIDTH, width);  // set the width of the image
  TIFFSetField(out, TIFFTAG_IMAGELENGTH, height);    // set the height of the image
  TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);   // set number of channels per pixel
  TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8);    // set the size of the channels
  TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
  //   Some other essential fields to set that you do not have to understand for now.
  TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

  tsize_t linebytes = sampleperpixel * w;     // length in memory of one row of pixel in the image.
  unsigned char *buf = NULL;        // buffer used to store the row of pixel information for writing to file
  //    Allocating memory to store the pixels of current row
  if (TIFFScanlineSize(out) != linebytes)
    buf =(unsigned char *)_TIFFmalloc(linebytes);
  else
    buf = (unsigned char *)_TIFFmalloc(TIFFScanlineSize(out));

  // We set the strip size of the file to be size of one row of pixels
  TIFFSetField(out, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out, w*sampleperpixel));

  //Now writing image to the file one strip at a time
  for (uint32 row = 0; row < h; row++)
  {
    memcpy(buf, &image[(h-row-1)*linebytes], linebytes);    // check the index here
    if (TIFFWriteScanline(out, buf, row, 0) < 0)
      break;
  }

  TIFFClose(out);
}
