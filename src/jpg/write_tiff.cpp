
/*
 * The TIFF writer from  http://paulbourke.net/dataformats/tiff/
 */

#include <sys/types.h>
#include <stdio.h>
#include <string.h>

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
//#include "write_tiff.h"

//Very simple TIFF writer based on
//  http://paulbourke.net/dataformats/tiff/

void WriteHexString(FILE *fptr,char *s)
{
   unsigned int i,c;
   char hex[3];
   for (i=0;i<strlen(s);i+=2) {
      hex[0] = s[i];
      hex[1] = s[i+1];
      hex[2] = '\0';
      sscanf(hex,"%X",&c);
      putc(c,fptr);
   }
} //end WriteHexString()

#define I_AM_LITTLE (((union { unsigned x; unsigned char c; }){1}).c)

bool isLittleEndian()
{
    short int number = 0x1;
    char *numPtr = (char*)&number;
    return (numPtr[0] == 1);
} //end isLittleEndian()

int normalize16 (uint16_t * img, int px){
//in general, TIFF readers tend to ignore the Maximum sample value tag
// they expect 16-bit images to be in the range 0..65535
// be warned this could screw up DICOM scale factors....
int mx = img[0];
for (int i = 0; i < px; i++)
	if (img[i] > mx) mx = img[i];
if (mx == 0) return mx;
double scale = 65535/mx;
for (int i = 0; i < px; i++)
	img[i] = scale * img[i];
return mx;
} //end normalize16()

unsigned char * swapImgBytes (unsigned char * img, int imgBytes) {
	unsigned char sw;
	int hi = 1;
	while (hi < imgBytes) {
		sw = img[hi];
		img[hi] = img[hi-1];
		img[hi-1] = sw;
		hi = hi + 2;
	}
	return img;
} //end swapImgBytes()

unsigned char * deplanar(unsigned char* img, int nx, int ny, int bits, int frames){
    //rrrr...rgggg...gbbbb...b to rgbrgbrgb...
    if ((nx < 1) || (ny < 1) || (bits != 8) || (frames != 3)) return img;
    int sliceBytes8 = nx*ny;
    unsigned char *slice24 = (unsigned char *)malloc( sliceBytes8*3 );
    int sliceOffsetR = 0;
    int sliceOffsetG = sliceBytes8;
    int sliceOffsetB = sliceBytes8+sliceBytes8;
	memcpy(&slice24[0], &img[0], sliceBytes8*3 );
	int i = 0;
	int j = 0;
	for (int rgb = 0; rgb < sliceBytes8; rgb++) {
		img[i] =slice24[sliceOffsetR+j];
		i++;
		img[i] =slice24[sliceOffsetG+j];
		i++;
		img[i] =slice24[sliceOffsetB+j];
		i++;
		j++;
	}
	free(slice24);
    return img;
} // nii_flipImgY()

int write_tiff_img (const char *fname, unsigned char *img, int nx, int ny, int bits, int frames, int isPlanar) {
	//if ((bits < 1) || (bits > 16)) return 1;
	int imgBytes = nx*ny*frames;
	if (imgBytes < 1) return 1;


	if (bits > 8) normalize16((uint16_t *) img, imgBytes);
	if (bits > 8) imgBytes = 2 * imgBytes;
	if (isLittleEndian() && (bits > 8)) img = swapImgBytes(img, imgBytes);
	if (isPlanar) img = deplanar(img, nx, ny, bits, frames);
	FILE* fptr = fopen(fname, "wb");
	if (!fptr) return 1;
	/* Write the header */
	WriteHexString(fptr,(char *)"4d4d002a");    /* Big endian & TIFF identifier */
	int offset = imgBytes + 8;
	putc((offset & 0xff000000) / 16777216,fptr);
	putc((offset & 0x00ff0000) / 65536,fptr);
	putc((offset & 0x0000ff00) / 256,fptr);
	putc((offset & 0x000000ff),fptr);
	/* Write the binary data */
	fwrite(img, 1, imgBytes, fptr);

	/* Write the footer */
	WriteHexString(fptr,(char *)"000e");  /* The number of directory entries (14) */
	/* Width tag, short int */
	WriteHexString(fptr,(char *)"0100000300000001");
	fputc((nx & 0xff00) / 256,fptr);    /* Image width */
	fputc((nx & 0x00ff),fptr);
	WriteHexString(fptr,(char *)"0000");
	/* Height tag, short int */
	WriteHexString(fptr,(char *)"0101000300000001");
	fputc((ny & 0xff00) / 256,fptr);    /* Image height */
	fputc((ny & 0x00ff),fptr);
	WriteHexString(fptr,(char *)"0000");

	/* Bits per sample tag, short int */
	if (frames == 3) {
		WriteHexString(fptr,(char *)"0102000300000003");
		offset = imgBytes + 182;
		putc((offset & 0xff000000) / 16777216,fptr);
		putc((offset & 0x00ff0000) / 65536,fptr);
		putc((offset & 0x0000ff00) / 256,fptr);
		putc((offset & 0x000000ff),fptr);
	} else {
		if (bits > 8)
			WriteHexString(fptr,(char *)"010200030000000100100000");
		else
			WriteHexString(fptr,(char *)"010200030000000100080000");
	}

	/* Compression flag, short int */
	WriteHexString(fptr,(char *)"010300030000000100010000");

	/* Photometric interpolation tag, short int */
	if (frames == 3)
		WriteHexString(fptr,(char *)"010600030000000100020000");//2 = RGB
	else //http://www.awaresystems.be/imaging/tiff/tifftags/photometricinterpretation.html
		WriteHexString(fptr,(char *)"010600030000000100010000");//1 = BlackIsZero
	/* Strip offset tag, long int */
	WriteHexString(fptr,(char *)"011100040000000100000008");

	/* Orientation flag, short int */
	WriteHexString(fptr,(char *)"011200030000000100010000");

	/* Sample per pixel tag, short int */
	if (frames == 3)
	WriteHexString(fptr,(char *)"011500030000000100030000");
	else
	WriteHexString(fptr,(char *)"011500030000000100010000");
	/* Rows per strip tag, short int */
	WriteHexString(fptr,(char *)"0116000300000001");
	fputc((ny & 0xff00) / 256,fptr);
	fputc((ny & 0x00ff),fptr);
	WriteHexString(fptr,(char *)"0000");

	/* Strip byte count flag, long int */
	WriteHexString(fptr,(char *)"0117000400000001");
	offset = imgBytes;
	putc((offset & 0xff000000) / 16777216,fptr);
	putc((offset & 0x00ff0000) / 65536,fptr);
	putc((offset & 0x0000ff00) / 256,fptr);
	putc((offset & 0x000000ff),fptr);

	/* Minimum sample value flag, short int */
	if (frames == 3)  {
		WriteHexString(fptr,(char *)"0118000300000003");
		offset = imgBytes + 188;
		putc((offset & 0xff000000) / 16777216,fptr);
		putc((offset & 0x00ff0000) / 65536,fptr);
		putc((offset & 0x0000ff00) / 256,fptr);
		putc((offset & 0x000000ff),fptr);
	} else
		WriteHexString(fptr,(char *)"011800030000000100000000");

	/* Maximum sample value tag, short int */
	if (frames == 3)  {
		WriteHexString(fptr,(char *)"0119000300000003");
		offset = imgBytes + 194;
		putc((offset & 0xff000000) / 16777216,fptr);
		putc((offset & 0x00ff0000) / 65536,fptr);
		putc((offset & 0x0000ff00) / 256,fptr);
		putc((offset & 0x000000ff),fptr);
	} else {
		if (bits > 8) {
			WriteHexString(fptr,(char *)"0119000300000001FFFF0000");
		}
		else
			WriteHexString(fptr,(char *)"011900030000000100FF0000");
	}

	/* Planar configuration tag, short int */
	WriteHexString(fptr,(char *)"011c00030000000100010000");

	/* Sample format tag, short int */
	if (frames == 3)  {
		WriteHexString(fptr,(char *)"0153000300000003");
		offset = imgBytes + 200; //200
		putc((offset & 0xff000000) / 16777216,fptr);
		putc((offset & 0x00ff0000) / 65536,fptr);
		putc((offset & 0x0000ff00) / 256,fptr);
		putc((offset & 0x000000ff),fptr);
	} else
		WriteHexString(fptr,(char *)"015300030000000100010000");

	/* End of the directory entry */
	WriteHexString(fptr,(char *)"00000000");
	if ( frames == 3  ) {
	   /* Bits for each colour channel */
	   WriteHexString(fptr,(char *)"000800080008");
	   /* Minimum value for each component */
	   WriteHexString(fptr,(char *)"000000000000");
	   /* Maximum value per channel */
	   WriteHexString(fptr,(char *)"00ff00ff00ff");
	   /* Samples per pixel for each channel */
	   WriteHexString(fptr,(char *)"000100010001");
	}
 return 0;
} //end write_tiff_img()
