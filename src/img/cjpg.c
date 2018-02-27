
/* // gcc -std=c11 -ljpeg ... */

/* typedef unsigned char Pixel; */

/* typedef struct { */
/*  unsigned int width, height; */
/*  Pixel *data; */
/* } RawImage; */
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "jpeglib.h"


struct jpeg_compress_struct  encInfo;
struct jpeg_decompress_struct decInfo;
struct jpeg_error_mgr errMgr;
J_COLOR_SPACE colorSpace = JCS_RGB;
int numComponents = 3;
int quality = 80;
struct stat file_info;
unsigned long jpg_size;
unsigned char *jpg_buffer;


void JpegEnc_init();
void JpegEnc_destroy();
void JpegDec_init();
void JpegDec_destroy();
/**
 * Sets the color space of the images to be worked with.
 *
 * @param colorSpace     the color space. Can be: JCS_GRAYSCALE, JCS_RGB, JCS_YCbCr, or JCS_CMYK.
 * @param numComponents  the number of components in the color space specified
 */
//    void setColorSpace(J_COLOR_SPACE colorSpace, int numComponents);

/**
 * Gets the current color space.
 *
 * @return  the color space identificator
 */
//    J_COLOR_SPACE getColorSpace() const;

/**
 * Gets the number of components in the current color space.
 *
 * @return  the number of color space components
 */
//  int getNumComponents() const;

/**
 * Sets a compression quality.
 *
 * @param quality  the quality in percents [0, 100]
 */
// void setQuality(int quality);

/**
 * Gets the current compression quality.
 *
 * @return  the quality value
 */
// int getQuality() const;

/**
 * Encodes a raw image to a JPEG file on disk using file descriptor.
 *
 * @param img     the pointer to buffer with the raw image
 * @param width   the width of the image
 * @param height  the height of the image
 * @param file    the file descriptor
 * @return        true if success, false if failed
 */
bool JpegEnc_encode1(uint8_t* img, int width, int height, FILE* file);

/**
 * Encodes a raw image to a JPEG file on disk using path.
 *
 * @param img     the pointer to buffer with the raw image
 * @param width   the width of the image
 * @param height  the height of the image
 * @param path    the path to the file
 * @return        true if success, false if failed
 */
bool JpegEnc_encode(uint8_t* img, int width, int height, const char* path);

/*private:*/
bool JpegEnc_encode0(uint8_t* img, int width, int height);
//};

/**
 * JPEG decoder.
 */
/* class JpegDec { */
/* public: */
/*     JpegDec(); */
/*     virtual ~JpegDec(); */

/*     JpegDec(const JpegDec&) = delete; */
/*     void operator=(const JpegDec&) = delete; */

/**
 * Sets the color space of the images to be worked with.
 *
 * @param colorSpace  the color space. Can be: JCS_GRAYSCALE, JCS_RGB, JCS_YCbCr, or JCS_CMYK.
 */
void setColorSpace(J_COLOR_SPACE colorSpace);

/**
 * Gets the current color space.
 *
 * @return  the color space identifier
 */
// J_COLOR_SPACE getColorSpace() const;

/**
 * Decodes a JPEG image from a memory buffer.
 *
 * @param[in]  buffer  the pointer to the buffer with encoded image
 * @param[in]  len     the size of the buffer
 * @param[out] img     the decoded raw image. The memory is allocated inside the function.
 * @param[out] width   the width of the decoded image
 * @param[out] height  the height of the decoded image
 * @return             true if success, false if failed
 */
//    bool JpegDec_decode2(uint8_t* buffer, size_t len, uint8_t*& img, int& width, int& height);

/**
 * Reads and decodes a JPEG image through file descriptor.
 *
 * @param[in]  file    the file descriptor
 * @param[out] img     the decoded raw image. The memory is allocated inside the function.
 * @param[out] width   the width of the decoded image
 * @param[out] height  the height of the decoded image
 * @return             true if success, false if failed
 */
bool JpegDec_decode1(FILE* file, uint8_t** img, int* width, int* height);

/**
 * Reads and decodes a JPEG file on disk.
 *
 * @param[in]  path    the path to the file
 * @param[out] img     the decoded raw image. The memory is allocated inside the function.
 * @param[out] width   the width of the decoded image
 * @param[out] height  the height of the decoded image
 * @return             true if success, false if failed
 */
bool JpegDec_decode(const char* path, uint8_t** img, int* width, int* height);

//private:
bool JpegDec_decode0(uint8_t** img, int* width, int* height);


//    jpeg_error_mgr errMgr;

//    J_COLOR_SPACE colorSpace = JCS_RGB;
//};

/**
 * Gets a number of color components for the given color space.
 *
 * @param colorSpace  the color space
 * @return            the number of components for supported color space; 0 for unsupported
 */


int getNumComponents(J_COLOR_SPACE colorSpace)
{
    switch(colorSpace) {
    case JCS_GRAYSCALE:
        return 1;
    case JCS_RGB:
    case JCS_YCbCr:
        return 3;
    case JCS_CMYK:
        return 4;
    default:
        return 0;
    }
}

void
JpegEnc_init()
{
    encInfo.err = jpeg_std_error(&errMgr);
    jpeg_create_compress(&encInfo);
}
void
JpegEnc_destroy()
{
    jpeg_destroy_compress(&encInfo);
}

//void JpegEnc_setColorSpace(J_COLOR_SPACE colorSpace_, int numComponents_) {
//    colorSpace = colorSpace_;
//    numComponents = numComponents_;
//}

J_COLOR_SPACE JpegEnc_getColorSpace()
{
    return colorSpace;
}

int JpegEnc_getNumComponents()
{
    return numComponents;
}

void JpegEnc_setQuality(int quality_)
{
    quality = quality_;
}

int JpegEnc_getQuality()
{
    return quality;
}

bool JpegEnc_encode1(uint8_t* img, int width, int height, FILE* file)
{
    if(file == NULL) {
        return false;
    }

    jpeg_stdio_dest(&encInfo, file);
    return JpegEnc_encode0(img, width, height);
}

bool JpegEnc_encode(uint8_t* img, int width, int height, const char* path)
{
    FILE* file = fopen(path, "wb");
    bool rv = JpegEnc_encode1(img, width, height, file);
    fclose(file);
    return rv;
}

bool JpegEnc_encode0(uint8_t* img, int width, int height)
{
    encInfo.image_width = width;
    encInfo.image_height = height;
    encInfo.in_color_space = colorSpace;
    encInfo.input_components = numComponents;
    jpeg_set_defaults(&encInfo);
    jpeg_set_quality(&encInfo, quality, true);

    jpeg_start_compress(&encInfo, true);

    int rowSize = encInfo.image_width * encInfo.input_components;
    JSAMPROW rows[1];
    while(encInfo.next_scanline < encInfo.image_height) {
        rows[0] = img + encInfo.next_scanline * rowSize;
        jpeg_write_scanlines(&encInfo, rows, 1);
    }

    jpeg_finish_compress(&encInfo);
    return true;
}

void JpegDec_init()
{
    decInfo.err = jpeg_std_error(&errMgr);
    jpeg_create_decompress(&decInfo);
}

void JpegDec_destroy()
{
    jpeg_destroy_decompress(&decInfo);
}

void JpegDec_setColorSpace(J_COLOR_SPACE colorSpace_)
{
    colorSpace = colorSpace_;
}

J_COLOR_SPACE JpegDec_getColorSpace()
{
    return colorSpace;
}

/* bool JpegDec_decode2(uint8_t* buffer,size_t length,uint8_t*& img, int& width, int& height) { */
/*     if (buffer == nullptr) { */
/*         return false; */
/*     } */

/*     jpeg_mem_src(&decInfo, buffer, length); */
/*     return  JpegDec_decode0(img, width, height); */
/* } */

bool JpegDec_decode1(FILE* file, uint8_t** img, int* width, int* height)
{
    if(file == NULL) {
        return false;
    }

    jpeg_stdio_src(&decInfo, file);
    return  JpegDec_decode0(img, width, height);
}

bool JpegDec_decode(const char* path, uint8_t** img, int* width, int* height)
{
    FILE* file = fopen(path, "rb");
    bool rv = JpegDec_decode1(file, img, width, height);
    fclose(file);
    return rv;
}

bool JpegDec_decode0(uint8_t** img, int *width, int* height)
{
    jpeg_read_header(&decInfo, true);
    decInfo.out_color_space = colorSpace;
    decInfo.raw_data_out = false;

    jpeg_start_decompress(&decInfo);

    *width = decInfo.image_width;
    *height = decInfo.image_height;
    size_t rowSize = decInfo.image_width * decInfo.num_components;
    *img = (uint8_t*)malloc(decInfo.image_height * rowSize * sizeof(uint8_t));

    JSAMPROW rows[1];
    size_t offset = 0;
    while(decInfo.output_scanline < decInfo.image_height) {
        rows[0] = *img + offset;
        jpeg_read_scanlines(&decInfo, rows, 1);
        offset += rowSize;
    }

    jpeg_finish_decompress(&decInfo);
    return true;
}


#define CWRITE_JPEG cwrite_jpeg_
#define CREAD_JPEG cread_jpeg_

void cwrite_jpeg(uint8_t** img,  const char* file_name, int* width, int* height, int* quality_, int* colorspec, bool* status)
{
    JpegEnc_init();
    *status = false;
    if(file_name == NULL) return;
    //switch(colorspec){
    //  case JCS_GRAYSCALE: j.setColorSpace(JCS_GRAYSCALE,j.getNumComponents());break;
    //  default: // j.setColorSpace(JCS_RGB,j.getNumComponents());
    // || colorspec == JCS_RGB ){ //  , JCS_YCbCr|| JCS_CMYK
    //  j.setColorSpace(colorspec);
    // }
    quality = *quality_;
    *status = JpegEnc_encode(*img, *width, *height, file_name);
    JpegEnc_destroy();
}

void cread_jpeg(const char* file_name, uint8_t* img, int* width, int* height, int* colorspec,  bool* status)
{
    *status = false;
    if(file_name == NULL) return;
    JpegDec_init();
    *status = JpegDec_decode(file_name, &img, width, height);
    *colorspec = getNumComponents(colorSpace);
    JpegDec_destroy();
}

#define WRITE_JPEG write_jpeg_
#define READ_JPEG read_jpeg_


#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
void read_jpeg(const char* file_name, int** out_buffer, int*width, int*height, int* colorspec,  bool* status)
{
    struct stat file_info;
    int rc, i, pixel_size;
    unsigned long bmp_size;
    unsigned char * buffer;
    int row_stride;
    int fd;
    *status = false;
    rc = stat(file_name, &file_info);
    fprintf(stderr, "Read_jpeg stat input file %s", file_name);

    if(rc){
      fprintf(stderr, "stat failed ");
      exit(EXIT_FAILURE);
    }
    jpg_size = file_info.st_size;
    jpg_buffer = (unsigned char*) malloc(jpg_size + 100);
    if (jpg_buffer == NULL){
      fprintf(stderr, "Failed to allocate bytes %d", (int)(jpg_size+100));
      return;
    }
    fprintf(stderr, "Read_jpeg stat input file size %d", (int)jpg_size);

    fd = open(file_name, O_RDONLY);
    if(fd<0) {
      fprintf(stderr, "Failed to open input file %s", file_name);
      return ;
    }
    i = 0;
    while(i < jpg_size) {
        rc = read(fd, jpg_buffer + i, jpg_size - i);
        i += rc;
    }
    close(fd);
    fprintf(stderr, "Read_jpeg decompressing ");

    decInfo.err = jpeg_std_error(&errMgr);
    jpeg_create_decompress(&decInfo);
    jpeg_mem_src(&decInfo, jpg_buffer, jpg_size);
    rc = jpeg_read_header(&decInfo, TRUE);
    fprintf(stderr, "Read_jpeg reading header %d ", rc);

    if(rc != 1) exit(EXIT_FAILURE);
    fprintf(stderr, "Read_jpeg start decompress");

    jpeg_start_decompress(&decInfo);

    *width = decInfo.output_width;
    *height = decInfo.output_height;
    pixel_size = decInfo.output_components;
    bmp_size = (*width) * (*height) * pixel_size;
    buffer = (unsigned char*) malloc(bmp_size);
    fprintf(stderr, "Read_jpeg buffer allocated %d", bmp_size);
    if(buffer==NULL){
      fprintf(stderr, "Read_jpeg buffer allocate failed");
      exit(EXIT_FAILURE);
    }
    // The row_stride is the total number of bytes it takes to store an
    // entire scanline (row).
    row_stride = (*width) * pixel_size;
    fprintf(stderr, "Read_jpeg buffer row size %d\n", row_stride);

    while(decInfo.output_scanline < decInfo.output_height) {
        unsigned char *buffer_array[1];
        buffer_array[0] = buffer + \
                          (decInfo.output_scanline) * row_stride;

        jpeg_read_scanlines(&decInfo, buffer_array, 1);
        fprintf(stderr, " %d", decInfo.output_scanline);

    }
    fprintf(stderr, "Read_jpeg scanlines complete %d\n", decInfo.output_scanline);

    jpeg_finish_decompress(&decInfo);
    fprintf(stderr, "Read_jpeg finish\n");
    jpeg_destroy_decompress(&decInfo);
    fprintf(stderr, "Read_jpeg destroy\n");

    // And free the input buffer
    free(jpg_buffer);
    *colorspec = pixel_size;
    *status = true;
    out_buffer = buffer;
}


void write_jpeg(void** img, const char* filename, int* width, int* height, int* quality_, int* colorspec,  bool*status)
{
    FILE* outfile;
    unsigned long pixel;
    int x, y, depth = 24;
    char* buffer;
    JSAMPROW row_pointer;
    int w, h;
    int ** imgint = (int**) img;
    unsigned char* jbuf;
    int numerator, newsize, size;
    size = 256;
    for(numerator = 1; numerator <= 8; numerator++) {
        newsize = (*width) * numerator / 8;
        if(newsize > size)
            break;
    }

    outfile = fopen(filename, "wb");
    if(!outfile) {
        fprintf(stderr, "Failed to open output file %s", filename);
        return ;
    }
    printf( "File %s opened ", filename);
    w = *width; h = *height;
    printf("Width  %d   Height %d ", w, h);
    /* collect separate RGB values to a buffer */
    buffer = malloc(sizeof(char) * 3 * w * h);

    for(y = 0; y < *height; y++) {
        for(x = 0; x < *width; x++) {
            pixel = (unsigned long) imgint[x][y];
            if (pixel > 0xffffff) {
              pixel = 0xffffff;
            }else{
            buffer[y * (*width) * 3 + x * 3 + 0] = (char)(pixel >> 16);
            buffer[y * (*width) * 3 + x * 3 + 1] = (char)((pixel & 0x00ff00) >> 8);
            buffer[y * (*width) * 3 + x * 3 + 2] = (char)(pixel & 0x0000ff);
            }
        }
    }

    encInfo.err = jpeg_std_error(&errMgr);
    jpeg_create_compress(&encInfo);
    jpeg_stdio_dest(&encInfo, outfile);
    encInfo.image_width = *width;
    encInfo.image_height = *height;
    encInfo.input_components = 1;
    encInfo.in_color_space = JCS_GRAYSCALE;
    jpeg_set_defaults(&encInfo);
    encInfo.scale_num = numerator;
    encInfo.scale_denom = 8;
    jpeg_set_quality(&encInfo, *quality_, TRUE);
    jpeg_start_compress(&encInfo, TRUE);
    while(encInfo.next_scanline < encInfo.image_height) {
        row_pointer = (JSAMPROW) &buffer[encInfo.next_scanline * (depth >> 3) * (*width)];
        jpeg_write_scanlines(&encInfo, &row_pointer, 1);
    }
    jpeg_finish_compress(&encInfo);
    jpeg_destroy_compress(&encInfo);
    free(buffer);
    fclose(outfile);

    return ;
}
