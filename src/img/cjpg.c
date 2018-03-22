
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
int Jpeg_quality = 80;
int Jpeg_width = 80;
int Jpeg_height = 80;
int* Jpeg_img_buffer;
struct stat file_info;
unsigned long jpg_size;
unsigned char *jpg_buffer;



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



int JpegEnc_encode(uint8_t* img, int width, int height, const char* path)
{
    FILE* file = fopen(path, "wb");
    if(file == NULL) {
        return 1;
    }

    jpeg_stdio_dest(&encInfo, file);

    encInfo.image_width = width;
    encInfo.image_height = height;
    encInfo.in_color_space = colorSpace;
    encInfo.input_components = numComponents;
    jpeg_set_defaults(&encInfo);
    jpeg_set_quality(&encInfo, Jpeg_quality, true);

    jpeg_start_compress(&encInfo, true);

    int rowSize = encInfo.image_width * encInfo.input_components;
    JSAMPROW rows[1];
    while(encInfo.next_scanline < encInfo.image_height) {
        rows[0] = img + encInfo.next_scanline * rowSize;
        jpeg_write_scanlines(&encInfo, rows, 1);
    }

    jpeg_finish_compress(&encInfo);
    fclose(file);
    return 0;
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



int JpegDec_decode(const char * path, uint8_t** imgbuf, int *width, int* height)
{
    FILE* file = fopen(path, "rb");
    if(file == NULL) {
        return 1;
    }

    jpeg_stdio_src(&decInfo, file);
    jpeg_read_header(&decInfo, true);
    decInfo.out_color_space = colorSpace;
    decInfo.raw_data_out = false;
    jpeg_start_decompress(&decInfo);
    *width = decInfo.image_width;
    *height = decInfo.image_height;


    size_t rowSize = decInfo.image_width * decInfo.num_components;
    *imgbuf = malloc(decInfo.image_height * rowSize * sizeof(uint8_t));

    JSAMPROW rows[1];
    size_t offset = 0;
    while(decInfo.output_scanline < decInfo.image_height) {
        rows[0] = *imgbuf + offset;
        jpeg_read_scanlines(&decInfo, rows, 1);
        offset += rowSize;
    }

    jpeg_finish_decompress(&decInfo);
    fclose(file);
    return 0;
}


#define CWRITE_JPEG cwrite_jpeg_
#define CREAD_JPEG cread_jpeg_

void cwrite_jpeg(void** imgptr,  const char* file_name, int* width, int* height, int* quality, int* colorspec, int* status)
{
    unsigned int pixel;
    JpegEnc_init();
    *status = 1;
    if(file_name == NULL) return;
    //switch(colorspec){
    //  case JCS_GRAYSCALE: j.setColorSpace(JCS_GRAYSCALE,j.getNumComponents());break;
    //  default: // j.setColorSpace(JCS_RGB,j.getNumComponents());
    // || colorspec == JCS_RGB ){ //  , JCS_YCbCr|| JCS_CMYK
    //  j.setColorSpace(colorspec);
    // }
    if((J_COLOR_SPACE)*colorspec != JCS_GRAYSCALE) return;
    Jpeg_quality = *quality;
    uint8_t *img = malloc(sizeof(char) * (*height) * (*width) * 3);
    for(int y = 0; y < *height; y++) {
        for(int x = 0; x < *width; x++) {
            pixel = (unsigned int)((int**)imgptr)[x][y];
            if(pixel > 0xffffff) {
                pixel = 0xffffff;
            } else {
                img[y * (*width) * 3 + x * 3 + 0] = (uint8_t)(pixel >> 16);
                img[y * (*width) * 3 + x * 3 + 1] = (uint8_t)((pixel & 0x00ff00) >> 8);
                img[y * (*width) * 3 + x * 3 + 2] = (uint8_t)(pixel & 0x0000ff);
            }
        }
    }
    *status =  JpegEnc_encode(img, *width, *height, file_name);
    JpegEnc_destroy();
    *status = 0;
}

void cread_jpeg(const char* file_name, void** outptr, int* width, int* height, int* colorspec, int* status)
{
    uint8_t** img;
    size_t rows, cols;
    *status = 1;
    if(file_name == NULL) {
        printf("cread_jpeg: Filename is empty.\n");
        *status = 1;
        return;
    }
    JpegDec_init();
    *status =  JpegDec_decode(file_name, img, width, height);
    *colorspec = getNumComponents(colorSpace);
    JpegDec_destroy();


    cols = *height; rows = *width;
    int (*image)[cols] = malloc(sizeof(*image) * rows);
    // *(*image) = (int**)malloc(sizeof(int**)*(*height));
    if(! image) {
        printf("cread_jpeg: Image buffer is not allocated.\n");
        *status = 1;
        return;
    }
    for(int j = 0; j < *height; j++)
        for(int i = 0; i < *width; i++)
            image[i][j] = (int) img[i][j];
    free(img);
    *outptr = (void*) image;
    *status = 0;
}

#define WRITE_JPEG write_jpeg_
#define READ_JPEG2 read_jpeg_
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

void read_jpeg2(const char* file_name, int** out_buffer, int*width, int*height, int* colorspec,  int* status)
{
    struct stat file_info;
    int rc, i, j, pixel_size;
    unsigned long bmp_size, pixel;
    unsigned char * buffer;
    int row_stride;
    int fd;
    *status = 1;
    rc = stat(file_name, &file_info);
    fprintf(stderr, "Read_jpeg stat input file %s", file_name);

    if(rc) {
        fprintf(stderr, "stat failed ");
        exit(EXIT_FAILURE);
    }
    jpg_size = file_info.st_size;
    jpg_buffer =  malloc(jpg_size + 100);
    if(jpg_buffer == NULL) {
        fprintf(stderr, "Failed to allocate bytes %d", (int)(jpg_size + 100));
        return;
    }
    fprintf(stderr, "Read_jpeg stat input file size %d", (int)jpg_size);

    fd = open(file_name, O_RDONLY);
    if(fd < 0) {
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
    buffer =  malloc(bmp_size);
    fprintf(stderr, "Read_jpeg buffer allocated %ld", bmp_size);
    if(buffer == NULL) {
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

    // And free the jpg buffer
    free(jpg_buffer);
    *colorspec = pixel_size;
    *status = 0;
    //Set the output buffer

    for(i = 0; i < *width; i++) {
        for(j = 0; j < *height; j++) {
            if(pixel_size == 1) {
                pixel = (int) buffer[i + j * (*width)];
            } else {
                pixel = ((unsigned long) buffer[j * (*width) * 3 + i * 3 + 0]) << 16 ;
                pixel += ((unsigned long) buffer[j * (*width) * 3 + i * 3 + 1]) << 8;
                pixel += (unsigned long) buffer[j * (*width) * 3 + i * 3 + 2] ;
            }
            (*out_buffer)[i + j * (*width)] = (int) pixel;
        }
    }

    for(i = 0; i < *width; i++) {
        fprintf(stderr, " %d ", out_buffer[0][i]);
    }
}


void write_jpeg(void *imgptr, const char* filename, int* width, int* height, int* quality, int* colorspec,  int*status)
{
    if(strlen(filename) > 0 && strcmp(filename, "-") != 0)
        fprintf(stderr, "cjpg::write_jpeg  opening output file %s", filename);
    fflush(stderr);
    FILE* outfile;
    unsigned long pixel;
    int x, y, depth = 24;
    char* buffer;
    JSAMPROW row_pointer;
    int ** img = imgptr;
    int w, h;
    printf("In write_jpeg in cjpg.c:  %d ", *img);
    //  unsigned char* jbuf;

    int size, numerator;
    size = 256;
    *status = 1;
    {
        for(numerator = 1; numerator <= 8; numerator++) {
            int newsize = (*width) * numerator / 8;
            if(newsize > size)
                break;
        }
    }

    outfile = fopen(filename, "wb");
    if(!outfile) {
        fprintf(stderr, "Failed to open output file %s", filename);
        return ;
    }
    printf("File %s opened ", filename);
    w = *width; h = *height;
    printf("Width  %d   Height %d ", w, h);
    /* collect separate RGB values to a buffer */
    buffer = malloc(sizeof(char) * 3 * w * h);
    if(*colorspec == 1) {
        for(y = 0; y < *height; y++) {
            for(x = 0; x < *width; x++) {
                pixel = (unsigned long) img[x + y * (*width)];
                if(pixel > 0xffffff) {
                    pixel = 0xffffff;
                } else {
                    buffer[y * (*width) * 3 + x * 3 + 0] = (char)(pixel & 0x0000ff);
                    buffer[y * (*width) * 3 + x * 3 + 1] = (char)(pixel & 0x0000ff);
                    buffer[y * (*width) * 3 + x * 3 + 2] = (char)(pixel & 0x0000ff);
                }
            }
        }
    } else {
        for(y = 0; y < *height; y++) {
            for(x = 0; x < *width; x++) {
                pixel = (unsigned long) img[x][y];
                if(pixel > 0xffffff) {
                    pixel = 0xffffff;
                } else {
                    buffer[y * (*width) * 3 + x * 3 + 0] = (char)(pixel >> 16);
                    buffer[y * (*width) * 3 + x * 3 + 1] = (char)((pixel & 0x00ff00) >> 8);
                    buffer[y * (*width) * 3 + x * 3 + 2] = (char)(pixel & 0x0000ff);
                }
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
    jpeg_set_quality(&encInfo, *quality, 1);
    jpeg_start_compress(&encInfo, TRUE);
    while(encInfo.next_scanline < encInfo.image_height) {
        row_pointer = (JSAMPROW) &buffer[encInfo.next_scanline * (depth >> 3) * (*width)];
        jpeg_write_scanlines(&encInfo, &row_pointer, 1);
    }
    jpeg_finish_compress(&encInfo);
    jpeg_destroy_compress(&encInfo);
    free(buffer);
    fclose(outfile);
    *status = 0;
    return ;
}



// memdjpeg - A super simple example of how to decode a jpeg in memory
// Kenneth Finnegan, 2012
// blog.thelifeofkenneth.com
//
// After installing jpeglib, compile with:
// cc memdjpeg.c -ljpeg -o memdjpeg
//
// Run with:
// ./memdjpeg filename.jpg
//
// Version     Date     Time          By
// -------  ----------  -----       ---------
// 0.01     2012-07-09  11:18       Kenneth Finnegan
//

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <syslog.h>
#include <sys/stat.h>

#include <jpeglib.h>


void read_jpeg(const char* file_name, void** img, int*width, int*height, int* colorspec,  int* status)
{
    int fd, rc, i, j;

    char *syslog_prefix =  malloc(1024);
    sprintf(syslog_prefix, "read_jpeg");
    openlog(syslog_prefix, LOG_PERROR | LOG_PID, LOG_USER);

    if(file_name == NULL) {
        fprintf(stderr, "read_jpeg failed\n");
        exit(EXIT_FAILURE);
    }

    // Variables for the decompressor itself
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    // Variables for the output buffer, and how long each row is
    unsigned long bmp_size;
    unsigned char * bmp_buffer;
    int row_stride, w, h, pixel_size;


    // Load the jpeg data from a file into a memory buffer for
    // the purpose of this demonstration.
    // Normally, if it's a file, you'd use jpeg_stdio_src, but just
    // imagine that this was instead being downloaded from the Internet
    // or otherwise not coming from disk
    rc = stat(file_name, &file_info);
    if(rc) {
        syslog(LOG_ERR, "FAILED to stat source jpg");
        exit(EXIT_FAILURE);
    }
    jpg_size = file_info.st_size;
    jpg_buffer = malloc(jpg_size + 100);

    fd = open(file_name, O_RDONLY);
    i = 0;
    while(i < jpg_size) {
        rc = read(fd, jpg_buffer + i, jpg_size - i);
        syslog(LOG_INFO, "Input: Read %d/%lu bytes", rc, jpg_size - i);
        i += rc;
    }
    close(fd);


    syslog(LOG_INFO, "Proc: Create Decompress struct");
    // Allocate a new decompress struct, with the default error handler.
    // The default error handler will exit() on pretty much any issue,
    // so it's likely you'll want to replace it or supplement it with
    // your own.
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);


    syslog(LOG_INFO, "Proc: Set memory buffer as source");
    // Configure this decompressor to read its data from a memory
    // buffer starting at unsigned char *jpg_buffer, which is jpg_size
    // long, and which must contain a complete jpg already.
    //
    // If you need something fancier than this, you must write your
    // own data source manager, which shouldn't be too hard if you know
    // what it is you need it to do. See jpeg-8d/jdatasrc.c for the
    // implementation of the standard jpeg_mem_src and jpeg_stdio_src
    // managers as examples to work from.
    jpeg_mem_src(&cinfo, jpg_buffer, jpg_size);


    syslog(LOG_INFO, "Proc: Read the JPEG header");
    // Have the decompressor scan the jpeg header. This won't populate
    // the cinfo struct output fields, but will indicate if the
    // jpeg is valid.
    rc = jpeg_read_header(&cinfo, TRUE);

    if(rc != 1) {
        syslog(LOG_ERR, "File does not seem to be a normal JPEG");
        exit(EXIT_FAILURE);
    }

    syslog(LOG_INFO, "Proc: Initiate JPEG decompression");
    // By calling jpeg_start_decompress, you populate cinfo
    // and can then allocate your output bitmap buffers for
    // each scanline.
    jpeg_start_decompress(&cinfo);

    w = cinfo.output_width;
    h = cinfo.output_height;
    pixel_size = cinfo.output_components;

    syslog(LOG_INFO, "Proc: Image is %d by %d with %d components",
           w, h, pixel_size);

    bmp_size = w * h * pixel_size;
    bmp_buffer =  malloc(bmp_size);

    // The row_stride is the total number of bytes it takes to store an
    // entire scanline (row).
    row_stride = w * pixel_size;


    syslog(LOG_INFO, "Proc: Start reading scanlines");
    //
    // Now that you have the decompressor entirely configured, it's time
    // to read out all of the scanlines of the jpeg.
    //
    // By default, scanlines will come out in RGBRGBRGB...  order,
    // but this can be changed by setting cinfo.out_color_space
    //
    // jpeg_read_scanlines takes an array of buffers, one for each scanline.
    // Even if you give it a complete set of buffers for the whole image,
    // it will only ever decompress a few lines at a time. For best
    // performance, you should pass it an array with cinfo.rec_outbuf_height
    // scanline buffers. rec_outbuf_height is typically 1, 2, or 4, and
    // at the default high quality decompression setting is always 1.
    while(cinfo.output_scanline < cinfo.output_height) {
        unsigned char *buffer_array[1];
        buffer_array[0] = bmp_buffer + \
                          (cinfo.output_scanline) * row_stride;

        jpeg_read_scanlines(&cinfo, buffer_array, 1);

    }
    syslog(LOG_INFO, "Proc: Done reading scanlines");


    // Once done reading *all* scanlines, release all internal buffers,
    // etc by calling jpeg_finish_decompress. This lets you go back and
    // reuse the same cinfo object with the same settings, if you
    // want to decompress several jpegs in a row.
    //
    // If you didn't read all the scanlines, but want to stop early,
    // you instead need to call jpeg_abort_decompress(&cinfo)
    jpeg_finish_decompress(&cinfo);

    // At this point, optionally go back and either load a new jpg into
    // the jpg_buffer, or define a new jpeg_mem_src, and then start
    // another decompress operation.

    // Once you're really really done, destroy the object to free everything
    jpeg_destroy_decompress(&cinfo);
    // And free the input buffer
    free(jpg_buffer);

    i = h - 1; { //  for (i=0;i<h; i++){
        for(j = 0; j < w; j++) printf("%4.4x ", bmp_buffer[i * w + j]);
        printf("\n");
    }

    *width = w;
    *height = h;

    int (*img_)[h] =  malloc(sizeof(*img_) * w);
    if(!img_) {
        syslog(LOG_ERR, "Image data buffer not allocated");
        exit(EXIT_FAILURE);
    }

    for(i = 0; i < h; i++) for(j = 0; j < w; j++) {
            img_[i][j] = (int) bmp_buffer[i * w + j];
            if(i == h - 1) printf("%d ", img_[i][j]);
        }
    //}else{
    free(bmp_buffer);
    img = ((void*) img_);
    syslog(LOG_INFO, "End of decompression");
    *status = 0;

    // Write the decompressed bitmap out to a ppm file, just to make sure
    // it worked.
    /* fd = open("output.ppm", O_CREAT | O_WRONLY, 0666); */
    /* char buf[1024]; */
    /* rc = sprintf(buf, "P6 %d %d 255\n", *width, *height); */
    /* write(fd, buf, rc); // Write the PPM image header before data */
    /* write(fd, bmp_buffer, bmp_size); // Write out all RGB pixel data */

    /* close(fd); */
}

#define SETUP_JPEG setup_jpeg_
void setup_jpeg(int*width, int*height, int** img)
{
    Jpeg_width = *width;
    Jpeg_height = *height;
    Jpeg_img_buffer = malloc(sizeof(int) * Jpeg_height * Jpeg_width);
    *img = Jpeg_img_buffer;
}
#define DESTROY_JPEG destroy_jpeg_
void destroy_jpeg()
{
    free(Jpeg_img_buffer);
}
