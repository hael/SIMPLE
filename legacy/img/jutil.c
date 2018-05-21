


#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <jpeglib.h>


namespace jpegutil
{
class JpegEnc
{
public:
    JpegEnc();
    virtual ~JpegEnc();

    JpegEnc(const JpegEnc&) = delete;
    void operator=(const JpegEnc&) = delete;

    /**
     * Sets the color space of the images to be worked with.
     *
     * @param colorSpace     the color space. Can be: JCS_GRAYSCALE, JCS_RGB, JCS_YCbCr, or JCS_CMYK.
     * @param numComponents  the number of components in the color space specified
     */
    void setColorSpace(J_COLOR_SPACE colorSpace, int numComponents);

    /**
     * Gets the current color space.
     *
     * @return  the color space identificator
     */
    J_COLOR_SPACE getColorSpace() const;

    /**
     * Gets the number of components in the current color space.
     *
     * @return  the number of color space components
     */
    int getNumComponents() const;

    /**
     * Sets a compression quality.
     *
     * @param quality  the quality in percents [0, 100]
     */
    void setQuality(int quality);

    /**
     * Gets the current compression quality.
     *
     * @return  the quality value
     */
    int getQuality() const;

    /**
     * Encodes a raw image to a JPEG file on disk using file descriptor.
     *
     * @param img     the pointer to buffer with the raw image
     * @param width   the width of the image
     * @param height  the height of the image
     * @param file    the file descriptor
     * @return        true if success, false if failed
     */
    bool encode(uint8_t* img, int width, int height, FILE* file);

    /**
     * Encodes a raw image to a JPEG file on disk using path.
     *
     * @param img     the pointer to buffer with the raw image
     * @param width   the width of the image
     * @param height  the height of the image
     * @param path    the path to the file
     * @return        true if success, false if failed
     */
    bool encode(uint8_t* img, int width, int height, const char* path);

private:
    bool encode(uint8_t* img, int width, int height);

    jpeg_compress_struct encInfo;
    jpeg_error_mgr errMgr;

    J_COLOR_SPACE colorSpace = JCS_RGB;
    int numComponents = 3;

    int quality = 80;
};

/**
 * JPEG decoder.
 */
class JpegDec
{
public:
    JpegDec();
    virtual ~JpegDec();

    JpegDec(const JpegDec&) = delete;
    void operator=(const JpegDec&) = delete;

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
    J_COLOR_SPACE getColorSpace() const;

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
    bool decode(uint8_t* buffer, size_t len, uint8_t*& img, int& width, int& height);

    /**
     * Reads and decodes a JPEG image through file descriptor.
     *
     * @param[in]  file    the file descriptor
     * @param[out] img     the decoded raw image. The memory is allocated inside the function.
     * @param[out] width   the width of the decoded image
     * @param[out] height  the height of the decoded image
     * @return             true if success, false if failed
     */
    bool decode(FILE* file, uint8_t*& img, int& width, int& height);

    /**
     * Reads and decodes a JPEG file on disk.
     *
     * @param[in]  path    the path to the file
     * @param[out] img     the decoded raw image. The memory is allocated inside the function.
     * @param[out] width   the width of the decoded image
     * @param[out] height  the height of the decoded image
     * @return             true if success, false if failed
     */
    bool decode(const char* path, uint8_t*& img, int& width, int& height);

private:
    bool decode(uint8_t*& img, int& width, int& height);

    jpeg_decompress_struct decInfo;
    jpeg_error_mgr errMgr;

    J_COLOR_SPACE colorSpace = JCS_RGB;
};

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

JpegEnc::JpegEnc()
{
    encInfo.err = jpeg_std_error(&errMgr);
    jpeg_create_compress(&encInfo);
}

JpegEnc::~JpegEnc()
{
    jpeg_destroy_compress(&encInfo);
}

void JpegEnc::setColorSpace(J_COLOR_SPACE colorSpace, int numComponents)
{
    this->colorSpace = colorSpace;
    this->numComponents = numComponents;
}

J_COLOR_SPACE JpegEnc::getColorSpace() const
{
    return colorSpace;
}

int JpegEnc::getNumComponents() const
{
    return numComponents;
}

void JpegEnc::setQuality(int quality)
{
    this->quality = quality;
}

int JpegEnc::getQuality() const
{
    return quality;
}

bool JpegEnc::encode(uint8_t* img, int width, int height, FILE* file)
{
    if(file == nullptr) {
        return false;
    }

    jpeg_stdio_dest(&encInfo, file);
    return encode(img, width, height);
}

bool JpegEnc::encode(uint8_t* img, int width, int height, const char* path)
{
    FILE* file = fopen(path, "wb");
    bool rv = encode(img, width, height, file);
    fclose(file);
    return rv;
}

bool JpegEnc::encode(uint8_t* img, int width, int height)
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

JpegDec::JpegDec()
{
    decInfo.err = jpeg_std_error(&errMgr);
    jpeg_create_decompress(&decInfo);
}

JpegDec::~JpegDec()
{
    jpeg_destroy_decompress(&decInfo);
}

void JpegDec::setColorSpace(J_COLOR_SPACE colorSpace)
{
    this->colorSpace = colorSpace;
}

J_COLOR_SPACE JpegDec::getColorSpace() const
{
    return colorSpace;
}

bool JpegDec::decode(uint8_t* buffer, size_t len, uint8_t*& img, int& width, int& height)
{
    if(buffer == nullptr) {
        return false;
    }

    jpeg_mem_src(&decInfo, buffer, len);
    return decode(img, width, height);
}

bool JpegDec::decode(FILE* file, uint8_t*& img, int& width, int& height)
{
    if(file == nullptr) {
        return false;
    }

    jpeg_stdio_src(&decInfo, file);
    return decode(img, width, height);
}

bool JpegDec::decode(const char* path, uint8_t*& img, int& width, int& height)
{
    FILE* file = fopen(path, "rb");
    bool rv = decode(file, img, width, height);
    fclose(file);
    return rv;
}

bool JpegDec::decode(uint8_t*& img, int& width, int& height)
{
    jpeg_read_header(&decInfo, true);
    decInfo.out_color_space = colorSpace;
    decInfo.raw_data_out = false;

    jpeg_start_decompress(&decInfo);

    width = decInfo.image_width;
    height = decInfo.image_height;
    size_t rowSize = decInfo.image_width * decInfo.num_components;
    img = new uint8_t[decInfo.image_height * rowSize];

    JSAMPROW rows[1];
    size_t offset = 0;
    while(decInfo.output_scanline < decInfo.image_height) {
        rows[0] = img + offset;
        jpeg_read_scanlines(&decInfo, rows, 1);
        offset += rowSize;
    }

    jpeg_finish_decompress(&decInfo);
    return true;
}

class JpegDecPlanarYCbCr
{
public:
    JpegDecPlanarYCbCr();
    virtual ~JpegDecPlanarYCbCr();

    JpegDecPlanarYCbCr(const JpegDecPlanarYCbCr&) = delete;
    void operator=(const JpegDecPlanarYCbCr&) = delete;

    bool decodeI420(FILE* file, uint8_t*& img, int& width, int& height);

private:
    jpeg_decompress_struct decInfo;
    jpeg_error_mgr errMgr;

    JSAMPLE* cache = nullptr;
    size_t cacheSize = 0;
};

JpegDecPlanarYCbCr::JpegDecPlanarYCbCr()
{
    decInfo.err = jpeg_std_error(&errMgr);
    jpeg_create_decompress(&decInfo);
}

JpegDecPlanarYCbCr::~JpegDecPlanarYCbCr()
{
    jpeg_destroy_decompress(&decInfo);
}

bool JpegDecPlanarYCbCr::decodeI420(FILE* file, uint8_t*& img, int& width, int& height)
{
    jpeg_stdio_src(&decInfo, file);

    jpeg_read_header(&decInfo, TRUE);
    decInfo.out_color_space = JCS_YCbCr;
    decInfo.raw_data_out = TRUE;

    jpeg_start_decompress(&decInfo);

    const size_t linesCount = decInfo.max_v_samp_factor * DCTSIZE;

    size_t newCacheSize = linesCount * decInfo.output_width * decInfo.out_color_components;
    if(newCacheSize != cacheSize) {
        delete[] cache;
        cache = new JSAMPLE[newCacheSize];
        cacheSize = newCacheSize;
    }

    JSAMPROW cols[3][linesCount];
    JSAMPLE* ptr = cache;
    for(int c = 0; c < 3; c++) {
        for(size_t i = 0; i < linesCount; i++) {
            cols[c][i] = ptr;
            ptr += decInfo.output_width;
        }
    }

    JSAMPARRAY yuvPtr[3] = {cols[0], cols[1], cols[2]};

    const size_t lumLen = decInfo.output_width * decInfo.output_height;
    const size_t chromWidth = decInfo.output_width / 2;
    const size_t chromHeight = decInfo.output_height / 2;
    const size_t chromLen = chromWidth * chromHeight;
    const size_t dstSize = lumLen + chromLen + chromLen;
    JSAMPLE* dst = new JSAMPLE[dstSize];

    JSAMPLE* y = dst;
    JSAMPLE* u = y + lumLen;
    JSAMPLE* v = u + chromLen;

    for(size_t row = 0; row < decInfo.output_height; row += linesCount) {
        jpeg_read_raw_data(&decInfo, yuvPtr, linesCount);

        // 4:4:4 -> I420

        for(size_t i = 0; i < linesCount; i++) {
            memcpy(y, cols[0][i], decInfo.output_width);
            y += decInfo.output_width;
        }

        for(size_t i = 0; i < linesCount; i += 2) {
            for(size_t j = 0; j != decInfo.output_width; j += 2) {
                *u = cols[1][i][j];
                u += 1;

                *v = cols[2][i][j];
                v += 1;
            }
        }
    }

    img = dst;
    width = decInfo.output_width;
    height = decInfo.output_height;

    jpeg_finish_decompress(&decInfo);
    return true;
}
}

extern "C"
{
    void write_jpeg(uint8_t*& img,  const char* file_name, int& width, int& height, int&quality, int& colorspec, bool* status)
    {
        jpegutil::JpegEnc j; J_COLOR_SPACE jcs;
        *status = false;
        if(file_name == NULL) return;
        //switch(colorspec){
        //  case JCS_GRAYSCALE: j.setColorSpace(JCS_GRAYSCALE,j.getNumComponents());break;
        //  default: // j.setColorSpace(JCS_RGB,j.getNumComponents());
        // || colorspec == JCS_RGB ){ //  , JCS_YCbCr|| JCS_CMYK
        //  j.setColorSpace(colorspec);
        // }
        j.setQuality(quality);
        *status = j.encode(img, width, height, file_name);
    }
    void read_jpeg(const char* file_name, uint8_t*& img, int& width, int& height, int* colorspec,  bool* status)
    {
        jpegutil::JpegDec j; J_COLOR_SPACE jcs;
        *status = false;
        if(file_name == NULL) return;

        *status = j.decode(file_name, img, width, height);
        jcs = j.getColorSpace();
        *colorspec = (int) jpegutil::getNumComponents(jcs);

    }


}
