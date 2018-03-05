#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG

#define STBI_READ_JPG stbi_read_jpg_
#define STBI_WRITE_JPG stbi_write_jpg_
#define STBI_WRITE_PNG stbi_write_png_

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "stb_image.h"
#include "stb_image_write.h"


int stbi_read_jpg (const char* file_name, void** data, int* width, int* height, int* components)
{
    FILE *file = fopen(file_name, "rb");
    fseek(file, 0, SEEK_END);
    long size = ftell(file);
    rewind(file);
    unsigned char *buffer = malloc(size);
    fread(buffer, size, 1, file);
    fclose(file);
    size_t len = size;
    bool failure = false;

    stbi_uc *result = stbi_load_from_memory(buffer, len, width, height, components, 0);
failure = failure || result == NULL;
if(failure==false) stbi_image_free((void *)result);
*data = result;


return failure? 1 : 0;
}
