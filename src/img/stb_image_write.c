
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG

#define STBI_READ_JPG stbi_read_jpg_
#define STBI_WRITE_JPG stbi_write_jpg_
#define STBI_WRITE_PNG stbi_write_png_

//#include <stdbool.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <errno.h>
#include "stb_image.h"
#include "stb_image_write.h"
//    #include <sys/types.h>
//       #include <sys/stat.h>
//       #include <unistd.h>
//      #include <fcntl.h>


int stbi_read_jpg(const char* file_name, void** result, int* width, int* height, int* components)
{
  fprintf(stderr, " In stbi_read_jpg -- reading 1 %s\n",file_name);fflush(stderr);
  //stbi_uc
  *result = stbi_load(file_name, width, height, components, 3);
  if(*result ==NULL) { fprintf(stderr, " stbi_load result invalid \n");return 0;}
  if(0){
   fprintf(stderr, " stbi_load result valid \n");
  fprintf(stderr, " stbi_load result valid \n");
  fprintf(stderr, " stbi_load result size %ld\n", sizeof(*result)/sizeof(stbi_uc) );
  fprintf(stderr, " stbi_load result[0] %4.4x\n", ((stbi_uc)*result));
  fprintf(stderr, " stbi_load result[rowend] %4.4x\n", ((stbi_uc) *result +(*width) -1));
  fprintf(stderr, " stbi_load result[end] %4.4x\n", ((stbi_uc)(*result) +(*width)*(*height)-1));
  }
  //*data = malloc(sizeof(unsigned char) * (*width) * (*height) * 3);
  //data = &result;
  //memcpy(*data, result, (*width) * (*height) * 3);
   //stbi_image_free((void *)result);
  return 1;
}


/*

  //fprintf(stderr, " In stbi_read_jpg -- reading 1 %s\n",file_name);fflush(stderr);
    bool failure = false;
    struct stat buf;
    if(stat(file_name, &buf) != 0 ) {
      fprintf(stderr," fstat returned non-zero\n");
      return 1;
    }
    fprintf(stderr, " In stbi_read_jpg --  total size in bytes %d\n", buf.st_size);
    FILE *file = NULL;
    file = fopen(file_name, "rb");
    if(!file){
       fprintf(stderr, " In stbi_read_jpg -- file_name not found %s\n",file_name);//fflush(stderr);
       return 1;
     }
      fprintf(stderr, " In stbi_read_jpg -- file found %s  %ld\n",file_name, *file);//fflush(stderr);
      unsigned char *buffer = malloc( buf.st_size);
     int fd = open( file_name, O_RDONLY);
     ssize_t ret = read(fd, buffer, buf.st_size );
     if(ret == -1)  {
       fprintf(stderr," read returned error \n");
       free(buffer); close(fd);
       return 1;
     }
     close(fd);
     if( fseek(file, 0L, SEEK_END) != 0 ){
       fprintf(stderr, " In stbi_read_jpg -- fseek failed\n");
       perror(" In stbi_read_jpg -- fseek failed\n");
       fclose(file);
       return 1;
     }
     fprintf(stderr, " In stbi_read_jpg -- reading ftell  %s\n",file_name);fflush(stderr);
     long size = ftell(file);
     fprintf(stderr, " In stbi_read_jpg -- rewind  %s\n",file_name);fflush(stderr);
    rewind(file);
    fprintf(stderr, " In stbi_read_jpg -- reading 6  %s\n",file_name);fflush(stderr);
    unsigned char *buffer = malloc( sizeof(unsigned char) * size);
    if(buffer == NULL) fprintf(stderr, " buffer invalid \n");
    size_t bytes_read = fread(buffer, size, 1, file);
     fclose(file);
     fprintf(stderr, " In %s line %d -- bytes read  %d\t%d\n", __FILE__, __LINE__, bytes_read, size );fflush(stderr);

      size_t len = buf.st_size;
     fprintf(stderr, " In stbi_read_jpg -- calling stbi_load_from_memory\n");
     stbi_uc *result = stbi_load_from_memory(buffer, len, width, height, components, 0);
     free(buffer);



    //    stbi_uc *result = stbi_load(file_name, width, height, components, 3);
     //    stbi_uc *result = stbi_load_from_file(file, width, height, components, 3);
    if(!result) fprintf(stderr, " result invalid \n");
    //fclose(file);

    failure = failure || result == NULL;
    if(failure == false) {
      *data = malloc(sizeof(unsigned char) * (*width) * (*height) * 3);
      //fprintf(stderr, " size of result sizeof(result)/sizeof(stbi_uc) %d\n", sizeof(result)/sizeof(stbi_uc));
      //fprintf(stderr, " size of data sizeof(data)/sizeof(char) %d\n", sizeof(data)/sizeof(char));
      //for (int i=0;i<sizeof(result)/sizeof(stbi_uc); i++) (*data)[i] = (unsigned char) result[i];
      memcpy(*data, result, sizeof(result));
      stbi_image_free((void *)result);
    }
    return failure ? 1 : 0;
}
*/
