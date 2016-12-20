/*
 *   -- SIMPLE addon
 *      Author: Frederic Bonnet, Date: 02nd of May 2016
 *      Monash University
 *      May 2016
 *
 *      Routine which test the conversion of a float towards a string
 *
 *      Non Special case
 * @precisions normal z -> s d c
*/

#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <iostream>

#if defined (CUDA) /* prepropressing for the cuda environment */
#include <cuda.h> // need CUDA_VERSION 7.0
#include <cudnn.h>
#endif

#include "simple.h"
#include "timming.h"

#define NINT 10

/* tester code */
int main(int argc, char *argv[])
{
  int rc = 0; //Return code
  int max_iter = NINT;
  float iint = 123.45;
  int length = 1;
  int indexed_int = 1;
  int size_int_in;
  int num;
  
  length = get_length_of_stringf(&iint);
  indexed_int = get_indexed_int(length);
  printf(ANSI_COLOR_BRIGHT_BLUE"The length of iint: %i\n",length);
  printf(ANSI_COLOR_BRIGHT_GREEN"The length of indexed_int: %i\n",indexed_int);

  size_int_in = (get_length_of_string(&indexed_int)) * sizeof(char);
  char *char_out = (char*)malloc(size_int_in);
 
  for (int i=0; i<NINT ; i++){
    //num = iint + i;
    num = get_indexed_int(length) + iint + i;
    printf(ANSI_COLOR_BRIGHT_BLUE"Inserting an integer[");
    printf(ANSI_COLOR_BRIGHT_YELLOW"%f",iint++);
    printf(ANSI_COLOR_BRIGHT_BLUE"] into a buffer,");
    
    rc = convert_int2char_pos(char_out,&num);
    printf("The buffer looks as: ");
    printf(ANSI_COLOR_BRIGHT_GREEN"%s\n",char_out);
  }

}
