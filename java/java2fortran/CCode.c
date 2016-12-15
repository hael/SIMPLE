// CCode with documentation

#include <stdio.h>
#include "JavaCode.h"  // Required header for JNI

// FORTRAN routines have to be prototyped as extern, and parameters are
// passed by reference.  Note also that for g77 the function name in C by
// default must be prefixed by a "_".

extern int sumsquaredf_(int *, int []);

// When calling C code from Java, main() must be replaced by a declaration
// similar to below, where the function name is given by "Java_" + the name
// of the class in the Java code that calls this C code, in this case
// "JavaCode", + "_" + the name of this C function called from Java, in this
// case "sumsquaredc".  This is followed by at least two parameters as below,
// plus possibly more if more are required.

JNIEXPORT jint JNICALL Java_JavaCode_sumsquaredc(JNIEnv *env,
                       jobject obj, jintArray ja) {

//  Data from any additional parameters are passed via special pointers as
//  shown here.

    jsize n = (*env)->GetArrayLength(env, ja);
    jint *a = (*env)->GetIntArrayElements(env, ja, 0);
    int i,result;

    printf("-- We are now in the C program CCode --\n");
    printf("Print contents of array a[] copied from arr[] in Java\n");
    for (i=0; i<n; i++)
        printf("%2d %5d\n",i,a[i]);

    printf("Call the FORTRAN code\n");

//  The data are passed to the FORTRAN program, then the results are returned
//  in a way similar to this.

    result = sumsquaredf_(&n,a);
    printf("-- We have now returned back to the C program --\n");

    printf("Print contents of array a[]\n");
    for (i=0; i<n; i++)
       printf("%2d %5d\n",i,a[i]);

    printf("Sum of squares in array = %d\n",result);

//  Instead of ending as a normal C program, the pointers must be cleared
//  before returning to Java.

    (*env)->ReleaseIntArrayElements(env, ja, a, 0);
    return result;
}
