//               Example of C Code with Conditional Compilation
//               ==============================================
//
// This illustrates an example of C code that can either be called from Java,
// or for debugging purposes is the main C code not called from Java.  In
// both cases the code calls code in FORTRAN.  The ending of each line has the
// following meanings:
// 1) "// J" - this line is specific to being called from Java,
// 2) "// c" - this line is specific to the main C program when not called
//             from Java,
// 3) Blank  - this line is common to both.

#define JAVA 1 // This flag is set when the code is called from Java.
#include <stdio.h>

// FORTRAN routines have to be prototyped as extern, and parameters are
// passed by reference.  Note also that for g77 the function name in C by
// default must have "_" as a suffix.

extern int sumsquaredf_(int *, int []);

#ifdef JAVA // ---- This block of code is compiled when called from Java ------
                                                                               // J
#include "JavaCode.h"  // Required header for JNI                              // J
                                                                               // J
// When calling C code from Java, main() must be replaced by a declaration     // J
// similar to below, where the function name is given by "Java_" + the name    // J
// of the class in the Java code that calls this C code, in this case          // J
// "JavaCode", + "_" + the name of this C function called from Java, in this   // J
// case "sumsquaredc".  This is followed by at least two parameters as below,  // J
// plus possibly more if they are required.                                    // J
                                                                               // J
JNIEXPORT jint JNICALL Java_JavaCode_sumsquaredc(JNIEnv *env,                  // J
                       jobject obj, jintArray ja) {                            // J
                                                                               // J
// Data from any additional parameters are passed via special pointers as      // J
// shown here.                                                                 // J
                                                                               // J
    jsize n = (*env)->GetArrayLength(env, ja);                                 // J
    jint *a = (*env)->GetIntArrayElements(env, ja, 0);                         // J
    int i,result;                                                              // J
                                                                               // J
#else // ---- This block of code is compiled when Java is not used ------------
                                                                               // C
#define LEN 10                                                                 // C
int main() {                                                                   // C
    int i,n,result,a[LEN];                                                     // C
    n = LEN;         // Do some initializations when Java is not used.         // C
    for (i=0; i<n; i++)                                                        // C
        a[i]=i;                                                                // C
                                                                               // C
#endif // ---- This block of code is common to C and Java ---------------------

    printf("-- We are now in the C program CCode --\n");
    printf("Print contents of array a[] copied from arr[] in Java\n");
    for (i=0; i<n; i++)
        printf("%2d %5d\n",i,a[i]); 

    printf("Call the FORTRAN code\n");

//  The data are passed to the FORTRAN program, then the results are
//  returned in a way similar to this.

    result = sumsquaredf_(&n,a);
    printf("-- We have now returned back to the C program --\n");

    printf("Print contents of array a[]\n");
    for (i=0; i<n; i++)
        printf("%2d %5d\n",i,a[i]);

    printf("Sum of squares in array = %d\n",result);

#ifdef JAVA // ---- This block of code is compiled when called from Java ------
                                                                               // J
//  Instead of ending as a normal C program, the pointers must be              // J
//  cleared before returning to Java.                                          // J
                                                                               // J
    (*env)->ReleaseIntArrayElements(env, ja, a, 0);                            // J
    return result;                                                             // J
                                                                               // J
#else // ---- Otherwise return as a normal C program --------------------------
                                                                               // C
    return;                                                                    // C
                                                                               // C
#endif
// ---- Exit C code then either return to Java or finish ----------------------
}
