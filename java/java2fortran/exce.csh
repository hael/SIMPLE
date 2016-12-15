rm *.o *.class

javac JavaCode.java

javah -jni JavaCode

gcc -c -D_REENTRANT -fPIC -I/opt/java/jdk/jdk1.6.0_45/include/ -I/opt/java/jdk/jdk1.6.0_45/include/linux/ -c CCode.c

gfortran -ffree-form -cpp -O3 -fPIC -c FortranCode.f90

gcc -shared CCode.o FortranCode.o -o libmycodeinc.so

LD_LIBRARY_PATH=/home/frederic/SourceCode/PhysSim/Java/Java_to_Fortran

export LD_LIBRARY_PATH


