#include <fstream>
#define NNN 2460

int main(int argc, char **argv)
{   
    ifstream fin;

    // setting MPI environment

    int rank, nprocs;
    MPI_File file;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // reading the initial file

    fin.open("pdf_SP500c001.dat");
    for (int i=0;i<NNN;i++)
    {  
        fin  >> res[i];
        cout << res[i] << endl; // to see, what I have in the file
    }  
    fin.close();

    // starting position in the "res" array as a function of "rank" of process
    int Pstart = (NNN / nprocs) * rank ;
    // specifying Offset for writing to file
    MPI_Offset offset = sizeof(double)*Pstart;
    MPI_File file;
    MPI_Status status;

    // opening one shared file
    MPI_File_open(MPI_COMM_WORLD, "final.txt", MPI_MODE_CREATE|MPI_MODE_WRONLY,
                          MPI_INFO_NULL, &file);

    // setting local for each node array

    double * localArray;
    localArray = new double [NNN/nprocs];

    // Performing some basic manipulation (squaring each element of array)
    for (int i=0;i<(NNN / nprocs);i++)
    {
        localArray[i] = res[Pstart+i]*res[Pstart+i];
    }

    // Writing the result of each local array to the shared final file:

    MPI_File_seek(file, offset, MPI_SEEK_SET);
    MPI_File_write(file, localArray, NNN/nprocs, MPI_DOUBLE, &status);
    MPI_File_close(&file);

    MPI_Finalize();

    return 0;
}
