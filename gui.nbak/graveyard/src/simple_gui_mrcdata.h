#ifndef SIMPLE_GUI_MRCDATA_H_INCLUDED
#define SIMPLE_GUI_MRCDATA_H_INCLUDED

#include <iostream>
#include <sys/stat.h>
#include <cstring>

struct MRCFile {
	float *data;
	int miniheader[4]; // MAKE TGHIS WHILE HEADER
};

void readMRCFile(MRCFile *mrcfile, char* filename, int frameid){
	
	std::ifstream 	input;
	int				datalength;
	int				datastart;
	
	input.open(filename, std::ifstream::binary);
	if (input.is_open()){
		input.read ((char*)mrcfile->miniheader, 16);// MAKE TGHIS WHILE HEADER
		datalength = mrcfile->miniheader[0] * mrcfile->miniheader[1];
		datastart = 1024 + (datalength * frameid * 4);
		input.seekg(datastart, std::fstream::beg);
		mrcfile->data = new float[datalength];
		input.read((char*)mrcfile->data, datalength * 4);
	}
	
	input.close();
}

void readMRCFileHeader(MRCFile *mrcfile, char* filename){
	
	std::ifstream 	input;
	int				datalength;
	
	input.open(filename, std::ifstream::binary);
	if (input.is_open()){
		input.read ((char*)mrcfile->miniheader, 16);// MAKE TGHIS WHILE HEADER
	}
	input.close();
}

void selectMRCFile(char* inputfilename,  char* outputfilename, char* selection){
	char*							selectionsplit;
	int								selectionid;
	struct stat						buffer;
	std::ifstream 					infile; 
	std::ofstream 					outfile;
	char*							readbuffer;
	int*							readbufferint;
	int								nx;
	int								ny;
	int								nz;
	int								location;
	
	if(stat(inputfilename, &buffer) == 0 ){
		outfile.open(outputfilename, std::ifstream::binary);
		infile.open(inputfilename, std::ifstream::binary);
		if (infile.is_open() && outfile.is_open()){
			readbuffer = new char [1024];
			infile.read(readbuffer, 1024);
			readbufferint = (int*)readbuffer;
			nx = readbufferint[0];
			ny = readbufferint[1];
			outfile.write(readbuffer, 1024);
			delete [] readbuffer;
			
			selectionsplit = std::strtok(selection, ",");
			readbuffer = new char [nx * ny * 4];
			nz = 0;
			while (selectionsplit != NULL) {
				selectionid = std::stoi(selectionsplit, nullptr);
				location = 1024 + (nx * ny * 4 * selectionid);
				infile.seekg(location, std::fstream::beg);
				infile.read(readbuffer, nx * ny * 4);
				outfile.write(readbuffer, nx * ny * 4);
				nz++;
				selectionsplit = strtok (NULL, ",");
			}
			delete [] readbuffer;
				
			readbuffer = new char [4];
			outfile.seekp(8, std::fstream::beg);
			readbufferint = (int*)readbuffer;
			readbufferint[0] = nz;
			outfile.write(readbuffer, 4);
			delete [] readbuffer;
			
			infile.close();
			outfile.close();
		}
	}
}
#endif
