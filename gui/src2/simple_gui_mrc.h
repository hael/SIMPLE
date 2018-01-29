// open include guard

#ifndef __MRC_H_INCLUDED__
#define __MRC_H_INCLUDED__

//=================================


// included external dependencies

#include <string>  
#include <fstream>  

//=================================


// struct declaration

struct MRCFileHeader {
	int 							nx;
	int								ny;
	int								nz;
	int								mode;
	int								nxstart;
	int								nystart;
	int								nzstart;	
	int								mx;
	int								my;
	int								mz;
	int								cellx;
	int								celly;
	int								cellz;
	int								cella;
	int								cellb;
	int								cellg;
	int								mapc;	
	int								mapr;
	int								maps;
	int								dmin;
	int								dmax;
	int								dmean;
	int								ispg;
	int								nsymbt;
	int								extra1;
	int								extra2;
	int								extra3;
	int								extra4;
	int								extra5;
	int								extra6;
	int								extra7;
	int								extra8;
	int								extra9;
	int								extra10;
	int								extra11;
	int								extra12;
	int								extra13;
	int								extra14;
	int								extra15;
	int								extra16;
	int								extra17;
	int								extra18;
	int								extra19;
	int								extra20;
	int								extra21;
	int								extra22;
	int								extra23;
	int								extra24;
	int								extra25;
	int								originx;
	int								originy;
	int								originz;
	char							map[4];
	char							machst[4];
	int								rms;
	int								nlabl;
	char							label1[80];
	char							label2[80];
	char							label3[80];
	char							label4[80];
	char							label5[80];
	char							label6[80];
	char							label7[80];
	char							label8[80];
	char							label9[80];
	char							label10[80];
};

//=================================


// included internal dependencies


//=================================


//  Read header

void readMRCHeader() {
	#include <chrono>
	std::chrono::system_clock::time_point		starttime;
	std::chrono::system_clock::time_point		endtime;
	std::chrono::duration<double> elapsed;
	
	char						buffer[1024];
	struct MRCFileHeader* 		header;
	std::ifstream 				input;
	
	starttime = std::chrono::system_clock::now();
	header = new MRCFileHeader();
	input.open("/tmp/test.mrc", std::ifstream::binary);
	if (input.is_open()){
		input.read ((char*)buffer, 1024); // can i do this in one?
		memcpy(header, &buffer[0], 1024);
		input.close();
	}
	delete header;
	endtime = std::chrono::system_clock::now();
	elapsed = endtime - starttime;
	log("Read header in  " + std::to_string(elapsed.count()) + " seconds");

	
	std::cout << header->nx << std::endl;
/*	std::ifstream 				input;
	int							datastart;
	int							datalength;
	int 						nx;
	int 						ny;
	int 						nz;
	int 						mode;
	int							extendedheader;
	short*						tmpint;
	int 						count;
	float						pixelvalue;
		
	input.open(filename.c_str(), std::ifstream::binary);
	
	if (input.is_open()){
		input.read ((char*)mrcfile->header, 1024);
		memcpy(&nx, &mrcfile->header[0], 4);
		memcpy(&ny, &mrcfile->header[4], 4);
		memcpy(&nz, &mrcfile->header[8], 4);
		memcpy(&mode, &mrcfile->header[12], 4);
		memcpy(&extendedheader, &mrcfile->header[92], 4);
		std::cout << filename << " " << nx << " " << ny << " " << nz << " " << mode << std::endl;
		std::cout << extendedheader << std::endl;

		mrcfile->data = new float[nx * ny];
		if(mode == 2){
			datalength = nx * ny * 4;
			datastart = 1024 + extendedheader + (datalength * frameid);
			input.seekg(datastart, std::fstream::beg);
			input.read((char*)mrcfile->data, datalength);
		} else if (mode == 6){
			tmpint = new short[nx * ny];
			datalength = nx * ny * 2;
			datastart = 1024 + extendedheader + (datalength * frameid);
			input.seekg(datastart, std::fstream::beg);
			input.read((char*)tmpint, datalength);
			for(count = 0; count < (nx * ny); count++){
				mrcfile->data[count] = tmpint[count];
			}			
			delete [] tmpint;
		}
		mrcfile->nx = nx;
		mrcfile->ny = ny;
		mrcfile->nz = nz;
		mrcfile->mode = 2;
		mrcfile->extendedheader = extendedheader;
		
	}
	
	input.close();*/
	
}

//=================================


// close include guard

#endif

//=================================
