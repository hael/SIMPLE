#ifndef SIMPLE_GUI_MRCDATA_H_INCLUDED
#define SIMPLE_GUI_MRCDATA_H_INCLUDED

struct MRCFile {
	float *data;
	int miniheader[4]; // MAKE TGHIS WHILE HEADER
};

void readMRCFile(MRCFile *mrcfile, std::string filename){
	
	std::ifstream 	input;
	int				datalength;
	
	input.open(filename.c_str(), std::ifstream::binary);
	if (input.is_open()){
		input.read ((char*)mrcfile->miniheader, 16);// MAKE TGHIS WHILE HEADER
		input.seekg(1024, std::fstream::beg);
		datalength = mrcfile->miniheader[0] * mrcfile->miniheader[1];
		mrcfile->data = new float[datalength];
		input.read((char*)mrcfile->data, datalength * 4);
	}
	
	input.close();
}



#endif
