#ifndef SIMPLE_GUI_VIEW_H_INCLUDED
#define SIMPLE_GUI_VIEW_H_INCLUDED

#include "simple_gui_mrcdata.h"


struct JPEGData {
	unsigned char 	*pixels;
	unsigned char 	*jpeg;
	unsigned long	jpegsize;
};


void viewStarFile(StarFile *starfile, std::string columnname, int *columnnumber){
	std::vector<std::string>::iterator 	it;
	int 								columncount = 0;
	
	for(it = starfile->columns.begin() ; it != starfile->columns.end(); ++it){
		if(*it == columnname){
			*columnnumber = columncount;
			return;
		}
		columncount++;
	}
}

void encodeJPEG(MRCFile *mrcfile, JPEGData *jpegdata){
	unsigned long 					jpegsize = 0;
	struct jpeg_compress_struct 	cinfo;
	struct jpeg_error_mgr 			jerr;
	
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_mem_dest(&cinfo, &jpegdata->jpeg, &jpegdata->jpegsize);
	cinfo.image_width = mrcfile->miniheader[0];
	cinfo.image_height = mrcfile->miniheader[1];
	cinfo.input_components = 1;
	cinfo.in_color_space = JCS_GRAYSCALE;
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, 50, TRUE);
	jpeg_start_compress(&cinfo, TRUE);
	JSAMPROW row_pointer;
	for(int i = 0; i < mrcfile->miniheader[1]; i++){
		row_pointer = (JSAMPROW) &jpegdata->pixels[i * mrcfile->miniheader[0]];
		jpeg_write_scanlines(&cinfo, &row_pointer, 1);
	}
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);	
}

void processImage(MRCFile *mrcfile, JPEGData *jpegdata, float contrast, float brightness){
	
	float 			datamean = 0;
	int 			datalength;
	float 			datadiff2 = 0;
	int 			i;
	float 			datasd;
	float 			datanorm;
	
	datalength = mrcfile->miniheader[0] * mrcfile->miniheader[1];

	for(i = 0; i < datalength; i++){
		datamean += mrcfile->data[i];
	}
	
	datamean = datamean / datalength;

	for(i = 0; i < datalength; i++){
		datadiff2 += pow((mrcfile->data[i] - datamean), 2);
	}
	
	datasd = sqrt(datadiff2 /datalength);
	
	jpegdata->pixels = new unsigned char[datalength];

	for(i = 0; i < datalength; i++){
		datanorm = (mrcfile->data[i] - datamean) / datasd;
		datanorm *= 128/contrast; //contrast
		datanorm += brightness; //brightness
		if(datanorm > 254){
			jpegdata->pixels[i] = (int) 254;
		} else if(datanorm < 0){
			jpegdata->pixels[i] = (int) 0;
		} else {
			jpegdata->pixels[i] = round(datanorm);
		}
	}

}

void returnJPEG(std::string &messagedata, JPEGData *jpegdata){
	std::string				contrast;
	std::string				brightness;
	std::string				filename;
	int 					stringcount = 0;
	
	while (stringcount < messagedata.size()) {
		if (messagedata.at(stringcount) == '&'){
			stringcount++;
			break;
		}
		stringcount++;
	}
	
	messagedata.erase(0, stringcount);
	
	stringcount = 0;
	while (stringcount < messagedata.size()) {
		if (messagedata.at(stringcount) == '&'){
			stringcount++;
			break;
		}
		stringcount++;
	}

	contrast = messagedata.substr(0, stringcount);
	messagedata.erase(0, stringcount);

	stringcount = 0;
	while (stringcount < messagedata.size()) {
		if (messagedata.at(stringcount) == '&'){
			stringcount++;
			break;
		}
		stringcount++;
	}
	
	brightness = messagedata.substr(0, stringcount);
	messagedata.erase(0, stringcount);
	
	filename = messagedata;
	
	MRCFile* mrcfile = new MRCFile();
		
	readMRCFile(mrcfile, "/tmp/TEST/59_unblur_ctffind/2_intg.mrc");
	processImage(mrcfile, jpegdata, strtof(contrast.c_str(), NULL), strtof(brightness.c_str(), NULL));
	encodeJPEG(mrcfile, jpegdata);
	std::cout << "JPEG READY" << std::endl;
	delete mrcfile;
}

void starFileViewer(struct mg_connection *nc, std::string &messagedata, std::stringstream &returnstring){
	
	std::string							filename;
	std::vector<std::string>::iterator 	it;

	StarFile *starfile = new StarFile();
	getKeyValFromString(&filename, messagedata, "filename");
	
	readStarFile(starfile, filename);
	
	returnstring << "starcolumns=";
	for(it = starfile->columns.begin() ; it != starfile->columns.end(); ++it){
		returnstring << *it << ";";
	}
	
	delete starfile;
	
}

void fileViewer(struct mg_connection *nc, std::string &messagedata){

	std::string			filename;
	std::stringstream 	returnstring;
	
	getKeyValFromString(&filename, messagedata, "filename");
	
	if(filename.find(".star") != std::string::npos){
		starFileViewer(nc, messagedata, returnstring);
	}
	
}


#endif
