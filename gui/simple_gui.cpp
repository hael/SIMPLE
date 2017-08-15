#define KEY(X) #X
#define VALUE(X) KEY(X)
#define SIMPLE_DIR_VALUE VALUE(SIMPLE_DIR)

#include <stdio.h>
#include <jpeglib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <sqlite3.h> 
#include <sys/types.h>
#include <pwd.h>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <ctime>
#include <vector>
#include <signal.h>

extern "C" {
    #include "mongoose.h"
}

static const char *s_http_port = "8084";
static struct mg_serve_http_opts s_http_server_opts;
static sig_atomic_t s_signal_received = 0;

struct StarFile {
	std::vector<std::string> columns;
	std::vector<std::string> data;
};

struct UniDoc {
	std::vector<std::string> data;
};

static void getKeyValFromString(std::string *output, std::string argstring, std::string key){ 				
	std::size_t location;
	char character;
	output->clear();
	key.append("=");
	location = argstring.find(key) + key.length();
	if(location != key.length() - 1){
		while ( location < argstring.length() ){
			character = argstring.c_str()[location];
			if (character == '^'){
				*output += " ";
			} else if (character != ' '){
				*output += character;
			} else {
				break;
			}
			location = location + 1;
		}
	}
}

static void getElementFromString(std::string *output, std::string argstring, int element){ 
	std::size_t location;
	char character;
	int currentelement;
	
	output->clear();
	location = 0;
	currentelement = 0;
	
	while(argstring.c_str()[location] == ' ' && location < argstring.length()){
		location = location + 1;
	}
	
	while ( location < argstring.length() ){
		character = argstring.c_str()[location];
		if ( character != ' '){
			if(currentelement == element){
				output->push_back(character);
			}
		} else {
			currentelement = currentelement + 1;
		}
		location = location + 1;
	}
}

static void writeStarFile(StarFile *starfile, std::string filename){
	std::ofstream output;
	std::vector<std::string>::iterator it;
	
	output.open(filename.c_str());
	if(output.is_open()){
		// Write data_
		output << "data_" << "\n";
		// Write loop_
		output << "loop_" << "\n";
		// Write column headers
		for (it = starfile->columns.begin() ; it != starfile->columns.end(); ++it){
			output << "_" << *it << "\n";
		}
		// Write data
		for (it = starfile->data.begin() ; it != starfile->data.end(); ++it){
			output << *it << "\n";
		}
		output.close();
	} else {
		
	}
	
}

static void readStarFile(StarFile *starfile, std::string filename){
	std::ifstream input;
	std::string line;
	std::cout << filename << std::endl;
	input.open(filename.c_str());
	if(input.is_open()){
		while(getline(input, line)){
			if(line.at(0) == '_'){
				line.erase(0,1);
				starfile->columns.push_back(line);
				std::cout << line << std::endl;
			} else if(line.at(0) == '#' || line.at(0) == ' '){
			} else if(line.find("data_") != std::string::npos || line.find("loop_") != std::string::npos) {
			} else {
				starfile->data.push_back(line);
				std::cout << line << std::endl;
			}
		}
		input.close();
	} else {
		std::cout << "no starfil" << std::endl;
	}
	
}

static void getStarColumn(StarFile *starfile, std::string columnname, int *columnnumber){
	std::vector<std::string>::iterator it;
	int columncount = 0;
	
	for(it = starfile->columns.begin() ; it != starfile->columns.end(); ++it){
		if(*it == columnname){
			*columnnumber = columncount;
			return;
		}
		columncount++;
	}
}

static void readUniDoc(UniDoc *unidoc, std::string filename){
	std::ifstream input;
	std::string line;
	input.open(filename.c_str());
	if(input.is_open()){
		while(getline(input, line)){
			unidoc->data.push_back(line);
		}
		input.close();
	} else {
		std::cout << "no unidoc" << std::endl;
	}
	
}

static void getUser(std::string *user, std::string messagedata){ 							
	struct passwd* pw;
	register uid_t uid;

	getKeyValFromString(user, messagedata, "usr");
	
	if(*user == ""){
		uid = geteuid();
		pw = getpwuid(uid);
		*user = pw->pw_name;
	}
}

static void getUserHome(std::string *userhome, std::string user){ 						
	struct passwd* pw;
	
	pw = getpwnam(user.c_str());
	*userhome = pw->pw_dir;
}

static void returnWSMessage(struct mg_connection *nc, std::string errorstring){
	std::string returnstring = "message=";
	for(unsigned int i = 0; i < errorstring.length(); i++){
		if(errorstring.c_str()[i] == ' '){
			returnstring += "^";
		}else{
			returnstring += errorstring.c_str()[i];
		}
	}
	mg_printf_websocket_frame(nc, WEBSOCKET_OP_TEXT, returnstring.c_str(), returnstring.length());
}

static void processImage(float *data, unsigned char *pixels, int datalength, int stackcount, float contrast, float brightness){
	float datamean = 0;
	for(int i = 0; i < datalength; i++){
		datamean += data[i];
	}
	datamean = datamean / datalength;
	
	float datadiff2 = 0;
	for(int i = 0; i < datalength; i++){
		datadiff2 += pow((data[i] - datamean),2);
	}
	
	float datasd = sqrt(datadiff2 /datalength);

	float datanorm;
	
	for(int i = 0; i < datalength; i++){
		datanorm = (data[i] - datamean) / datasd;
		datanorm *= 128/contrast; //contrast
		datanorm += brightness; //brightness
		if(datanorm > 254){
			pixels[i + (stackcount * datalength)] = (int) 254;
		} else if(datanorm < 0){
			pixels[i + (stackcount * datalength)] = (int) 0;
		} else {
			pixels[i + (stackcount * datalength)] = round(datanorm);
		}
	}
}

static void binaryImage(float *data, unsigned char *pixels, int nx, int ny){
	// Messing with binarization
	int backgroundtilesize = 50; 
	int tilex = 0;
	int tiley = 0;
	int x,y;
	float datamean = 0;
	float datasd = 0;
	float datanorm = 0;
	
	int *model = new int[nx * ny];
	
	while(tilex < (nx - backgroundtilesize)){
		tiley = 0;
		while(tiley < (ny - backgroundtilesize)){
				datamean = 0;
				datasd = 0;
				datanorm = 0;
				for(x = tilex; x < (tilex + backgroundtilesize); x++){
					for(y = tiley; y < (tiley + backgroundtilesize); y++){
						datamean += data[(y * nx) + x];
					}
				}
				datamean = datamean / (backgroundtilesize * backgroundtilesize);
				
				for(x = tilex; x < (tilex + backgroundtilesize); x++){
					for(y = tiley; y < (tiley + backgroundtilesize); y++){
						datasd += pow((data[(y * nx) + x] - datamean),2);
					}
				}
				datasd = sqrt(datasd / (backgroundtilesize * backgroundtilesize));
				
				for(x = tilex; x < (tilex + backgroundtilesize); x++){
					for(y = tiley; y < (tiley + backgroundtilesize); y++){
						
						datanorm = (data[(y * nx) + x] - datamean) / datasd;
						data[(y * nx) + x] = datanorm;
						
						//if(datanorm > 1){
						//	model[(y * nx) + x] = 1;
						//} else {
						//	model[(y * nx) + x] = 0;
						//}
						
					/*	datanorm *= 12.8;
						datanorm += 128;
						if(round(datanorm) > 128){
							data[(y * nx) + x] = 255;
						} else if(round(datanorm) < 0){
							data[(y * nx) + x] = 0;
						} else {	
							data[(y * nx) + x] =  round(datanorm);
						}*/
					//	} else if(datanorm < -1){
							//pixels[(y * nx) + x] = (int) 0;
						//
							//pixels[(y * nx) + x] = (int)0;
						//	data[(y * nx) + x] = ;
						//}
					}
				}
			tiley += backgroundtilesize;
		}
		tilex += backgroundtilesize;
	}
	
	tilex = 1;
	
	while(tilex < (nx - 1)){
		tiley = 1;
		while(tiley < (ny - 1)){
			datanorm = 0;
			datanorm += data[((tiley - 1) * nx) + (tilex - 1)];
			datanorm += data[((tiley - 1) * nx) + tilex];
			datanorm += data[((tiley - 1) * nx) + (tilex +1)];
			datanorm += data[(tiley * nx) + (tilex - 1)];
			datanorm += data[(tiley * nx) + tilex];
			datanorm += data[(tiley * nx) + (tilex + 1)];
			datanorm += data[((tiley + 1) * nx) + (tilex - 1)];
			datanorm += data[((tiley + 1) * nx) + tilex];
			datanorm += data[((tiley + 1) * nx) + (tilex +1)];
			datanorm /= 9;
			data[(tiley * nx) + tilex] = datanorm;
						
			tiley++;
		}
		tilex++;
	}
	
	tilex = 1;
	while(tilex < (nx - 1)){
		tiley = 1;
		while(tiley < (ny - 1)){
			datanorm = 0;
			datanorm += data[((tiley - 2) * nx) + (tilex - 2)] / 32;
			datanorm += data[((tiley - 2) * nx) + (tilex - 1)] / 32;;
			datanorm += data[((tiley - 2) * nx) + tilex] / 32;;
			datanorm += data[((tiley - 2) * nx) + (tilex +1)] / 32;;
			datanorm += data[((tiley - 2) * nx) + (tilex +2)] / 32;;
			datanorm += data[((tiley - 1) * nx) + (tilex - 2)] / 32;;
			datanorm += data[((tiley - 1) * nx) + (tilex - 1)] / 18;;
			datanorm += data[((tiley - 1) * nx) + tilex] / 18;;
			datanorm += data[((tiley - 1) * nx) + (tilex +1)] / 18;;
			datanorm += data[((tiley - 1) * nx) + (tilex +2)] / 32;;
			datanorm += data[(tiley * nx) + (tilex - 2)] / 32;;
			datanorm += data[(tiley * nx) + (tilex - 1)] / 18;;
			datanorm += data[(tiley * nx) + tilex] / 18;;
			datanorm += data[(tiley * nx) + (tilex + 1)] / 18;;
			datanorm += data[(tiley * nx) + (tilex + 2)] / 32;;
			datanorm += data[((tiley + 1) * nx) + (tilex - 2)] / 32;;
			datanorm += data[((tiley + 1) * nx) + (tilex - 1)] / 18;;
			datanorm += data[((tiley + 1) * nx) + tilex] / 18;;
			datanorm += data[((tiley + 1) * nx) + (tilex +1)] / 18;;
			datanorm += data[((tiley + 1) * nx) + (tilex +2)] / 32;;
			datanorm += data[((tiley + 2) * nx) + (tilex - 2)] / 32;;
			datanorm += data[((tiley + 2) * nx) + (tilex - 1)] / 32;;
			datanorm += data[((tiley + 2) * nx) + tilex] / 32;;
			datanorm += data[((tiley + 2) * nx) + (tilex +1)] / 32;;
			datanorm += data[((tiley + 2) * nx) + (tilex +2)] / 32;;
			
			if(datanorm < 0){
				pixels[(tiley * nx) + tilex] = (int) 0;
			} else if(datanorm > 0){
				pixels[(tiley * nx) + tilex] = (int) 255;	
			}
						
			tiley++;
		}
		tilex++;
	}
	
	
/*	
	while(tilex < (nx - 1)){
		tiley = 1;
		while(tiley < (ny - 1)){
			datanorm = 0;
			datanorm += data[((tiley - 1) * nx) + (tilex - 1)];
			datanorm += data[((tiley - 1) * nx) + tilex];
			datanorm += data[((tiley - 1) * nx) + (tilex +1)];
			datanorm += data[(tiley * nx) + (tilex - 1)];
			datanorm += data[(tiley * nx) + tilex];
			datanorm += data[(tiley * nx) + (tilex + 1)];
			datanorm += data[((tiley + 1) * nx) + (tilex - 1)];
			datanorm += data[((tiley + 1) * nx) + tilex];
			datanorm += data[((tiley + 1) * nx) + (tilex +1)];
			datanorm /= 9;
			

			if(round(datanorm) > 200){
				pixels[(tiley * nx) + tilex] = 255;
			} else if(round(datanorm) < 50){
				pixels[(tiley * nx) + tilex] = 0;
			} else {
				pixels[(tiley * nx) + tilex] = round(datanorm);
			}
						
			tiley++;
		}
		tilex++;
	}	

	*/
	/*
	while(tilex < (nx - neigbourtilesize)){
		tiley = 0;
		while(tiley < (ny - neigbourtilesize)){
				for(x = tilex; x < (tilex + neigbourtilesize); x++){
					for(y = tiley; y < (tiley + neigbourtilesize); y++){
						datamean += data[(y * nx) + x];
					}
				}
				datamean = datamean / (neigbourtilesize * neigbourtilesize);
				
				for(x = tilex; x < (tilex + neigbourtilesize); x++){
					for(y = tiley; y < (tiley + neigbourtilesize); y++){
						
						if(datamean > 200){
							//pixels[(y * nx) + x] = (int) 255;
							data[(y * nx) + x] = 255;
					//	} else if(datanorm < -1){
							//pixels[(y * nx) + x] = (int) 0;
						} else {
							//pixels[(y * nx) + x] = (int)0;
							data[(y * nx) + x] = 0;
						}
					}
				}
			tiley += neigbourtilesize;
		}
		tilex += neigbourtilesize;
	}
	
	neigbourtilesize = 3; 
	tilex = 0;
	
	while(tilex < (nx - neigbourtilesize)){
		tiley = 0;
		while(tiley < (ny - neigbourtilesize)){
				for(x = tilex; x < (tilex + neigbourtilesize); x++){
					for(y = tiley; y < (tiley + neigbourtilesize); y++){
						datamean += data[(y * nx) + x];
					}
				}
				datamean = datamean / (neigbourtilesize * neigbourtilesize);
				
				for(x = tilex; x < (tilex + neigbourtilesize); x++){
					for(y = tiley; y < (tiley + neigbourtilesize); y++){
						
						if(datamean > 200){
							pixels[(y * nx) + x] = (int) 255;
							//data[(y * nx) + x] = 255;
					//	} else if(datanorm < -1){
							//pixels[(y * nx) + x] = (int) 0;
						} else {
							pixels[(y * nx) + x] = (int)0;
							//data[(y * nx) + x] = 0;
						}
					}
				}
			tiley += neigbourtilesize;
		}
		tilex += neigbourtilesize;
	}
	*/
	delete[] model;
}

static void sendJPEG(struct mg_connection *nc, unsigned char *header, unsigned char *pixels, int nx, int ny, int framecount){
	unsigned char *jpeg = NULL;
	unsigned long jpegsize = 0;
	unsigned long messagesize = 560;
	unsigned char *message;
	
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_mem_dest(&cinfo, &jpeg, &jpegsize);
	cinfo.image_width = nx;
	cinfo.image_height = ny * framecount;
	cinfo.input_components = 1;
	cinfo.in_color_space = JCS_GRAYSCALE;
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, 50, TRUE);
	jpeg_start_compress(&cinfo, TRUE);
	JSAMPROW row_pointer;
	for(int i = 0; i < ny * framecount; i++){
		row_pointer = (JSAMPROW) &pixels[i * nx];
		jpeg_write_scanlines(&cinfo, &row_pointer, 1);
	}
	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
	
	
	messagesize += jpegsize;
	message = new unsigned char[messagesize];
	
	memcpy(&message[0], &header[0], 560);
	memcpy(&message[560], &jpeg[0], jpegsize);
	mg_send_websocket_frame(nc, WEBSOCKET_OP_BINARY, &message[0], messagesize);
	delete[] message;
	free(jpeg);
	
}

static void createHeader(unsigned char *header, std::string headertext){
	unsigned int i;
	for(i = 0; i < headertext.length(); i++){
		header[i] = headertext[i];
	}
	while(i < 560){
		header[i] = ' ';
		i++;
	}
}

static void getMRCFile(struct mg_connection *nc, std::string filename, std::string callback, float contrast, float brightness){
	std::ifstream mrcfile;
	int attributes[4];
	std::stringstream returnstringss;
	std::string returnstring;
	std::stringstream headertext;
	int datalength;
	
	unsigned char *pixels;
	unsigned char *header;
	float *data;
	header = new unsigned char[560];
	
	mrcfile.open(filename.c_str(), std::ifstream::binary);
	
	if (mrcfile.is_open()){
		
		mrcfile.read ((char*)attributes, 16);
		mrcfile.seekg(1024, std::fstream::beg);
	
		// Case when MRC mode 2 and single image
		if(attributes[3] == 2 && attributes[2] == 1){
			datalength = attributes[0] * attributes[1];
			data = new float[datalength];
			pixels = new unsigned char[datalength];
			mrcfile.read((char*)data, datalength * 4);
			processImage(data, pixels, datalength, 0, contrast, brightness);
			delete [] data;
			headertext << "callback=" << callback << " ";
			headertext << "filename=" << filename << " ";
			headertext << "nx=" << attributes[0] << " ";
			headertext << "ny=" << attributes[1] << " ";
			headertext << "nz=" << attributes[2] << " ";
			createHeader(header,headertext.str());
			sendJPEG(nc, header, pixels, attributes[0], attributes[1], 1);
			delete [] header;
			delete [] pixels;
			
		} else if(attributes[3] == 2 && attributes[2] > 1){	
			datalength = attributes[0] * attributes[1];
			data = new float[datalength];
			int maxframes = floor(65000 / attributes[2]);
			
			if(attributes[2] < maxframes){
				maxframes = attributes[2];
			}

			pixels = new unsigned char[datalength * maxframes];
			
			for(int stackcount = 0; stackcount < maxframes; stackcount++){
				mrcfile.read((char*)data, datalength * 4);
				processImage(data, pixels, datalength, stackcount, contrast, brightness);
			}
			delete [] data;
			headertext << "callback=" << callback << " ";
			headertext << "nx=" << attributes[0] << " ";
			headertext << "ny=" << attributes[1] << " ";
			headertext << "nz=" << maxframes << " ";
			
			createHeader(header,headertext.str());
			sendJPEG(nc, header, pixels, attributes[0], attributes[1], maxframes);
			delete [] header;
			delete [] pixels;
		}
	}
	mrcfile.close();
}

static void getTextFile(struct mg_connection *nc, std::string filename, std::string callback){
	std::ifstream textfile;
	std::string line, returntextstring;
	std::stringstream returntext;
	
	textfile.open(filename.c_str());
	
	returntext << "callback=" << callback << " ";
	returntext << "filename=" << filename << " ";
	returntext << "text=";
	
	if (textfile.is_open()){
		while(getline(textfile, line)){
			returntext << line << "\n";
		}
	}
	textfile.close();
	returntextstring = returntext.str();
	mg_printf_websocket_frame(nc, WEBSOCKET_OP_TEXT, returntextstring.c_str(), returntextstring.length());
}

static void getFile(struct mg_connection *nc, std::string messagedata){
	std::string filename, callback, brightness, contrast;
	
	getKeyValFromString(&filename, messagedata, "filename");
	getKeyValFromString(&callback, messagedata, "callback");
	getKeyValFromString(&contrast, messagedata, "contrast");
	getKeyValFromString(&brightness, messagedata, "brightness");

	if(filename != ""){
		if (filename.find(".mrc") != std::string::npos) {
			getMRCFile(nc, filename, callback, strtof(contrast.c_str(), NULL), strtof(brightness.c_str(), NULL));
		} else if(filename.find(".log") != std::string::npos) {
			getTextFile(nc, filename, callback);
		}
	}
}

static void saveStackMRC(struct mg_connection *nc, std::string messagedata){
	std::string stackfile, filename, directory, selection;
	std::ifstream infile; 
	std::ofstream outfile;
	struct stat buffer;
	char *readbuffer;
	int *readbufferint;
	int nx, ny, nz;
	std::stringstream error, message;
	
	getKeyValFromString(&stackfile, messagedata, "stackfile");
	getKeyValFromString(&filename, messagedata, "filename");
	getKeyValFromString(&directory, messagedata, "directory");
	getKeyValFromString(&selection, messagedata, "selection");
	
	if(stat(directory.c_str(), &buffer) == 0 ){
		if(filename != ""){
			directory.append("/");
			directory.append(filename);
			//outfile.open(directory.c_str(), std::ifstream::binary|std::ifstream::out|std::ifstream::in);
			outfile.open(directory.c_str(), std::ifstream::binary);
			infile.open(stackfile.c_str(), std::ifstream::binary);
			
			if (infile.is_open() && outfile.is_open()){
				readbuffer = new char [1024];
				infile.read(readbuffer, 1024);
				readbufferint = (int*)readbuffer;
				nx = readbufferint[0];
				ny = readbufferint[1];
				
				outfile.write(readbuffer, 1024);
				delete[] readbuffer;
				
				readbuffer = new char [nx * ny * 4];
				infile.seekg(1024, std::fstream::beg);
				infile.read(readbuffer, nx * ny * 4);
				outfile.write(readbuffer, nx * ny * 4);
				delete [] readbuffer;
				
				readbuffer = new char [4];
				outfile.seekp(8, std::fstream::beg);
				readbufferint = (int*)readbuffer;
				readbufferint[0] = 1;
				outfile.write(readbuffer, 4);
				delete [] readbuffer;
				
				infile.close();
				outfile.close();
			}else{
				error << "Failed to save selection. Couldn't open files";
				returnWSMessage(nc, error.str());
				return;
			}
		}else{
			error << "Failed to save selection. Filename not given";
			returnWSMessage(nc, error.str());
			return;
		}
	} else {
		error << "Failed to save selection. Directory does not exist";
		returnWSMessage(nc, error.str());
		return;
	}
	
	message << "Succesfully saved selection to " << directory;
	returnWSMessage(nc, message.str());
}

static void SQLQuery(struct mg_connection *nc, std::string messagedata){
	sqlite3 *db;
	sqlite3_stmt *stmt;
	int rc;
	std::string dbpath, datastring, user, query, returnstring, callback;
	std::stringstream error, returnstringss, tmpss;

	getUser(&user, messagedata);
	
	getUserHome(&dbpath, user);
	dbpath += "/.simple.sqlite";
	
	getKeyValFromString(&query, messagedata, "sql");
	getKeyValFromString(&callback, messagedata, "callback");

	rc = sqlite3_open(dbpath.c_str(), &db);
	if(rc){
		error << "Failed to open database " << sqlite3_errmsg(db);
		returnWSMessage(nc, error.str());
		return;
	}
	
	rc = sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, NULL);
	if( rc != SQLITE_OK ){
		error << "Failed to query database " << sqlite3_errmsg(db);
		returnWSMessage(nc, error.str());
		return;
	}
	
	returnstringss << "callback=" << callback << " data=";
	
	while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
		returnstringss << "{";
		int id = sqlite3_column_int(stmt, 0);
		returnstringss << '"' << "id" << '"' << ":" << '"' << id << '"';
		for(int i = 1; i < sqlite3_column_count(stmt); i++){
			const char *columnname = sqlite3_column_name(stmt, i);
			const unsigned char *columnvalue = sqlite3_column_text(stmt, i);
			returnstringss << "," << '"' << columnname << '"' << ":" << '"' << columnvalue << '"';
		}
		returnstringss << "};";
	}
	
	if (rc != SQLITE_DONE) {
		error << "Failed to parse query " << sqlite3_errmsg(db);
		returnWSMessage(nc, error.str());
		return;
	}

	sqlite3_finalize(stmt);
	sqlite3_close(db);
	
	returnstring = returnstringss.str();
	mg_printf_websocket_frame(nc, WEBSOCKET_OP_TEXT, returnstring.c_str(), returnstring.length());
	return;
}

static void createJob(struct mg_connection *nc, std::string messagedata, int *outjobid, std::string *jobdir){
	std::string user, projectfolder, projecttable, dbpath, query, jobname, returnstring, jobtype, jobdescription, jobdirectory;
	std::stringstream error, jobidss;
	sqlite3 *db;
	sqlite3_stmt *stmt;
	sqlite3_int64 jobid;
	int rc;
	struct stat buffer;
	
	getUser(&user, messagedata);
	getUserHome(&dbpath, user);
	dbpath += "/.simple.sqlite";
	
	getKeyValFromString(&projectfolder, messagedata, "projectfolder");
	getKeyValFromString(&projecttable, messagedata, "projecttable");
	getKeyValFromString(&jobname, messagedata, "jobname");
	getKeyValFromString(&jobtype, messagedata, "jobtype");
	getKeyValFromString(&jobdescription, messagedata, "jobdescription");
	getKeyValFromString(&jobdirectory, messagedata, "jobdirectory");
	
	rc = sqlite3_open(dbpath.c_str(), &db);
	if(rc){
		error << "Failed to open database " << sqlite3_errmsg(db);
		returnWSMessage(nc, error.str());
		return;
	}
	
	query = "INSERT INTO " + projecttable + "(JOBNAME,JOBSTATUS,JOBTYPE,JOBDESCRIPTION) VALUES ('" + jobname + "', 'Running','" + jobtype + "','" + jobdescription + "');";

	rc = sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, NULL);
	if( rc != SQLITE_OK ){
		error << "Failed to query database " << sqlite3_errmsg(db);
		returnWSMessage(nc, error.str());
		return;
	}
	
	while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
	}
	
	if (rc != SQLITE_DONE) {
		error << "Failed to parse query " << sqlite3_errmsg(db);
		returnWSMessage(nc, error.str());
		return;
	}
	
	jobid = sqlite3_last_insert_rowid(db);
	
	*outjobid = (long long)jobid;

	jobidss << *outjobid;
	*jobdir = projectfolder + "/" + jobidss.str() + "_" + jobtype;
	
	if(stat(jobdir->c_str(), &buffer) != 0 ){
		if(jobdirectory == ""){
			mkdir(jobdir->c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		}else {
			if (symlink(jobdirectory.c_str(), jobdir->c_str()) != 0){
				error << "Failed to link external folder" << jobdirectory.c_str() << " " << jobdir->c_str();
				returnWSMessage(nc, error.str());
				return;
			}
		}
	} else {
		error << "Failed to create job directory";
		returnWSMessage(nc, error.str());
		return;
	}
	
	query = "UPDATE " + projecttable + " SET JOBFOLDER='" + *jobdir +"' WHERE id=" + jobidss.str() + ";";

	rc = sqlite3_prepare_v2(db, query.c_str(), -1, &stmt, NULL);
	if( rc != SQLITE_OK ){
		error << "Failed to query database " << sqlite3_errmsg(db);
		returnWSMessage(nc, error.str());
		return;
	}
	
	while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
	}
	
	if (rc != SQLITE_DONE) {
		error << "Failed to parse query " << sqlite3_errmsg(db);
		returnWSMessage(nc, error.str());
		return;
	}
	
	sqlite3_finalize(stmt);
	
	sqlite3_close(db);
}

static void deleteJob(struct mg_connection *nc, std::string messagedata){
	std::stringstream message;
	std::string removedir;
	SQLQuery(nc, messagedata);
	
	getKeyValFromString(&removedir, messagedata, "removedir");
	std::cout << removedir <<std::endl; // implememt folder deletion
	message << "Job deleted";
	returnWSMessage(nc, message.str());
	
}

static void killJob(struct mg_connection *nc, std::string messagedata){
	std::stringstream message;
	std::string pid;
	SQLQuery(nc, messagedata);
	
	getKeyValFromString(&pid, messagedata, "pid");
	kill(atoi(pid.c_str()), 9);
	
	message << "Job killed";
	returnWSMessage(nc, message.str());
	
}

static void pipeLineView(struct mg_connection *nc, std::string messagedata){
	std::string unidoc, callback, line, tmp, returnstring, selection;
	std::stringstream error, returnstringss;
	std::ifstream unidocstream, selectionstream;
	struct stat buffer;
	
	getKeyValFromString(&unidoc, messagedata, "dir");
	getKeyValFromString(&callback, messagedata, "callback");
	
	if(unidoc == ""){
		error << "Directory not set";
		returnWSMessage(nc, error.str());
		return;
	}
	
	if(callback == ""){
		error << "Callback not set";
		returnWSMessage(nc, error.str());
		return;
	}
	
	unidoc.append("/simple_unidoc.txt");
	
	returnstringss << "callback=" << callback << " data=";
	
	if(stat(unidoc.c_str(), &buffer) == 0){
		unidocstream.open(unidoc.c_str());
		while(getline(unidocstream, line)){
			returnstringss << "{";
			getKeyValFromString(&tmp, line, "intg");
			returnstringss << '"' << "intg" << '"' << ":" << '"' << tmp << '"' << ',';
			getKeyValFromString(&tmp, line, "pspec");
			returnstringss << '"' << "pspec" << '"' << ":" << '"' << tmp << '"' << ',';
			getKeyValFromString(&tmp, line, "thumb");
			returnstringss << '"' << "thumb" << '"' << ":" << '"' << tmp << '"' << ',';
			getKeyValFromString(&tmp, line, "dfx");
			returnstringss << '"' << "dfx" << '"' << ":" << '"' << tmp << '"' << ',';
			getKeyValFromString(&tmp, line, "dfy");
			returnstringss << '"' << "dfy" << '"' << ":" << '"' << tmp << '"' << ',';
			getKeyValFromString(&tmp, line, "angast");
			returnstringss << '"' << "angast" << '"' << ":" << '"' << tmp << '"';
			returnstringss << "};";	
		}
		unidocstream.close();
		
	} else {
		error << "Unidoc doesn't exist";
		returnWSMessage(nc, error.str());
		return;
	}
	
	getKeyValFromString(&selection, messagedata, "dir");
	selection.append("/.selection");
	
	if(stat(selection.c_str(), &buffer) == 0){
		returnstringss << " selected=";
		selectionstream.open(selection.c_str());
		while(getline(selectionstream, line)){
			returnstringss << line << ";";
		}
		selectionstream.close();
	}
	
	returnstring = returnstringss.str();
	mg_printf_websocket_frame(nc, WEBSOCKET_OP_TEXT, returnstring.c_str(), returnstring.length());
}

static void pipeLineSave(struct mg_connection *nc, std::string messagedata){
	std::string filename, directory, data, selection;
	struct stat buffer;
	std::stringstream error, message;
	std::ofstream selectionfile;
	size_t pos = 0;
	
	getKeyValFromString(&filename, messagedata, "filename");
	getKeyValFromString(&directory, messagedata, "directory");
	getKeyValFromString(&data, messagedata, "data");
	
	if(stat(directory.c_str(), &buffer) == 0 ){
		if(filename != ""){
			directory.append("/");
			directory.append(filename);
			selectionfile.open(directory.c_str());
			while ((pos = data.find(";")) != std::string::npos) {
				selection = data.substr(0, pos);
				if (selection != ""){
					if (selectionfile.is_open()){
						selectionfile << selection << "\n";
					}
				}
				data.erase(0, pos + 1);
			}
			selectionfile.close();
		}else{
			error << "Failed to save selection. Filename not given";
			returnWSMessage(nc, error.str());
			return;
		}
	} else {
		error << "Failed to save selection. Directory does not exist";
		returnWSMessage(nc, error.str());
		return;
	}
	
	message << "Succesfully saved selection to " << directory;
	returnWSMessage(nc, message.str());
}

static void import(struct mg_connection *nc, std::string messagedata){
	int jobid;
	std::string jobdir, callback, returnstring, directory, suffix, filetype, projecttable;
	std::stringstream returnmessage, sqlss;
	DIR *directorycontents;
	struct dirent *entry;
	StarFile *starfile = new StarFile();
	
	createJob(nc, messagedata, &jobid, &jobdir);
	
	getKeyValFromString(&projecttable, messagedata, "projecttable");
	getKeyValFromString(&directory, messagedata, "dir");
	getKeyValFromString(&suffix, messagedata, "suffix");
	getKeyValFromString(&filetype, messagedata, "filetype");
	
	returnmessage << "Job started";
	returnWSMessage(nc, returnmessage.str());
	
	if(filetype == "movies"){
		jobdir += "/movies.star";
		starfile->columns.push_back("splmovie");
	}else if(filetype == "micrographs"){
		jobdir += "/micrographs.star";
		starfile->columns.push_back("splmicrograph");
	}else{
		returnmessage << "Error: No file type specified";
		returnWSMessage(nc, returnmessage.str());
		return;
	}
	
	directorycontents = opendir(directory.c_str());
	entry = readdir(directorycontents);
	
	while (entry != NULL){
		if(strstr(entry->d_name, suffix.c_str())){
			starfile->data.push_back(directory + "/" + entry->d_name);
		}
		entry = readdir(directorycontents);
	}
	closedir(directorycontents);
	
	writeStarFile(starfile, jobdir);
	
	sqlss << "UPDATE^" << projecttable << "^SET^JOBSTATUS='Complete'^WHERE^id=" << jobid << ";";
	SQLQuery(nc, sqlss.str());
}

static void externalJob(struct mg_connection *nc, std::string messagedata){
	int jobid;
	std::string jobdir, callback, returnstring, jobdirectory, sqlstr, projecttable;
	std::stringstream returnmessage, returnstringss, error, sql;
	
	createJob(nc, messagedata, &jobid, &jobdir);
	
	getKeyValFromString(&projecttable, messagedata, "projecttable");
	
	getKeyValFromString(&callback, messagedata, "callback");

	returnmessage << "External Job Added";
	returnWSMessage(nc, returnmessage.str());
	
	returnstringss << "callback=" << callback;
	returnstring = returnstringss.str();
	sql << "sql=UPDATE^" << projecttable << "^SET^JOBSTATUS='External'^WHERE^id=" << jobid << ";";
	sqlstr = sql.str();
	std::cout << sqlstr << std::endl;
	SQLQuery(nc, sqlstr.c_str());
	mg_printf_websocket_frame(nc, WEBSOCKET_OP_TEXT, returnstring.c_str(), returnstring.length());
}

static void fileSelector(struct mg_connection *nc, std::string messagedata){
	std::string directory, user, callback, returnstring, tmpstring;
	std::stringstream error, returnstringss;
	struct stat buffer;
	DIR *directorycontents;
	struct dirent *entry;
	std::vector<std::string> directories, files;
	std::size_t location;
	char character;
	
	getKeyValFromString(&directory, messagedata, "directory");
	getKeyValFromString(&callback, messagedata, "callback");
	
	if(directory == ""){
		getUser(&user, messagedata);
		getUserHome(&directory, user);
	}
	
	if(stat(directory.c_str(), &buffer) != 0){
		error << "Directory doesn't exist";
		returnWSMessage(nc, error.str());
		return;
	}
	
	directorycontents = opendir(directory.c_str());
	entry = readdir(directorycontents);
	while (entry != NULL){
		if (entry->d_type == DT_DIR){
			if(entry->d_name[0] != '.'){
				directories.push_back(entry->d_name);
			} 
		}else{
			if(entry->d_name[0] != '.'){
				files.push_back(entry->d_name);
			}
		}
		entry = readdir(directorycontents);
	}
	closedir(directorycontents);
	
	returnstringss << "callback=" << callback << " directory=" << directory;
	
	if(callback == "fileSelectorFileCallback"){
		returnstringss << " files=";
		std::sort (files.begin(), files.end());
		
		for (std::vector<std::string>::iterator it=files.begin(); it!=files.end(); ++it){
			returnstringss << directory << "/";
			tmpstring = *it;
			location = 0;
			while ( location < tmpstring.length() ){
				character = tmpstring.c_str()[location];
				if (character == ' '){
					returnstringss << "^"; 
				} else if (character != ' '){
					returnstringss << character;
				} else {
					break;
				}
				location = location + 1;
			}
			returnstringss << ";";
		}
	}
	
	std::sort (directories.begin(), directories.end());
	returnstringss << " directories=";
	for (std::vector<std::string>::iterator it=directories.begin(); it!=directories.end(); ++it){
		returnstringss << directory << "/";
		tmpstring = *it;
		location = 0;
		while ( location < tmpstring.length() ){
			character = tmpstring.c_str()[location];
			if (character == ' '){
				returnstringss << "^"; 
			} else if (character != ' '){
				returnstringss << character;
			} else {
				break;
			}
			location = location + 1;
		}
			returnstringss << ";";
	}
	
	returnstring = returnstringss.str();
	mg_printf_websocket_frame(nc, WEBSOCKET_OP_TEXT, returnstring.c_str(), returnstring.length());
	
}

static void createSimpleEnv(std::string jobdir, std::string messagedata){
	std::ofstream envfile;
	std::string tmp;
	
	jobdir.append("/simple_distr_config.env");
	envfile.open(jobdir.c_str());
	envfile << "simple_path = " << SIMPLE_DIR_VALUE << "\n";
	
	getKeyValFromString(&tmp, messagedata, "time_per_image");
	if(tmp != ""){
		envfile << "time_per_image = " << tmp << "\n";
	}
	
	getKeyValFromString(&tmp, messagedata, "qsys_qos");
	if(tmp != ""){
		envfile << "qsys_qos = " << tmp << "\n";
	}
	
	getKeyValFromString(&tmp, messagedata, "qsys_name");
	if(tmp != ""){
		envfile << "qsys_name = " << tmp << "\n";
	}
	
	getKeyValFromString(&tmp, messagedata, "qsys_partition");
	if(tmp != ""){
		envfile << "qsys_partition = " << tmp << "\n";
	}
	
	getKeyValFromString(&tmp, messagedata, "job_ntasks");
	if(tmp != ""){
		envfile << "job_ntasks = " << tmp << "\n";
	}
	
	getKeyValFromString(&tmp, messagedata, "job_memory_per_task");
	if(tmp != ""){
		envfile << "job_memory_per_task = " << tmp << "\n";
	}
	
	getKeyValFromString(&tmp, messagedata, "job_ntasks_per_socket");
	if(tmp != ""){
		envfile << "job_ntasks_per_socket = " << tmp << "\n";
	}
	
	getKeyValFromString(&tmp, messagedata, "job_cpus_per_task");
	if(tmp != ""){
		envfile << "job_cpus_per_task = " << tmp << "\n";
	}
	
	envfile.close();
}

static void simpleJob(struct mg_connection *nc, std::string messagedata){
	FILE* pipe;
	int pid, pid2, jobid, returnstatus;
	char buffer[512];
	std::stringstream ss, sspid, sqlss, error;
	std::string argscommand, line, command, arg, jobdir, path, prg, projecttable, moviesdir, starfilename, unidocfilename;
	std::vector<std::string> commandargs, unidocs;
	std::vector<std::string>::iterator it, it2;
	int columnid;
	std::string element;
	std::ofstream outfile;
	struct stat buffer2;
	DIR *dir;
	struct dirent *entry;
	bool unidocprocessed;
	
	createJob(nc, messagedata, &jobid, &jobdir);
	getKeyValFromString(&projecttable, messagedata, "projecttable");
	
	pid = fork();
	
	if(pid == 0){
		path = std::getenv("PATH");
		path.append(":");
		path.append(SIMPLE_DIR_VALUE);
		path.append("/bin");
		setenv("PATH", path.c_str(), 1);
		 
		getKeyValFromString(&argscommand, messagedata, "exec");
		argscommand.append(" prg="); 
		getKeyValFromString(&prg, messagedata, "prg");
		argscommand.append(prg);
	
		pipe = popen(argscommand.c_str(), "r");
		while (!feof(pipe)) {
			if (fgets(buffer, sizeof(buffer), pipe) != NULL){
				ss << buffer;
			}
		}
	
		pclose(pipe);
	
		while(std::getline(ss,line,'\n')){
			getElementFromString(&arg, line, 0);
			commandargs.push_back(arg);
		}
	
		createSimpleEnv(jobdir, messagedata);
		
		for (it=commandargs.begin(); it!=commandargs.end(); ++it){
			getKeyValFromString(&arg, messagedata, *it);
			if(arg != "" && *it != ""){
				argscommand.append(" ");
				argscommand.append(*it);
				argscommand.append("=");
				argscommand.append(arg);
			}
		}
		
		chdir(jobdir.c_str());
		
		// Per program setup
		if(prg == "unblur"){
			getKeyValFromString(&starfilename, messagedata, "starfile");
			StarFile* starfile = new StarFile();
			readStarFile(starfile, starfilename);
			outfile.open("movies.txt");
			if(outfile.is_open()){
				getStarColumn(starfile, "splMovie", &columnid);
				for(it = starfile->data.begin() ; it != starfile->data.end(); ++it){
					getElementFromString(&element, *it, columnid);
					outfile << element << "\n";
				}
				outfile.close();
				argscommand.append(" filetab=movies.txt");
			} else {
				return;
			}
			delete starfile;
		}else if(prg == "unblur_ctffind"){
			getKeyValFromString(&starfilename, messagedata, "starfile");
			StarFile* starfile = new StarFile();
			readStarFile(starfile, starfilename);
			outfile.open("movies.txt");
			if(outfile.is_open()){
				getStarColumn(starfile, "splMovie", &columnid);
				for(it = starfile->data.begin() ; it != starfile->data.end(); ++it){
					getElementFromString(&element, *it, columnid);
					outfile << element << "\n";
				}
				outfile.close();
				argscommand.append(" filetab=movies.txt");
			} else {
				return;
			}
			delete starfile;
		}else if(prg == "ctffind"){
			getKeyValFromString(&starfilename, messagedata, "starfile");
			StarFile* starfile = new StarFile();
			readStarFile(starfile, starfilename);
			outfile.open("micrographsforctf.txt");
			if(outfile.is_open()){
				getStarColumn(starfile, "splForCTFName", &columnid); // NEEDS TO CHECK FOR EXISTENCE, ELSE USE MICROGRASPHNAME
				for(it = starfile->data.begin() ; it != starfile->data.end(); ++it){
					getElementFromString(&element, *it, columnid);
					outfile << element << "\n";
				}
				outfile.close();
				argscommand.append(" filetab=micrographsforctf.txt");
			} else {
				return;
			}
			delete starfile;
		}else if(prg == "preproc"){
			getKeyValFromString(&moviesdir, messagedata, "moviesdir");
			if (symlink(moviesdir.c_str(), "movies") != 0){ // SOMETHING AMISS HERE
				error << "Failed to link movies directory " << moviesdir;
				returnWSMessage(nc, error.str());
				//return; 
			}
			argscommand.append(" dir_movies=movies dir_target=pipeline");
			pid2 = fork();
			if(pid2 == 0){
				StarFile* starfile = new StarFile();
				
				starfile->columns.push_back("splMovie");
				starfile->columns.push_back("splMicrographName");
				starfile->columns.push_back("splForCTFName");
				starfile->columns.push_back("splThumbName");
				starfile->columns.push_back("splPspecName");
				starfile->columns.push_back("splSamplingDistance");
				starfile->columns.push_back("rlnMagnification");
				starfile->columns.push_back("rlnDetectorPixelSize");
				starfile->columns.push_back("splVoltage");
				starfile->columns.push_back("splDefocusX");
				starfile->columns.push_back("splDefocusY");
				starfile->columns.push_back("rlnDefocusU");
				starfile->columns.push_back("rlnDefocusV");
				starfile->columns.push_back("splSphericalAberration");
				starfile->columns.push_back("splDefocusAngle");
				starfile->columns.push_back("splCtfMaxResolution");
				starfile->columns.push_back("splAmplitudeContrast");
				
				while(true){
					//LIST UNIDOC PARTS THEN TEST TO SEE IF PROCESSED. IF NOT ADD TO STARFILE
					dir = opendir("pipeline");
					entry = readdir(dir);
					while (entry != NULL){
						unidocfilename = entry->d_name;
						if(unidocfilename.find("unidoc_output_part") != std::string::npos){
							unidocprocessed = true;
							for(it = unidocs.begin() ; it != unidocs.end(); ++it){
								if(*it == unidocfilename){
									unidocprocessed = false;
									break;
								}
							}
							if(unidocprocessed){
								unidocs.push_back(unidocfilename);
								UniDoc* unidoc = new UniDoc(); 
								readUniDoc(unidoc, unidocfilename);
								// ADD INFO FROM UNIDOC TO STARFUILE
								starfile->data.push_back(unidoc->data[0]);
								delete unidoc;
								writeStarFile(starfile, "micrographs_preproc.star");
							}
						}
					}
					usleep(60000000);
				}
			}
		}
		
		argscommand.append(" 2>&1 > job.log");
		
		std::cout << argscommand << std::endl;
		pipe = popen(argscommand.c_str(), "r");
		returnstatus = pclose(pipe);
		
		sqlss << "callback=newJobCallback sql=";
		
		// Per Job cleanup and post processing
		if(prg == "unblur"){
			if(stat("simple_unidoc.txt", &buffer2) == 0){
				getKeyValFromString(&starfilename, messagedata, "starfile");
				StarFile* starfile = new StarFile();
				UniDoc* unidoc = new UniDoc();
				std::stringstream ss;
				
				readStarFile(starfile, starfilename);
				readUniDoc(unidoc, "simple_unidoc.txt");
				
				starfile->columns.push_back("splMicrographName");
				starfile->columns.push_back("splForCTFName");
				starfile->columns.push_back("splThumbName");
				starfile->columns.push_back("splPspecName");
				starfile->columns.push_back("splSamplingDistance");
				starfile->columns.push_back("rlnMagnification");
				starfile->columns.push_back("rlnDetectorPixelSize");
				
				int moviecolumn;
				getStarColumn(starfile, "splMovie", &moviecolumn);
				std::string moviename, unidocline, intg, forctf, thumb, pspec, smpd;
				
				for(it = starfile->data.begin() ; it != starfile->data.end(); ++it){
					getElementFromString(&moviename, *it, moviecolumn);
					for(it2 = unidoc->data.begin() ; it2 != unidoc->data.end(); ++it2){
						unidocline = *it2;
						if(unidocline.find(moviename) != std::string::npos){
							getKeyValFromString(&intg, unidocline, "intg");
							getKeyValFromString(&forctf, unidocline, "forctf");
							getKeyValFromString(&thumb, unidocline, "thumb");
							getKeyValFromString(&pspec, unidocline, "pspec");
							getKeyValFromString(&smpd, unidocline, "smpd");
							*it = *it + " " + jobdir + "/" + intg + " " + jobdir + "/" + forctf + " " + jobdir + "/" + thumb + " " + jobdir + "/" + pspec + " " + smpd + " 10000 " + smpd;
							break;
						}
					}
				}
				writeStarFile(starfile, "micrographs.star");
				delete starfile;
				delete unidoc;
				sqlss << "UPDATE^" << projecttable << "^SET^JOBSTATUS='Complete'^WHERE^id=" << jobid << ";";
				SQLQuery(nc, sqlss.str());
			} else {
				sqlss << "UPDATE^" << projecttable << "^SET^JOBSTATUS='Failed'^WHERE^id=" << jobid << ";";
				SQLQuery(nc, sqlss.str());
			}
		}else if(prg == "unblur_ctffind"){
			if(stat("simple_unidoc.txt", &buffer2) == 0){
				getKeyValFromString(&starfilename, messagedata, "starfile");
				StarFile* starfile = new StarFile();
				UniDoc* unidoc = new UniDoc();
				
				readStarFile(starfile, starfilename);
				readUniDoc(unidoc, "simple_unidoc.txt");
				
				starfile->columns.push_back("splMicrographName");
				starfile->columns.push_back("splForCTFName");
				starfile->columns.push_back("splThumbName");
				starfile->columns.push_back("splPspecName");
				starfile->columns.push_back("splSamplingDistance");
				starfile->columns.push_back("rlnMagnification");
				starfile->columns.push_back("rlnDetectorPixelSize");
				starfile->columns.push_back("splVoltage");
				starfile->columns.push_back("splDefocusX");
				starfile->columns.push_back("splDefocusY");
				starfile->columns.push_back("rlnDefocusU");
				starfile->columns.push_back("rlnDefocusV");
				starfile->columns.push_back("splSphericalAberration");
				starfile->columns.push_back("splDefocusAngle");
				starfile->columns.push_back("splCtfMaxResolution");
				starfile->columns.push_back("splAmplitudeContrast");

				int moviecolumn;
				getStarColumn(starfile, "splMovie", &moviecolumn);
				std::string moviename, unidocline, intg, forctf, thumb, pspec, smpd, voltage, dfx, dfy, angast, cs, ctfres, fraca;
				
				for(it = starfile->data.begin() ; it != starfile->data.end(); ++it){
					getElementFromString(&moviename, *it, moviecolumn);
					for(it2 = unidoc->data.begin() ; it2 != unidoc->data.end(); ++it2){
						unidocline = *it2;
						if(unidocline.find(moviename) != std::string::npos){
							getKeyValFromString(&intg, unidocline, "intg");
							getKeyValFromString(&forctf, unidocline, "forctf");
							getKeyValFromString(&thumb, unidocline, "thumb");
							getKeyValFromString(&pspec, unidocline, "pspec");
							getKeyValFromString(&smpd, unidocline, "smpd");
							getKeyValFromString(&voltage, unidocline, "voltage");
							getKeyValFromString(&dfx, unidocline, "dfx");
							getKeyValFromString(&dfy, unidocline, "dfy");
							getKeyValFromString(&angast, unidocline, "angast");
							getKeyValFromString(&cs, unidocline, "cs");
							getKeyValFromString(&ctfres, unidocline, "ctfres");
							getKeyValFromString(&fraca, unidocline, "fraca");
							*it = *it + " " + jobdir + "/" + intg + " " + jobdir + "/" + forctf + " " + jobdir + "/" + thumb + " " + jobdir + "/" + pspec + " " + smpd + " 10000 " + smpd + " " + voltage + " " + dfx + " " + dfx + " " + dfy + " " + dfy + " " + cs + " " + angast + " " + ctfres + " " + fraca;
							break;
						}
					}
				}
				writeStarFile(starfile, "micrographs_ctf.star");
				delete starfile;
				delete unidoc;
				sqlss << "UPDATE^" << projecttable << "^SET^JOBSTATUS='Complete'^WHERE^id=" << jobid << ";";
				SQLQuery(nc, sqlss.str());
			} else {
				sqlss << "UPDATE^" << projecttable << "^SET^JOBSTATUS='Failed'^WHERE^id=" << jobid << ";";
				SQLQuery(nc, sqlss.str());
			}
		}else if(prg == "ctffind"){
			if(stat("simple_unidoc.txt", &buffer2) == 0){
				/*
				getKeyValFromString(&starfilename, messagedata, "starfile");
				StarFile* starfile = new StarFile();
				UniDoc* unidoc = new UniDoc();
				
				readStarFile(starfile, starfilename);
				readUniDoc(unidoc, "simple_unidoc.txt");
				
				starfile->columns.push_back("splMicrographName");
				starfile->columns.push_back("splForCTFName");
				starfile->columns.push_back("splThumbName");
				starfile->columns.push_back("splPspecName");
				starfile->columns.push_back("splSamplingDistance");
				starfile->columns.push_back("rlnMagnification");
				starfile->columns.push_back("rlnDetectorPixelSize");
				starfile->columns.push_back("splVoltage");
				starfile->columns.push_back("splDefocusX");
				starfile->columns.push_back("splDefocusY");
				starfile->columns.push_back("rlnDefocusU");
				starfile->columns.push_back("rlnDefocusV");
				starfile->columns.push_back("splSphericalAberration");
				starfile->columns.push_back("splDefocusAngle");
				starfile->columns.push_back("splCtfMaxResolution");
				starfile->columns.push_back("splAmplitudeContrast");

				int moviecolumn;
				getStarColumn(starfile, "splMovie", &moviecolumn);
				std::string moviename, unidocline, intg, forctf, thumb, pspec, smpd, voltage, dfx, dfy, angast, cs, ctfres, fraca;
				
				for(it = starfile->data.begin() ; it != starfile->data.end(); ++it){
					getElementFromString(&moviename, *it, moviecolumn);
					for(it2 = unidoc->data.begin() ; it2 != unidoc->data.end(); ++it2){
						unidocline = *it2;
						if(unidocline.find(moviename) != std::string::npos){
							getKeyValFromString(&intg, unidocline, "intg");
							getKeyValFromString(&forctf, unidocline, "forctf");
							getKeyValFromString(&thumb, unidocline, "thumb");
							getKeyValFromString(&pspec, unidocline, "pspec");
							getKeyValFromString(&smpd, unidocline, "smpd");
							getKeyValFromString(&voltage, unidocline, "voltage");
							getKeyValFromString(&dfx, unidocline, "dfx");
							getKeyValFromString(&dfy, unidocline, "dfy");
							getKeyValFromString(&angast, unidocline, "angast");
							getKeyValFromString(&cs, unidocline, "cs");
							getKeyValFromString(&ctfres, unidocline, "ctfres");
							getKeyValFromString(&fraca, unidocline, "fraca");
							*it = *it + " " + jobdir + "/" + intg + " " + jobdir + "/" + forctf + " " + jobdir + "/" + thumb + " " + jobdir + "/" + pspec + " " + smpd + " 10000 " + smpd + " " + voltage + " " + dfx + " " + dfx + " " + dfy + " " + dfy + " " + cs + " " + angast + " " + ctfres + " " + fraca;
							break;
						}
					}
				}
				writeStarFile(starfile, "micrographs_ctf.star");
				delete starfile;
				delete unidoc;
				sqlss << "UPDATE^" << projecttable << "^SET^JOBSTATUS='Complete'^WHERE^id=" << jobid << ";";
				SQLQuery(nc, sqlss.str());
				*/ 
			} else {
				sqlss << "UPDATE^" << projecttable << "^SET^JOBSTATUS='Failed'^WHERE^id=" << jobid << ";";
				SQLQuery(nc, sqlss.str());
			}
		}
		
		/*if(WEXITSTATUS(returnstatus) == 0){
			sqlss << "UPDATE^" << projecttable << "^SET^JOBSTATUS='Complete'^WHERE^id=" << jobid << ";";
			
			
			
			
		}else{
			sqlss << "UPDATE^" << projecttable << "^SET^JOBSTATUS='Failed'^WHERE^id=" << jobid << ";";
		}*/

		_exit(0);
		
	}else if(pid > 0){
		signal(SIGCHLD, SIG_IGN);
		
		sqlss << "callback=refreshHistory sql=UPDATE^" << projecttable << "^SET^JOBPID=" << pid << "^WHERE^id=" << jobid << ";";
		std::cout << sqlss.str() << std::endl;
		SQLQuery(nc, sqlss.str());
	}
}

static void syncJob(struct mg_connection *nc, std::string messagedata){
	std::string sourcedir, destinationdir, destinationdirgluster, jobdir, projecttable, command, line, stoptime, killtime, ncunits, email, logfile;
	std::stringstream sqlss, ss;
	int jobid, pid, returnstatus, srccount, srccountlast, lastsync, killtimeepoch;
	FILE* pipe;
	char buffer[512];
	size_t pos;
	std::ofstream logoutput;
	
	createJob(nc, messagedata, &jobid, &jobdir);
	getKeyValFromString(&sourcedir, messagedata, "sourcedir");
	getKeyValFromString(&destinationdir, messagedata, "destinationdir");
	getKeyValFromString(&projecttable, messagedata, "projecttable");
	getKeyValFromString(&stoptime, messagedata, "stoptime");
	getKeyValFromString(&killtime, messagedata, "killtime");
	getKeyValFromString(&ncunits, messagedata, "ncunits");
	getKeyValFromString(&email, messagedata, "email");
	
	pid = fork();
	
	logfile = jobdir + "/job.log";
	
	destinationdir.append("/raw");
	destinationdirgluster = destinationdir;
	
	pos = destinationdir.find("beegfs");
	if(pos != std::string::npos){
		destinationdirgluster.replace(pos, 6, "glusterfs");
		command = "mkdir -p " + destinationdirgluster;
		system(command.c_str());
		symlink(destinationdirgluster.c_str(), destinationdir.c_str());
		destinationdir = destinationdirgluster;
	} else {
		command = "mkdir -p " + destinationdir;
		system(command.c_str());
	}
	
	srccountlast = 0;
	lastsync = std::time(0);
	
	if(pid == 0){
		killtimeepoch = std::time(0) + std::atoi(killtime.c_str());
		while(std::time(0) < killtimeepoch){
			command = "rsync -rltpv " + sourcedir + " " + destinationdir + " | wc -l";
			logoutput.open(logfile.c_str(), std::ofstream::out | std::ofstream::app);
			if(logoutput.is_open()){
				logoutput << std::time(0) << " " << "Starting new sync cycle ..." << "/n";
				logoutput <<"\t Sync command : " << command << "/n";
				logoutput <<"\t Seconds until kill time : " << killtimeepoch - std::time(0) << "/n";
				logoutput.close();
			}
			
			pipe = popen(command.c_str(), "r");	
			while (!feof(pipe)) {
				if (fgets(buffer, sizeof(buffer), pipe) != NULL){
					srccount = std::atoi(buffer);
				}
			}
			logoutput.open(logfile.c_str(), std::ofstream::out | std::ofstream::app);
			if(logoutput.is_open()){
				logoutput <<"\t Files synced in cycle : " << srccount - 4 << "/n";
				logoutput.close();
			}
			returnstatus = pclose(pipe);
			usleep(1000000);
			
			if(WEXITSTATUS(returnstatus) != 0){
				sqlss << "UPDATE^" << projecttable << "^SET^JOBSTATUS='Failed'^WHERE^id=" << jobid << ";";
				SQLQuery(nc, sqlss.str());
				command = "mail -s 'CRESCENT: Your sync job has encountered an error' " + email +" <<< 'An error was encountered'";
				system(command.c_str());
				logoutput.open(logfile.c_str(), std::ofstream::out | std::ofstream::app);
				if(logoutput.is_open()){
					logoutput <<"Critical error: Sync has shutdown" << "/n";
					logoutput.close();
				}
				_exit(0);
			}
			
			if (srccount > srccountlast){
				srccountlast = srccount;
				lastsync = std::time(0);
				logoutput.open(logfile.c_str(), std::ofstream::out | std::ofstream::app);
				if(logoutput.is_open()){
					logoutput <<"\t Sleeping for 30 seconds before next sync cycle" << "/n";
					logoutput.close();
				}
				usleep(30000000);
			} else {
				if(std::time(0) > (std::atoi(stoptime.c_str()) + lastsync)){
					break;
				}
				logoutput.open(logfile.c_str(), std::ofstream::out | std::ofstream::app);
				if(logoutput.is_open()){
					logoutput <<"\t No new files. Sleeping for 120 seconds before next sync cycle" << "/n";
					logoutput.close();
				}
				usleep(120000000);
			}	
		}
		logoutput.open(logfile.c_str(), std::ofstream::out | std::ofstream::app);
		if(logoutput.is_open()){
			logoutput <<"No new files seen within stop time. Terminating sync" << "/n";
			logoutput.close();
		}
		sqlss << "UPDATE^" << projecttable << "^SET^JOBSTATUS='Complete'^WHERE^id=" << jobid << ";";
		SQLQuery(nc, sqlss.str());
		command = "mail -s 'CRESCENT: Your sync job has completed' " + email + " <<< 'Sync job completed'";
		system(command.c_str());
		_exit(0);
	} else if (pid > 0){
		signal(SIGCHLD, SIG_IGN);
		sqlss << "callback=refreshHistory sql=UPDATE^" << projecttable << "^SET^JOBPID=" << pid << "^WHERE^id=" << jobid << ";";
		SQLQuery(nc, sqlss.str());
	}	
}

static void processWebsocketMessage(struct mg_connection *nc, std::string messagedata){
	std::string command; 
	getKeyValFromString(&command, messagedata, "cmd");
	if(command == "getfile"){
		getFile(nc,messagedata);
	}else if(command == "shutdown"){
		s_signal_received = 1;
	}else if(command == "sql"){
		SQLQuery(nc, messagedata);
	}else if(command == "pipelineview"){
		pipeLineView(nc, messagedata);
	}else if(command == "import"){
		import(nc, messagedata);
	}else if(command == "external"){
		externalJob(nc, messagedata);
	}else if(command == "pipelinesave"){
		pipeLineSave(nc, messagedata);
	}else if(command == "fileselector"){
		fileSelector(nc, messagedata);
	}else if(command == "simplejob"){
		simpleJob(nc, messagedata);
	}else if(command == "syncjob"){
		syncJob(nc, messagedata);
	}else if(command == "deletejob"){
		deleteJob(nc, messagedata);
	}else if(command == "killjob"){
		killJob(nc, messagedata);
	}else if(command == "savestack"){
		saveStackMRC(nc, messagedata);
	}
}

static void eventHandler(struct mg_connection *connection, int event, void *eventdata) {
	switch (event) {
		case MG_EV_HTTP_REQUEST: {
			printf("Serving HTTP data on port %s\n", s_http_port);
			mg_serve_http(connection, (struct http_message *) eventdata, s_http_server_opts);
			break;
		}
		case MG_EV_WEBSOCKET_HANDSHAKE_DONE: {
			printf("Websocket opened on port %s\n", s_http_port);
			break;
		}
		case MG_EV_WEBSOCKET_FRAME: {
			char buf[1000000];
			std::string messagedata;
			struct websocket_message *message = (struct websocket_message *) eventdata;
			struct mg_str messagestring = {(char *) message->data, message->size};
			snprintf(buf, sizeof(buf), "%.*s", (int) messagestring.len, messagestring.p);
			messagedata = buf;
			processWebsocketMessage(connection, messagedata);
			break;
		}
		case MG_EV_CLOSE: {
			printf("Websocket closed on port %s\n", s_http_port);
			break;
		}
	}
}

int main(void) {
	struct mg_mgr manager;	// Define Mongoose manager										
	struct mg_connection *connection;	// Define Mongoose connection
	  
	mg_mgr_init(&manager, NULL);	// Initialise Mongoose manager
	  
	printf("Starting Simple server on port %s\n", s_http_port);	// Local output
	  
	connection = mg_bind(&manager, s_http_port, eventHandler);	// Create Mongoose connection

	if (connection == NULL) {
		printf("Failed to create listener\n");
		return 1;
	}

	mg_set_protocol_http_websocket(connection);	// Set Mongoose protocol
	  
	s_http_server_opts.document_root = "www";
	s_http_server_opts.enable_directory_listing = "no";

	while (s_signal_received == 0) {
		mg_mgr_poll(&manager, 200);
	}
	mg_mgr_free(&manager);
	
	return 0;
}
