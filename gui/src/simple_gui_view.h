#ifndef SIMPLE_GUI_VIEW_H_INCLUDED
#define SIMPLE_GUI_VIEW_H_INCLUDED
/*
#include "simple_gui_mrcdata.h"
*/

#include "simple_gui_stardata.h"
#include "simple_gui_textdata.h"
#include "simple_gui_mrcdata.h"

#include <string>
#include <cstring>
#include <jpeglib.h>
#include <libgen.h>


extern "C" {
    #include "ext/mongoose/mongoose.h"
}

struct JPEGData {
	unsigned char 	*pixels;
	unsigned char 	*jpeg;
	unsigned long	jpegsize;
};

void encodeMRCJPEG(MRCFile *mrcfile, JPEGData *jpegdata, int size){
	struct jpeg_compress_struct 	cinfo;
	struct jpeg_error_mgr 			jerr;
	int								numerator;
	int 							newsize;
	
	for(numerator = 1; numerator <= 8; numerator++){
		newsize = mrcfile->miniheader[0] * numerator / 8;
		if(newsize > size){
			std::cout << size << " " << newsize << " " << numerator << std::endl;
			break;
		}
	}
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);
	jpeg_mem_dest(&cinfo, &jpegdata->jpeg, &jpegdata->jpegsize);
	cinfo.image_width = mrcfile->miniheader[0];
	cinfo.image_height = mrcfile->miniheader[1];
	cinfo.input_components = 1;
	cinfo.in_color_space = JCS_GRAYSCALE;
	jpeg_set_defaults(&cinfo);
	cinfo.scale_num = numerator;
	cinfo.scale_denom = 8;
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

void processMRCImage(MRCFile *mrcfile, JPEGData *jpegdata, float contrast, float brightness){
	
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

void MRCToJPEG(struct mg_connection* http_connection, struct http_message* message){
	char 	filename[2056];
	char	contrast[8];
	char	brightness[8];
	char	size[8];
	char	frameid[8];
	struct 	stat buffer;   
 
	mg_get_http_var(&message->query_string, "filename", filename, sizeof(filename));
	mg_get_http_var(&message->query_string, "contrast", contrast, sizeof(contrast));
	mg_get_http_var(&message->query_string, "brightness", brightness, sizeof(brightness));
	mg_get_http_var(&message->query_string, "size", size, sizeof(size));
	mg_get_http_var(&message->query_string, "frameid", frameid, sizeof(frameid));
	
	MRCFile* mrcfile = new MRCFile();
	JPEGData* jpegdata = new JPEGData();
	if(stat(filename, &buffer) == 0){
		readMRCFile(mrcfile, filename, atoi(frameid));
		processMRCImage(mrcfile, jpegdata, atof(contrast), atof(brightness));
		encodeMRCJPEG(mrcfile, jpegdata, atoi(size));
	}
	delete [] mrcfile->data;
	delete mrcfile;
	mg_send_head(http_connection, 200, jpegdata->jpegsize, "Content-Type: image/jpeg");
	mg_send(http_connection, jpegdata->jpeg, jpegdata->jpegsize);
	delete [] jpegdata->pixels;
	delete [] jpegdata->jpeg;
	delete jpegdata;
}

std::string starFileViewer(char* filename){
	// NEEDS WORK TO ONLY SEND DATA REQUIRED NOT WHOLE STAR FILE!
	std::string JSONstring;
	std::vector<std::string>::iterator 	columnsit;
	std::vector<std::string>::iterator 	datait;
	bool								prefixcolumnscomma;
	bool								prefixdatacomma;
	std::string							column;
	int									dataid;
	int									frameid;
	StarFile*    						starfile;
 
	starfile = new StarFile();
	readStarFile(starfile, filename);
	
	JSONstring.append("\"rootdirectory\":\"");
	JSONstring.append(dirname(filename));
	JSONstring.append("\",");
	JSONstring.append("\"displayoptions\":[");
	prefixcolumnscomma = false;
	for(columnsit = starfile->columns.begin() ; columnsit != starfile->columns.end(); ++columnsit){
		column = *columnsit;
		if (std::strstr(column.c_str(), "MicrographName") ||
		std::strstr(column.c_str(), "ForCTFName") ||
		std::strstr(column.c_str(), "ThumbName") ||
		std::strstr(column.c_str(), "CtfImage") ||
		std::strstr(column.c_str(), "CTFFIND") ||
		std::strstr(column.c_str(), "PspecName") ||
		std::strstr(column.c_str(), "ReferenceImage")) {
			if (prefixcolumnscomma) {
				JSONstring.append(",");
			} else {
				prefixcolumnscomma = true;
			}
			JSONstring.append("\"" + *columnsit + "\"");
		}
	}
	JSONstring.append("],");
	JSONstring.append("\"sortoptions\":[");
	prefixcolumnscomma = false;
	for(columnsit = starfile->columns.begin() ; columnsit != starfile->columns.end(); ++columnsit){
		if (prefixcolumnscomma) {
			JSONstring.append(",");
		} else {
			prefixcolumnscomma = true;
		}
		JSONstring.append("\"" + *columnsit + "\"");
	}

	JSONstring.append("],");
	
	JSONstring.append("\"snapshots\":[");
	prefixdatacomma = false;
	dataid = 0;
	for(datait = starfile->data.begin() ; datait != starfile->data.end(); ++datait){
		if (prefixdatacomma) {
			JSONstring.append(",");
		} else {
			prefixdatacomma = true;
		}
		JSONstring.append("{");

		for(unsigned int i = 0; i < starfile->columns.size(); i++){
			if (getElementFromString(*datait, i).find('@') != std::string::npos){
				JSONstring.append("\"frameid\":\"");
				frameid = atoi(getElementFromString(*datait, i).substr(0, getElementFromString(*datait, i).find('@')).c_str()) -1 ;
			//	JSONstring.append(getElementFromString(*datait, i).substr(0, getElementFromString(*datait, i).find('@')));
				JSONstring.append(std::to_string(frameid));
				JSONstring.append("\",\"" + starfile->columns[i] + "\":\"");
				JSONstring.append(getElementFromString(*datait, i).substr(getElementFromString(*datait, i).find('@') + 1 , getElementFromString(*datait, i).length()));
				JSONstring.append("\",");
			} else {
				JSONstring.append("\"" + starfile->columns[i] + "\":\"" + getElementFromString(*datait, i) + "\",");
			}
		}
		JSONstring.append("\"id\":\"");
		JSONstring.append(std::to_string(dataid));
		JSONstring.append("\"");
		dataid++;
		JSONstring.append("}");
	}
	
	JSONstring.append("]");
	
	delete starfile;
	
	return JSONstring;
	
}

std::string mrcFileViewer(char* filename){
	std::string 						JSONstring;
	bool								prefixdatacomma;
	int									dataid;
	MRCFile*    						mrcfile;
	char*								dirc;
	char*								basec;
	char*								dir;
	char*								base;
	
	mrcfile = new MRCFile();
	readMRCFileHeader(mrcfile, filename);
	dirc = strdup(filename);
	basec = strdup(filename);
	dir = dirname(dirc);
	base = basename(basec);
	
	JSONstring.append("\"displayoptions\":[\"filename\"],");
	JSONstring.append("\"sortoptions\":[\"frameid\", \"filename\"],");
	JSONstring.append("\"snapshots\":[");
	prefixdatacomma = false;
	dataid = 0;
	for(dataid = 0; dataid < mrcfile->miniheader[2]; dataid++){
		if (prefixdatacomma) {
			JSONstring.append(",");
		} else {
			prefixdatacomma = true;
		}
		JSONstring.append("{");
		JSONstring.append("\"filename\":\"");
		JSONstring.append(base);
		JSONstring.append("\",");
		JSONstring.append("\"frameid\":\"");
		JSONstring.append(std::to_string(dataid));
		JSONstring.append("\",");
		JSONstring.append("\"id\":\"");
		JSONstring.append(std::to_string(dataid));
		JSONstring.append("\"");
		JSONstring.append("}");
	}
	JSONstring.append("],");
	JSONstring.append("\"rootdirectory\":\"");
	JSONstring.append(dir);
	JSONstring.append("\"");
	

	delete mrcfile;
	
	return JSONstring;
	
}

std::string view2D (char* folder){
	std::string JSONstring;
	DIR 								*dir;
	struct dirent 						*entry;
	
	JSONstring.append("\"rootdirectory\":\"");
	JSONstring.append(folder);
	JSONstring.append("\",");
	
	JSONstring.append("\"iterations\":[ ");
	dir = opendir(folder);
	entry = readdir(dir);
	while (entry != NULL){
		if (strstr(entry->d_name, "classdoc") && strstr(entry->d_name, ".star") ) {
			JSONstring.append("\"");
			JSONstring.append(entry->d_name);
			JSONstring.append("\",");
		}
		entry = readdir(dir);
	}
	JSONstring.pop_back();
	JSONstring.append("]");
	return JSONstring;
}

std::string pipelineview (char* folder){
	std::string 						JSONstring;
	char								starfilename[2048];
	std::string							column;
	DIR 								*dir;
	struct dirent 						*entry;
	struct stat 						buffer;
	StarFile*    						starfile;
	std::vector<std::string>::iterator 	datait;
	std::vector<std::string>::iterator 	columnsit;
	bool								prefixsnapshotcomma;
	bool								prefixdatacomma;
	
	strcpy(starfilename, folder);
	strcat(starfilename, "/micrographs_preproc.star");
	JSONstring.append("\"rootdirectory\":\"");
	JSONstring.append(folder);
	JSONstring.append("\",");

	prefixsnapshotcomma = false;
	if(stat(starfilename, &buffer) == 0){
		starfile = new StarFile();
		readStarFile(starfile, starfilename);
		JSONstring.append("\"sortoptions\":[ ");
		for(columnsit = starfile->columns.begin() ; columnsit != starfile->columns.end(); ++columnsit){
			JSONstring.append("\"" + *columnsit + "\",");
		}
		JSONstring.pop_back();
		JSONstring.append("],");
		JSONstring.append("\"snapshots\":[ "); 
	/*	for(datait = starfile->data.begin() ; datait != starfile->data.end(); ++datait){
			if(prefixsnapshotcomma){
				JSONstring.append(",");
			} else {
				prefixsnapshotcomma = true;
			             }
			JSONstring.append("{");
	
			prefixdatacomma = false;
			for(unsigned int i = 0; i < starfile->columns.size(); i++){
				if(prefixdatacomma){
					JSONstring.append(",");
				} else {
					prefixdatacomma = true;
				}
				JSONstring.append("\"" + starfile->columns[i] + "\":\"" + getElementFromString(*datait, i) + "\"");
			}

			JSONstring.append("}");
		}*/
		for(unsigned int i = 0; i < starfile->data.size(); i++){
			if(prefixsnapshotcomma){
				JSONstring.append(",");
			} else {
				prefixsnapshotcomma = true;
			}
			JSONstring.append("{");
			JSONstring.append("\"id\":\"" + std::to_string(i) + "\"");
			for(unsigned int j = 0; j < starfile->columns.size(); j++){
				JSONstring.append(",");
				JSONstring.append("\"" + starfile->columns[j] + "\":\"" + getElementFromString(starfile->data[i], j) + "\"");
			}
			JSONstring.append("}");
		}
		delete starfile;
	}
	JSONstring.append("]");
	return JSONstring;
}

std::string particleView (char* particlestar, char* classnumber){
	std::string 						JSONstring;
	StarFile*							starfile;
	struct stat 						buffer;
	int									classcolumn;
	int									particlecolumn;
	std::string							testclassnumber;
	
	JSONstring.append("\"snapshots\":[ ");
	if(stat(particlestar, &buffer) == 0){
		starfile = new StarFile();
		readStarFile(starfile, particlestar);
		getStarColumn(starfile, "splReferenceImage", &particlecolumn);
		getStarColumn(starfile, "splClass", &classcolumn);
		for(unsigned int i = 0; i < starfile->data.size(); i++){
			testclassnumber = getElementFromString(starfile->data[i], classcolumn);
			std::cout << testclassnumber << " " << classnumber << std::endl;
			if(testclassnumber == classnumber){
				JSONstring.append("\"" + getElementFromString(starfile->data[i], particlecolumn) + "\",");
			}
		}
		delete starfile;
	}
	JSONstring.pop_back();
	JSONstring.append("] ");
	return JSONstring;
}

void starFileSaveSelection(char* inputfilename, char* outputfilename, char* inverseselection, std::string &JSONstring){
	char*							inverseselectionsplit;
	int								inverseselectionid;
	StarFile*    					starfile;

	starfile = new StarFile();
	readStarFile(starfile, inputfilename);
	inverseselectionsplit = std::strtok(inverseselection, ",");
	while (inverseselectionsplit != NULL) {
		inverseselectionid = std::stoi(inverseselectionsplit, nullptr);
		inverseselectionsplit = strtok (NULL, ",");
		starfile->data.erase(starfile->data.begin() + inverseselectionid);
	}
	writeStarFile(starfile, outputfilename);
	delete starfile;
	
	JSONstring = "{\"status\":\"success\"}";
}

void fileViewerInit (struct mg_connection* http_connection, struct http_message* message) {
 std::string  JSONstring;
 char         filename[1024];
 const char*	JSONchar;
 
 mg_get_http_var(&message->query_string, "filename", filename, sizeof(filename));
 JSONstring = "{";
 if (std::strstr(filename, ".star") != NULL) {
  JSONstring.append(starFileViewer(filename));
	}else if (std::strstr(filename, ".mrc") != NULL){
	JSONstring.append(mrcFileViewer(filename));	
 }else{
  JSONstring += "some:json";
 }
 JSONstring += "}";
	JSONchar = JSONstring.c_str();
	mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
	mg_send(http_connection, JSONchar, JSONstring.length());
}

void saveSelection (struct mg_connection* http_connection, struct http_message* message) {
	char         inputfilename[1024];
	char         outputfilename[1024];
	char         inverseselection[32768];
	std::string  JSONstring;
	const char*	JSONchar;
	 
	if (mg_get_http_var(&message->query_string, "inputfilename", inputfilename, sizeof(inputfilename)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "outputfilename", outputfilename, sizeof(outputfilename)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "inverseselection", inverseselection, sizeof(inverseselection)) <= 0){
		//error
	}
	if (std::strstr(inputfilename,".star") != NULL) {
		starFileSaveSelection(inputfilename, outputfilename, inverseselection, JSONstring);
	}
	JSONchar = JSONstring.c_str();
	mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
	mg_send(http_connection, JSONchar, JSONstring.length());

}

void save2DSelection (struct mg_connection* http_connection, struct http_message* message) {
	char         inputfilename[2048];
	char         outputfilename[2048];
	char         selection[32768];
	std::string  JSONstring;
	const char*	JSONchar;
	 
	if (mg_get_http_var(&message->query_string, "inputfilename", inputfilename, sizeof(inputfilename)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "outputfilename", outputfilename, sizeof(outputfilename)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "selection", selection, sizeof(selection)) <= 0){
		//error
	}
	if (std::strstr(inputfilename, ".mrc") != NULL) {
		selectMRCFile(inputfilename, outputfilename, selection);
	}
	JSONchar = JSONstring.c_str();
	mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
	mg_send(http_connection, JSONchar, JSONstring.length());

}

void save2DSelectionParticles (struct mg_connection* http_connection, struct http_message* message) {
	char         inputfilename[2048];
	char         outputfilename[2048];
	char         inverseselection[32768];
	std::string  JSONstring;
	const char*	 JSONchar;
	StarFile*		starfile;
	int			classcolumn;
	int			statecolumn;
	char*							inverseselectionsplit;
	std::string		classstring;


	if (mg_get_http_var(&message->query_string, "inputfilename", inputfilename, sizeof(inputfilename)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "outputfilename", outputfilename, sizeof(outputfilename)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "inverseselection", inverseselection, sizeof(inverseselection)) <= 0){
		//error
	}
	if (std::strstr(inputfilename, ".star") != NULL) {
		starfile = new StarFile();
		readStarFile(starfile, inputfilename);
		getStarColumn(starfile, "splClass", &classcolumn);
		getStarColumn(starfile, "splState", &statecolumn);
		inverseselectionsplit = std::strtok(inverseselection, ",");
		while (inverseselectionsplit != NULL) {
			for(unsigned int i = 0; i < starfile->data.size(); i++){
				classstring = getElementFromString(starfile->data[i], classcolumn);
				if(classstring == inverseselectionsplit){
					replaceStarColumn(starfile, i, statecolumn, "0");
				}
			}
			inverseselectionsplit = strtok (NULL, ",");
		}
		writeStarFile(starfile, outputfilename);
		delete starfile;
	}
	JSONchar = JSONstring.c_str();
	mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
	mg_send(http_connection, JSONchar, JSONstring.length());

}

void view2DInit (struct mg_connection* http_connection, struct http_message* message) {
	std::string  JSONstring;
	char         folder[2048];
	const char*	JSONchar;
 
	mg_get_http_var(&message->query_string, "folder", folder, sizeof(folder));
	JSONstring = "{";
	JSONstring.append(view2D(folder));
	JSONstring += "}";
	JSONchar = JSONstring.c_str();
	mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
	mg_send(http_connection, JSONchar, JSONstring.length());
}

void pipelineviewInit (struct mg_connection* http_connection, struct http_message* message) {
	std::string  JSONstring;
	char         folder[2048];
	const char*	JSONchar;
 
	mg_get_http_var(&message->query_string, "folder", folder, sizeof(folder));
	JSONstring = "{";
	JSONstring.append(pipelineview(folder));
	JSONstring += "}";
	JSONchar = JSONstring.c_str();
	mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
	mg_send(http_connection, JSONchar, JSONstring.length());
}

void particleViewInit (struct mg_connection* http_connection, struct http_message* message) {
	std::string  JSONstring;
	char         particlestar[2048];
	char         classnumber[32];
	const char*	JSONchar;
 
	mg_get_http_var(&message->query_string, "particlestar", particlestar, sizeof(particlestar));
	mg_get_http_var(&message->query_string, "class", classnumber, sizeof(classnumber));
	JSONstring = "{";
	JSONstring.append(particleView(particlestar, classnumber));
	JSONstring += "}";
	JSONchar = JSONstring.c_str();
	mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
	mg_send(http_connection, JSONchar, JSONstring.length());
}

void boxFileData (struct mg_connection* http_connection, struct http_message* message) {
	std::string  JSONstring;
	char         filename[2048];
	const char*	 JSONchar;
	std::ifstream input;
	std::string line;
 
	mg_get_http_var(&message->query_string, "filename", filename, sizeof(filename));
	JSONstring = "{ \"boxes\": [ ";
	
	input.open(filename);
	if(input.is_open()){
		
		while(getline(input, line)){
			JSONstring += "[";
			JSONstring += "\"" + getElementFromString(line, 0) +"\",";
			JSONstring += "\"" + getElementFromString(line, 1) +"\",";
			JSONstring += "\"" + getElementFromString(line, 2) +"\"";
			JSONstring += "],";
		}
		input.close();
		JSONstring.pop_back();
	} else {
		std::cout << "no boxfile" << std::endl;
	}
	
	JSONstring += "]}";
	JSONchar = JSONstring.c_str();
	mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
	mg_send(http_connection, JSONchar, JSONstring.length());
}

#endif
