#ifndef SIMPLE_GUI_POSTRUN_H_INCLUDED
#define SIMPLE_GUI_POSTRUN_H_INCLUDED

#include "simple_gui_textdata.h"
#include "simple_gui_stardata.h"
#include "simple_gui_unidocdata.h"

void postRunPreproc (struct mg_connection *nc, std::string &messagedata, std::string &jobdirectory) {

	DIR 								*dir;
	struct dirent 						*entry;
	bool								processed;
	std::string							unidocfilename;
	std::string							newstarfileline;
	std::string							unidocelement;
	bool								starfilelinecorrect;
	
	StarFile* starfile = new StarFile();
	
	starfile->columns.push_back("splDirectory");
	starfile->columns.push_back("splMovie");
	starfile->columns.push_back("splMicrographName");
	starfile->columns.push_back("splForCTFName");
	starfile->columns.push_back("splThumbName");
	starfile->columns.push_back("splPspecName");
	starfile->columns.push_back("splSamplingDistance");
	//starfile->columns.push_back("rlnMagnification");
	//starfile->columns.push_back("rlnDetectorPixelSize");
	starfile->columns.push_back("splVoltage");
	starfile->columns.push_back("splDefocusX");
	starfile->columns.push_back("splDefocusY");
	//starfile->columns.push_back("rlnDefocusU");
	//starfile->columns.push_back("rlnDefocusV");
	starfile->columns.push_back("splSphericalAberration");
	starfile->columns.push_back("splDefocusAngle");
	starfile->columns.push_back("splCtfMaxResolution");
	starfile->columns.push_back("splAmplitudeContrast");
	
	dir = opendir("pipeline");
	entry = readdir(dir);
	while (entry != NULL){
		unidocfilename = "pipeline/";
		unidocfilename += entry->d_name;
		if (unidocfilename.find("unidoc") != std::string::npos) {
			UniDoc* unidoc = new UniDoc();
			readUniDoc(unidoc, unidocfilename);
			newstarfileline = "";
			starfilelinecorrect = true;
			if (jobdirectory != ""){
				newstarfileline += jobdirectory + " ";
			}else{
				starfilelinecorrect = false;
			}
			getUniDocLineValue(unidoc, 0, "movie", unidocelement);
			if (unidocelement != ""){
				newstarfileline +=  unidocelement + " ";
			}else{
				starfilelinecorrect = false;
			}
			getUniDocLineValue(unidoc, 0, "intg", unidocelement);
			if (unidocelement != ""){
				newstarfileline +=  unidocelement + " ";
			}else{
				starfilelinecorrect = false;
			}
			getUniDocLineValue(unidoc, 0, "forctf", unidocelement);
			if (unidocelement != ""){
				newstarfileline +=  unidocelement + " ";
			}else{
				starfilelinecorrect = false;
			}
			getUniDocLineValue(unidoc, 0, "thumb", unidocelement);
			if (unidocelement != ""){
				newstarfileline +=  unidocelement + " ";
			}else{
				starfilelinecorrect = false;
			}
			getUniDocLineValue(unidoc, 0, "pspec", unidocelement);
			if (unidocelement != ""){
				newstarfileline +=  unidocelement + " ";
			}else{
				starfilelinecorrect = false;
			}
			getUniDocLineValue(unidoc, 0, "smpd", unidocelement);
			if (unidocelement != ""){
				newstarfileline +=  unidocelement + " ";
			}else{
				starfilelinecorrect = false;
			}
			getUniDocLineValue(unidoc, 0, "kv", unidocelement);
			if (unidocelement != ""){
				newstarfileline +=  unidocelement + " ";
			}else{
				starfilelinecorrect = false;
			}
			getUniDocLineValue(unidoc, 0, "dfx", unidocelement);
			if (unidocelement != ""){
				newstarfileline +=  unidocelement + " ";
			}else{
				starfilelinecorrect = false;
			}
			getUniDocLineValue(unidoc, 0, "dfy", unidocelement);
			if (unidocelement != ""){
				newstarfileline +=  unidocelement + " ";
			}else{
				starfilelinecorrect = false;
			}
			getUniDocLineValue(unidoc, 0, "cs", unidocelement);
			if (unidocelement != ""){
				newstarfileline +=  unidocelement + " ";
			}else{
				starfilelinecorrect = false;
			}
			getUniDocLineValue(unidoc, 0, "angast", unidocelement);
			if (unidocelement != ""){
				newstarfileline +=  unidocelement + " ";
			}else{
				starfilelinecorrect = false;
			}
			getUniDocLineValue(unidoc, 0, "ctfres", unidocelement);
			if (unidocelement != ""){
				newstarfileline +=  unidocelement + " ";
			}else{
				starfilelinecorrect = false;
			}
			getUniDocLineValue(unidoc, 0, "fraca", unidocelement);
			if (unidocelement != ""){
				newstarfileline +=  unidocelement + " ";
			}else{
				starfilelinecorrect = false;
			}
			delete unidoc;
			if (starfilelinecorrect){
				starfile->data.push_back(newstarfileline);
				std::cout << newstarfileline << std::endl;
			}
		}
		entry = readdir(dir);
	}
	closedir(dir);
	writeStarFile(starfile, "micrographs_preproc.star");
	delete starfile;	
}

#endif
