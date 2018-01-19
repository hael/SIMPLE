#ifndef SIMPLE_GUI_MIDRUN_H_INCLUDED
#define SIMPLE_GUI_MIDRUN_H_INCLUDED

#include "simple_gui_postrun.h"

void watchRunPreproc (struct mg_connection *nc, std::string &messagedata, int iterations) {

	int									iteration;
	DIR 								*dir;
	struct dirent 						*entry;
	bool								processed;
	std::string							unidocfilename;
	std::vector<std::string>			processedunidocs;
	std::vector<std::string>::iterator	it;
	
	StarFile* starfile = new StarFile();
	
	iteration = 0;
	
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
	
	if (iterations == 0){
		iterations = 100;
	}
	
	while (iteration <= iterations){
		dir = opendir("pipeline");
		entry = readdir(dir);
		while (entry != NULL){
			unidocfilename = "pipeline/";
			unidocfilename += entry->d_name;
			if (unidocfilename.find("unidoc") != std::string::npos) {
				processed = false;
				for (it = processedunidocs.begin() ; it != processedunidocs.end(); ++it) {
					if (*it == unidocfilename) {
						processed = true;
						break;
					}
				}
				if (!processed) {
					UniDoc* unidoc = new UniDoc();
					readUniDoc(unidoc, unidocfilename);
					starfile->data.push_back(unidoc->data[0]);
					delete unidoc;
					
				}
			}
			entry = readdir(dir);
		}
		closedir(dir);
		writeStarFile(starfile, "micrographs_preproc.star");
		usleep(6000000);
		iteration++;	
	}
}

#endif
