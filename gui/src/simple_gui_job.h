#ifndef SIMPLE_GUI_JOB_H_INCLUDED
#define SIMPLE_GUI_JOB_H_INCLUDED

#include "simple_gui_stardata.h"
#include "simple_gui_textdata.h"
#include "simple_gui_unidocdata.h"
#include "simple_gui_filesystem.h"
#include "simple_gui_sql.h"

#include <unistd.h>
#include <cstring>
#include <string>
#include <libgen.h>

extern "C" {
    #include "ext/mongoose/mongoose.h"
}


void addJob (struct mg_connection* http_connection, struct http_message* message, int &jobid) {
	std::string  JSONstring;
	std::string  sqlresult;
	const char*	 JSONchar;
	char		 databasepath[1024];
	char         table[128];
	char         jobname[128];
	char         jobdescription[2048];
	char		 jobfolder[1024];
	char		 jobtype[64];
	
	if (mg_get_http_var(&message->query_string, "table", table, sizeof(table)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "jobname", jobname, sizeof(jobname)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "jobdescription", jobdescription, sizeof(jobdescription)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "jobfolder", jobfolder, sizeof(jobfolder)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "jobtype", jobtype, sizeof(jobtype)) <= 0){
		//error
	}
	
	
	
	databaseLocation(message, databasepath);
	
	SQLQuery(databasepath, JSONstring, "INSERT INTO " + std::string(table) + "(\
										jobname,\
										jobdescription,\
										jobfolder,\
										jobtype,\
										jobstatus,\
										jobpid)\
										VALUES('"
										+ std::string(jobname) + "','"
										+ std::string(jobdescription) + "','"
										+ std::string(jobfolder) + "','"
										+ std::string(jobtype) + "','"
										+ "External" + "','"
										+ std::to_string(0) + "');", jobid);
}

void updateStatus (int &jobid, char* table, char* databasepath, std::string status) {
	std::string 	JSONstring;
	int				buffer;
	SQLQuery(databasepath, JSONstring, "UPDATE " + std::string(table) + " SET  jobstatus = \'" + status + "\' WHERE jobid=" + std::to_string(jobid) +";" , buffer);
}

void preprocPost (char* directory){
	//MAY NEED TO MODIFY TO UPDATE RATHER THAN RECREATE FOR SELECTION DURING PIPELINE
	DIR 								*dir;
	struct dirent 						*entry;
	std::string							folder;
	std::string							unidocfilename;
	std::string							newstarfileline;
	std::string							unidocelement;
	std::string							newstarfilename;
	bool								starfilelinecorrect;
	
	StarFile* starfile = new StarFile();
	starfile->columns.push_back("splMovie");
	starfile->columns.push_back("splMicrographName");
	starfile->columns.push_back("splBoxFileName");
	starfile->columns.push_back("splForCTFName");
	starfile->columns.push_back("splCTFFINDFit");
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
	
	folder = directory;
	folder += "/output/";
	dir = opendir(folder.c_str());
	entry = readdir(dir);
	
	while (entry != NULL){
		unidocfilename = folder;
		unidocfilename += entry->d_name;
		if (unidocfilename.find("unidoc") != std::string::npos) {
			UniDoc* unidoc = new UniDoc();
			readUniDoc(unidoc, unidocfilename);
			newstarfileline = "";
			starfilelinecorrect = true;
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
			unidocelement.replace(unidocelement.begin() + unidocelement.find(".mrc"), unidocelement.end(), ".box");
			unidocelement.replace(unidocelement.begin(), unidocelement.begin() + 7, "");
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
			unidocelement.replace(unidocelement.begin() + unidocelement.find(".mrc"), unidocelement.end(), "_ctffind_diag.mrc");
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
			}
		}
		entry = readdir(dir);
	}
	closedir(dir);
	newstarfilename = directory;
	newstarfilename += "/micrographs_preproc.star";
	writeStarFile(starfile, newstarfilename.c_str());
	delete starfile;	
}

void prime2DPostOld (char* directory){
	DIR 								*dir;
	struct dirent 						*entry;
	std::string							classdocfilename;
	std::string							cavgsfilename;
	std::string							primedocfilename;
	std::string							newstarfileline;
	std::string							unidocelement;
	std::string							newstarfilename;
	int									unidoclinecount;
	bool								starfilelinecorrect;
	
	StarFile* starfile = new StarFile();
	starfile->columns.push_back("splReferenceImage");
	starfile->columns.push_back("splParticleStar");
	starfile->columns.push_back("splClassPopulation");
	starfile->columns.push_back("splClassResolution");
	starfile->columns.push_back("splClassCorrelation");
	starfile->columns.push_back("splClassWeight");

	classdocfilename = directory;
	classdocfilename += "/classdoc.txt";
	cavgsfilename = directory;
	cavgsfilename += "/cavgs_final.mrc";
	primedocfilename = directory;
	primedocfilename += "/prime2Ddoc_final.txt";
	newstarfilename = directory;
	newstarfilename += "/cavgs_final.star";
	
	UniDoc* unidoc = new UniDoc();
	readUniDoc(unidoc, classdocfilename);
	
	for(unidoclinecount = 0; unidoclinecount < unidoc->data.size(); unidoclinecount++){
		newstarfileline = std::to_string(unidoclinecount + 1);
		newstarfileline += "@cavgs_final.mrc prime2Ddoc_final.star "; 
		starfilelinecorrect = true;
		
		getUniDocLineValue(unidoc, unidoclinecount, "class", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		getUniDocLineValue(unidoc, unidoclinecount, "pop", unidocelement);
                if (unidocelement != ""){
                        newstarfileline +=  unidocelement + " ";
                }else{
                        starfilelinecorrect = false;
                }
		getUniDocLineValue(unidoc, unidoclinecount, "res", unidocelement);
                if (unidocelement != ""){
                        newstarfileline +=  unidocelement + " ";
                }else{
                        starfilelinecorrect = false;
                }
		getUniDocLineValue(unidoc, unidoclinecount, "corr", unidocelement);
                if (unidocelement != ""){
                        newstarfileline +=  unidocelement + " ";
                }else{
                        starfilelinecorrect = false;
                }
		getUniDocLineValue(unidoc, unidoclinecount, "w", unidocelement);
                if (unidocelement != ""){
                        newstarfileline +=  unidocelement + " ";
                }else{
                        starfilelinecorrect = false;
                }
		if (starfilelinecorrect){
			starfile->data.push_back(newstarfileline);
		}
	}
	
	delete unidoc;
	writeStarFile(starfile, newstarfilename.c_str());
	delete starfile;
	
	newstarfilename = directory;
	newstarfilename += "/prime2Ddoc_final.star";
	
	starfile = new StarFile();
	starfile->columns.push_back("splReferenceImage");
	starfile->columns.push_back("splE1");
	starfile->columns.push_back("splE2");
	starfile->columns.push_back("splE3");
	starfile->columns.push_back("splX");
	starfile->columns.push_back("splY");
	starfile->columns.push_back("splDist");
	starfile->columns.push_back("splState");
	starfile->columns.push_back("splStateBalance");
	starfile->columns.push_back("splFrac");
	starfile->columns.push_back("splEo");
	starfile->columns.push_back("splSmpd");
	starfile->columns.push_back("splKv");
	starfile->columns.push_back("splCs");
	starfile->columns.push_back("splFraca");
	starfile->columns.push_back("splDfx");
	starfile->columns.push_back("splDfy");
	starfile->columns.push_back("splAngast");
	starfile->columns.push_back("splW");
	starfile->columns.push_back("splLp");
	starfile->columns.push_back("splClass");
	starfile->columns.push_back("splCorr");
	starfile->columns.push_back("splSpecscore");
	starfile->columns.push_back("splDistInplane");
	starfile->columns.push_back("splMiClass");
	starfile->columns.push_back("splMiInpl");
	starfile->columns.push_back("splMiJoint");

	unidoc = new UniDoc();
	readUniDoc(unidoc, primedocfilename);
	
	for(unidoclinecount = 0; unidoclinecount < unidoc->data.size(); unidoclinecount++){
		newstarfileline = std::to_string(unidoclinecount + 1);
		newstarfileline += "@sumstack.mrc"; 
		starfilelinecorrect = true;
		
		getUniDocLineValue(unidoc, unidoclinecount, "e1", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "e2", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "e3", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "x", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "y", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "dist", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "state", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "state_balance", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "frac", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "eo", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "smpd", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "kv", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "cs", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "fraca", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "dfx", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "dfy", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "angast", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "w", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "lp", unidocelement);
        if (unidocelement != ""){
            newstarfileline +=  unidocelement + " ";
        }else{
            starfilelinecorrect = false;
        }
        
        getUniDocLineValue(unidoc, unidoclinecount, "class", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "corr", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "specscore", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "dist_inpl", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "mi_class", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "mi_inpl", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		getUniDocLineValue(unidoc, unidoclinecount, "mi_joint", unidocelement);
		if (unidocelement != ""){
			newstarfileline +=  unidocelement + " ";
		}else{
			starfilelinecorrect = false;
		}
		
		if (starfilelinecorrect){
			starfile->data.push_back(newstarfileline);
		}
	}
	
	delete unidoc;
	writeStarFile(starfile, newstarfilename.c_str());
	delete starfile;	
}

void prime2DPost (char* directory){
	std::cout << "prime2DPost" << std::endl;
	DIR 								*dir;
	struct dirent 						*entry;
	std::string							classdocfilename;
	std::string							cavgsfilename;
	std::string							primedocfilename;
	std::string							newstarfileline;
	std::string							unidocelement;
	std::string							classdocstarfilename;	
	std::string							primedocstarfilename;	
	std::string							cavgsfilenamebase;
	std::string							primedocstarfilenamebase;
	std::string							stacktab;
	int									unidoclinecount;
	bool								starfilelinecorrect;
	struct stat							buffer;
	UniDoc* 							unidoc;
	StarFile* 							starfile;
	UniDoc*								stacktable;
	
	dir = opendir(directory);
	entry = readdir(dir);
	
	while (entry != NULL){
		if (strstr(entry->d_name, "cavgs_iter") && strstr(entry->d_name, "_even") == NULL && strstr(entry->d_name, "_odd") == NULL) {
			cavgsfilename = std::string(entry->d_name);
			cavgsfilenamebase = cavgsfilename;
			primedocfilename = cavgsfilename;
			primedocfilename.replace(0, 10, "prime2Ddoc_");
			primedocfilename.replace(primedocfilename.end()-4, primedocfilename.end(), ".txt");
			primedocstarfilename = primedocfilename;
			primedocstarfilename.replace(primedocstarfilename.end()-4, primedocstarfilename.end(), ".star");
			primedocstarfilenamebase = primedocstarfilename;
			classdocfilename = cavgsfilename;
			classdocfilename.replace(0, 10, "classdoc_");
			classdocfilename.replace(classdocfilename.end()-4, classdocfilename.end(), ".txt");
			classdocstarfilename = classdocfilename;
			classdocstarfilename.replace(classdocstarfilename.end()-4, classdocstarfilename.end(), ".star");
			cavgsfilename.replace(0, 0, "/");
			cavgsfilename.replace(0, 0, directory);
			primedocfilename.replace(0, 0, "/");
			primedocfilename.replace(0, 0, directory);
			classdocfilename.replace(0, 0, "/");
			classdocfilename.replace(0, 0, directory);
			primedocstarfilename.replace(0, 0, "/");
			primedocstarfilename.replace(0, 0, directory);
			classdocstarfilename.replace(0, 0, "/");
			classdocstarfilename.replace(0, 0, directory);
			stacktab = std::string(directory) + "/stktab_info.txt";
			if(stat(classdocstarfilename.c_str(), &buffer) != 0 && stat(classdocfilename.c_str(), &buffer) == 0){
				std::cout << "processing " << std::endl;
				starfile = new StarFile();
				starfile->columns.push_back("splReferenceImage");
				starfile->columns.push_back("splParticleStar");
				starfile->columns.push_back("splClassNumber");
				starfile->columns.push_back("splClassPopulation");
				starfile->columns.push_back("splClassResolution");
				starfile->columns.push_back("splClassCorrelation");
				starfile->columns.push_back("splClassWeight");
		
				unidoc = new UniDoc();
				readUniDoc(unidoc, classdocfilename);
		
				for(unidoclinecount = 0; unidoclinecount < unidoc->data.size(); unidoclinecount++){
					newstarfileline = std::to_string(unidoclinecount + 1);
					newstarfileline += "@" + cavgsfilenamebase + " " + primedocstarfilenamebase + " ";
					starfilelinecorrect = true;
					
					getUniDocLineValue(unidoc, unidoclinecount, "class", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					getUniDocLineValue(unidoc, unidoclinecount, "pop", unidocelement);
							if (unidocelement != ""){
									newstarfileline +=  unidocelement + " ";
							}else{
									starfilelinecorrect = false;
							}
					getUniDocLineValue(unidoc, unidoclinecount, "res", unidocelement);
							if (unidocelement != ""){
									newstarfileline +=  unidocelement + " ";
							}else{
									starfilelinecorrect = false;
							}
					getUniDocLineValue(unidoc, unidoclinecount, "corr", unidocelement);
							if (unidocelement != ""){
									newstarfileline +=  unidocelement + " ";
							}else{
									starfilelinecorrect = false;
							}
					getUniDocLineValue(unidoc, unidoclinecount, "w", unidocelement);
							if (unidocelement != ""){
									newstarfileline +=  unidocelement + " ";
							}else{
									starfilelinecorrect = false;
							}
					if (starfilelinecorrect){
						starfile->data.push_back(newstarfileline);
					}
				}
				delete unidoc;
				writeStarFile(starfile, classdocstarfilename.c_str());
				delete starfile;
			}
			
			if(stat(stacktab.c_str(), &buffer) == 0){
				stacktable = new UniDoc();
				readStackTab(stacktable, stacktab);
			}
			
			if(stat(primedocstarfilename.c_str(), &buffer) != 0){ 
				starfile = new StarFile();
				starfile->columns.push_back("splReferenceImage");
				starfile->columns.push_back("splE1");
				starfile->columns.push_back("splE2");
				starfile->columns.push_back("splE3");
				starfile->columns.push_back("splX");
				starfile->columns.push_back("splY");
				starfile->columns.push_back("splDist");
				starfile->columns.push_back("splState");
				starfile->columns.push_back("splStateBalance");
				starfile->columns.push_back("splFrac");
				starfile->columns.push_back("splEo");
				starfile->columns.push_back("splSmpd");
				starfile->columns.push_back("splKv");
				starfile->columns.push_back("splCs");
				starfile->columns.push_back("splFraca");
				starfile->columns.push_back("splDfx");
				starfile->columns.push_back("splDfy");
				starfile->columns.push_back("splAngast");
				starfile->columns.push_back("splW");
				starfile->columns.push_back("splLp");
				starfile->columns.push_back("splClass");
				starfile->columns.push_back("splCorr");
				starfile->columns.push_back("splSpecscore");
				starfile->columns.push_back("splDistInplane");
				starfile->columns.push_back("splMiClass");
				starfile->columns.push_back("splMiInpl");
				starfile->columns.push_back("splMiJoint");

				unidoc = new UniDoc();
				readUniDoc(unidoc, primedocfilename);
				
				for(unidoclinecount = 0; unidoclinecount < unidoc->data.size(); unidoclinecount++){
					//newstarfileline = std::to_string(unidoclinecount + 1);
					//newstarfileline += "@sumstack.mrc "; 
					newstarfileline = getElementFromString(stacktable->data[unidoclinecount], 0) + "@" + getElementFromString(stacktable->data[unidoclinecount], 1) + " ";
					starfilelinecorrect = true;
					
					getUniDocLineValue(unidoc, unidoclinecount, "e1", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "e2", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "e3", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "x", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "y", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "dist", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "state", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "state_balance", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "frac", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "eo", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "smpd", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "kv", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "cs", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "fraca", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "dfx", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "dfy", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "angast", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "w", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "lp", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "class", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "corr", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "specscore", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "dist_inpl", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "mi_class", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "mi_inpl", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					getUniDocLineValue(unidoc, unidoclinecount, "mi_joint", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					
					if (starfilelinecorrect){
						starfile->data.push_back(newstarfileline);
					}
				}
				
				delete unidoc;
				writeStarFile(starfile, primedocstarfilename.c_str());
				delete starfile;
				delete stacktable;
			}
		}
		entry = readdir(dir);
	}
}

void streamPost (char* directory){
	DIR 								*dir;
	struct dirent 						*entry;
	std::string							classdocfilename;
	std::string							cavgsfilename;
	std::string							primedocfilename;
	std::string							newstarfileline;
	std::string							unidocelement;
	std::string							classdocstarfilename;	
	std::string							primedocstarfilename;	
	std::string							cavgsfilenamebase;
	std::string							primedocstarfilenamebase;
	int									unidoclinecount;
	bool								starfilelinecorrect;
	UniDoc* 							unidoc;
	StarFile* 							starfile;
	struct stat							buffer;
	
	dir = opendir(directory);
	entry = readdir(dir);
	
	while (entry != NULL){
		if (strstr(entry->d_name, "cavgs_iter")&& strstr(entry->d_name, "_even") == NULL && strstr(entry->d_name, "_odd") == NULL) {
			cavgsfilename = entry->d_name;
			cavgsfilenamebase = cavgsfilename;
			primedocfilename = cavgsfilename;
			primedocfilename.replace(0, 10, "prime2Ddoc_");
			primedocfilename.replace(primedocfilename.end()-4, primedocfilename.end(), ".txt");
			primedocstarfilename = primedocfilename;
			primedocstarfilename.replace(primedocstarfilename.end()-4, primedocstarfilename.end(), ".star");
			primedocstarfilenamebase = primedocstarfilename;
			classdocfilename = cavgsfilename;
			classdocfilename.replace(0, 10, "classdoc_");
			classdocfilename.replace(classdocfilename.end()-4, classdocfilename.end(), ".txt");
			classdocstarfilename = classdocfilename;
			classdocstarfilename.replace(classdocstarfilename.end()-4, classdocstarfilename.end(), ".star");
			cavgsfilename.replace(0, 0, "/");
			cavgsfilename.replace(0, 0, directory);
			primedocfilename.replace(0, 0, "/");
			primedocfilename.replace(0, 0, directory);
			classdocfilename.replace(0, 0, "/");
			classdocfilename.replace(0, 0, directory);
			primedocstarfilename.replace(0, 0, "/");
			primedocstarfilename.replace(0, 0, directory);
			classdocstarfilename.replace(0, 0, "/");
			classdocstarfilename.replace(0, 0, directory);
			
			if(stat(classdocstarfilename.c_str(), &buffer) != 0 && stat(classdocfilename.c_str(), &buffer) == 0){
				starfile = new StarFile();
				starfile->columns.push_back("splReferenceImage");
				starfile->columns.push_back("splParticleStar");
				starfile->columns.push_back("splClassNumber");
				starfile->columns.push_back("splClassPopulation");
				starfile->columns.push_back("splClassResolution");
				starfile->columns.push_back("splClassCorrelation");
				starfile->columns.push_back("splClassWeight");
		
				unidoc = new UniDoc();
				readUniDoc(unidoc, classdocfilename);
				std::cout << classdocfilename << std::endl;
		
				for(unidoclinecount = 0; unidoclinecount < unidoc->data.size(); unidoclinecount++){
					newstarfileline = std::to_string(unidoclinecount + 1);
					newstarfileline += "@" + cavgsfilenamebase + " " + primedocstarfilenamebase + " ";
					starfilelinecorrect = true;
					
					getUniDocLineValue(unidoc, unidoclinecount, "class", unidocelement);
					if (unidocelement != ""){
						newstarfileline +=  unidocelement + " ";
					}else{
						starfilelinecorrect = false;
					}
					getUniDocLineValue(unidoc, unidoclinecount, "pop", unidocelement);
							if (unidocelement != ""){
									newstarfileline +=  unidocelement + " ";
							}else{
									starfilelinecorrect = false;
							}
					getUniDocLineValue(unidoc, unidoclinecount, "res", unidocelement);
							if (unidocelement != ""){
									newstarfileline +=  unidocelement + " ";
							}else{
									starfilelinecorrect = false;
							}
					getUniDocLineValue(unidoc, unidoclinecount, "corr", unidocelement);
							if (unidocelement != ""){
									newstarfileline +=  unidocelement + " ";
							}else{
									starfilelinecorrect = false;
							}
					getUniDocLineValue(unidoc, unidoclinecount, "w", unidocelement);
							if (unidocelement != ""){
									newstarfileline +=  unidocelement + " ";
							}else{
									starfilelinecorrect = false;
							}
					if (starfilelinecorrect){
						starfile->data.push_back(newstarfileline);
					}
				}
				
				delete unidoc;
				writeStarFile(starfile, classdocstarfilename.c_str());
				delete starfile;
			}
			
			/*starfile = new StarFile();
			starfile->columns.push_back("splReferenceImage");
			starfile->columns.push_back("splE1");
			starfile->columns.push_back("splE2");
			starfile->columns.push_back("splE3");
			starfile->columns.push_back("splX");
			starfile->columns.push_back("splY");
			starfile->columns.push_back("splDist");
			starfile->columns.push_back("splState");
			starfile->columns.push_back("splStateBalance");
			starfile->columns.push_back("splFrac");
			starfile->columns.push_back("splEo");
			starfile->columns.push_back("splSmpd");
			starfile->columns.push_back("splKv");
			starfile->columns.push_back("splCs");
			starfile->columns.push_back("splFraca");
			starfile->columns.push_back("splDfx");
			starfile->columns.push_back("splDfy");
			starfile->columns.push_back("splAngast");
			starfile->columns.push_back("splW");
			starfile->columns.push_back("splLp");
			starfile->columns.push_back("splClass");
			starfile->columns.push_back("splCorr");
			starfile->columns.push_back("splSpecscore");
			starfile->columns.push_back("splDistInplane");
			starfile->columns.push_back("splMiClass");
			starfile->columns.push_back("splMiInpl");
			starfile->columns.push_back("splMiJoint");

			unidoc = new UniDoc();
			readUniDoc(unidoc, primedocfilename);
			
			for(unidoclinecount = 0; unidoclinecount < unidoc->data.size(); unidoclinecount++){
				newstarfileline = std::to_string(unidoclinecount + 1);
				newstarfileline += "@sumstack.mrc "; 
				starfilelinecorrect = true;
				
				getUniDocLineValue(unidoc, unidoclinecount, "e1", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "e2", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "e3", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "x", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "y", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "dist", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "state", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "state_balance", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "frac", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "eo", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "smpd", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "kv", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "cs", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "fraca", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "dfx", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "dfy", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "angast", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "w", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "lp", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "class", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "corr", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "specscore", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "dist_inpl", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "mi_class", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "mi_inpl", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				getUniDocLineValue(unidoc, unidoclinecount, "mi_joint", unidocelement);
				if (unidocelement != ""){
					newstarfileline +=  unidocelement + " ";
				}else{
					starfilelinecorrect = false;
				}
				
				if (starfilelinecorrect){
					starfile->data.push_back(newstarfileline);
				}
			}
			
			delete unidoc;
			writeStarFile(starfile, primedocstarfilename.c_str());
			delete starfile;*/
		}
		entry = readdir(dir);
	}

}

void externalJob (struct mg_connection* http_connection, struct http_message* message) {
	char         jobname[128];
	char         jobdescription[2048];
	char         jobfolder[2048];
	char         jobtype[32];
	std::string  JSONstring;
	const char*	JSONchar;
	char		databasepath[1024];
	pid_t 		pid;
	int			jobid;
	char         table[128];
		
	if (mg_get_http_var(&message->query_string, "jobname", jobname, sizeof(jobname)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "jobdescription", jobdescription, sizeof(jobdescription)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "jobfolder", jobfolder, sizeof(jobfolder)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "jobtype", jobtype, sizeof(jobtype)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "table", table, sizeof(table)) <= 0){
		//error
	}
	
	databaseLocation(message, databasepath);
	addJob(http_connection, message, jobid);
	
	pid = fork();
	
	if (pid == 0){
		updateStatus (jobid,table, databasepath,"Running");
		if (std::strstr(jobtype, "preproc") ) {
			preprocPost (jobfolder);
		}else if (std::strstr(jobtype, "prime2d") ) {
			prime2DPost (jobfolder);
		}else if (std::strstr(jobtype, "stream") ) {
			streamPost (jobfolder);
		}
		updateStatus (jobid,table, databasepath,"External");
	} else if (pid > 0) {
		JSONchar = JSONstring.c_str();
		mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
		mg_send(http_connection, JSONchar, JSONstring.length());
	} else {
		std::cout << "Fork Failed" << std::endl;
	}
}

void getJobs (struct mg_connection* http_connection, struct http_message* message) {
	std::string  JSONstring;
	std::string  sqlresult;
	const char*	 JSONchar;
	char		 databasepath[1024];
	char         table[128];
	int			jobid;
	
	if (mg_get_http_var(&message->query_string, "table", table, sizeof(table)) <= 0){
		//error
	}
	
	databaseLocation(message, databasepath);
	
	SQLQuery(databasepath, JSONstring, "CREATE TABLE IF NOT EXISTS " + std::string(table) + " (\
										jobid integer PRIMARY KEY,\
										jobname text NOT NULL,\
										jobdescription text NOT NULL,\
										jobfolder text NOT NULL,\
										jobtype text NOT NULL,\
										jobstatus text NOT NULL,\
										jobpid integer NOT NULL);", jobid);
										
	SQLQuery(databasepath, sqlresult, "SELECT * FROM " + std::string(table), jobid);	
								
	JSONstring = "{\"jobs\": [";
	JSONstring += sqlresult;
	JSONstring += "]}";
	JSONchar = JSONstring.c_str();

	mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
	mg_send(http_connection, JSONchar, JSONstring.length());	
}



#endif
