#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstring>
#include <map>
#include <sqlite3.h> 
#include <pwd.h>
#include <map>
#include <sys/wait.h>
#include <jpeglib.h>
#include <libgen.h>
#include <algorithm>
#include <signal.h>
#include <cstdio>

extern "C" {
    #include "../ext/mongoose/mongoose.h"
}

struct mg_serve_http_opts 	http_server_opts;

struct JPEGResponse {
	
	unsigned char* 			pixels;
	unsigned char* 			jpeg;
	unsigned long			jpegsize;
	int						xdim;
	int						ydim;
	
};

struct JSONResponse {

	std::vector< std::map<std::string, std::string> >			snapshots;
	std::vector< std::map<std::string, std::string> >			projects;
	std::vector< std::map<std::string, std::string> >			jobs;
	std::vector< std::map<std::string, std::string> >			iterations;
	std::vector< std::map<std::string, std::string> >			boxes;
	std::vector<std::string>									files;
	std::vector<std::string>									directories;
	std::string													rootdirectory;
	std::string													prime2ddoc;
	std::string													inputfilename;
	std::string													logfile;
	std::string													error;
	std::vector< std::map<std::string, std::string> >			rerunjob;
	std::string													JSONstring;
	std::string													particlecount;
	std::string													jobfolder;
};

struct UniDoc {
	
	std::vector<std::string> data;
	
};

struct MRCFile {

	float*							data;
	char							header[1024];
	int 							nx;
	int								ny;
	int								nz;
	int								mode;
	int								extendedheader;
};

bool getRequestVariable (struct http_message* message, std::string variable, std::string& value){
	
	char		buffer[524288];
	int			returncode;
	
	returncode = mg_get_http_var(&message->query_string, variable.c_str(), buffer, sizeof(buffer));
	
	if (returncode > 0){
		
		value = std::string(buffer);
		return true;
		
	} else {
		
		return false;
	
	}
}

void getElementFromString(std::string argstring, std::string& returnstring, int element){ 
	std::size_t 		location;
	char 				character;
	int 				currentelement;
	int 				nextelement;
	
	returnstring.clear();

	location = argstring.find("  ");
	while (location != std::string::npos) {
		argstring.replace (location, 2, " ");
		location = argstring.find("  ");
	}

	location = 0;
	
	while(argstring.c_str()[0] == ' '){
		argstring.erase (0, 1);
	}
	
	currentelement = 0;
	
	while (currentelement  < element){
		argstring.erase(0, argstring.find(" ") + 1);
		currentelement += 1;
	}
	
	returnstring = argstring.substr(0, argstring.find(" "));
}

void listFilesInDirectory (std::vector<std::string>& files, std::string directory) {
	DIR 								*dir;
	struct dirent 						*entry;
	
	dir = opendir(directory.c_str());
	entry = readdir(dir);
	
	while (entry != NULL){
		if (entry->d_type != DT_DIR){
			files.push_back(entry->d_name);
		}
		entry = readdir(dir);
	}
	
	closedir(dir);
	
}

void listDirectory (std::vector<std::string>& files, std::vector<std::string>& directories, std::string directory) {
	DIR 								*dir;
	struct dirent 						*entry;
	
	dir = opendir(directory.c_str());
	entry = readdir(dir);
	
	while (entry != NULL){
		if (entry->d_type != DT_DIR){
			files.push_back(entry->d_name);
		}else{
			directories.push_back(entry->d_name);
		}
		entry = readdir(dir);
	}
	
	closedir(dir);
	
}

bool fileExists(std::string filename) {
	
	struct stat 		buf;

	if(stat(filename.c_str(), &buf) == 0) {
		return(true);
	} else {
		return(false);
	}
	
}

std::string zeroPadNumber(int num) {
	std::stringstream ss;
	std::string ret;
	
	ss << num; 
	ss >> ret;
	int str_length = ret.length();
	for (int i = 0; i < 3 - str_length; i++)
		ret = "0" + ret;
	return ret;
}

void getDatabaseLocation (std::string& databasepath) {
	
	struct passwd* 		pw;
	uid_t 				uid;
	
	uid = geteuid();
	pw = getpwuid(uid);
	
	databasepath = std::string(pw->pw_dir) + "/.simple.sqlite";
}

void SQLQuery(std::string& databasepath, std::vector< std::map<std::string, std::string> >& result, std::string sql){
	sqlite3*											db;
	sqlite3_stmt*										stmt;
	int 												rc;
	int													id;
	std::map<std::string, std::string>					resultlinekeyval;
	
	result.clear();
	
	if(sqlite3_open(databasepath.c_str(), &db)){
		std::cout << "Failed to open database " << sqlite3_errmsg(db) << std::endl;
	}
	sqlite3_busy_timeout(db, 5000);
	if(sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, NULL) != SQLITE_OK ){
		std::cout << "Failed to query database " << sqlite3_errmsg(db) << std::endl;
	}

	while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
		resultlinekeyval.clear();
		int id = sqlite3_column_int(stmt, 0);
		resultlinekeyval["id"] = std::to_string(id);
		for(int i = 1; i < sqlite3_column_count(stmt); i++){
			const char *columnname = sqlite3_column_name(stmt, i);
			const unsigned char *columnvalue = sqlite3_column_text(stmt, i);
			resultlinekeyval[std::string(columnname)] = std::string((char*)columnvalue);
		}
		result.push_back(resultlinekeyval);
	}

	if (rc != SQLITE_DONE) {
		std::cout << "Failed to parse query " << sqlite3_errmsg(db) << std::endl;
	}

	sqlite3_finalize(stmt);
	sqlite3_close(db);
}

void SQLQuery(std::string& databasepath, std::vector< std::map<std::string, std::string> >& result, std::string sql, int& rowid){
	sqlite3*											db;
	sqlite3_stmt*										stmt;
	int 												rc;
	int													id;
	std::map<std::string, std::string>					resultlinekeyval;
	
	std::cout << sql << std::endl;
	
	if(sqlite3_open(databasepath.c_str(), &db)){
		std::cout << "Failed to open database " << sqlite3_errmsg(db) << std::endl;
	}
	sqlite3_busy_timeout(db, 5000);
	if(sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, NULL) != SQLITE_OK ){
		std::cout << "Failed to query database " << sqlite3_errmsg(db) << std::endl;
	}

	while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
		resultlinekeyval.clear();
		int id = sqlite3_column_int(stmt, 0);
		resultlinekeyval["id"] = std::to_string(id);
		for(int i = 1; i < sqlite3_column_count(stmt); i++){
			const char *columnname = sqlite3_column_name(stmt, i);
			const unsigned char *columnvalue = sqlite3_column_text(stmt, i);
			resultlinekeyval[std::string(columnname)] = std::string((char*)columnvalue);
		}
		result.push_back(resultlinekeyval);
	}

	if (rc != SQLITE_DONE) {
		std::cout << "Failed to parse query " << sqlite3_errmsg(db) << std::endl;
	}

	if(sqlite3_prepare_v2(db, "SELECT last_insert_rowid();", -1, &stmt, NULL)  != SQLITE_OK ){
		std::cout << "Failed to query database " << sqlite3_errmsg(db) << std::endl;
	}	
	
	while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
		rowid = sqlite3_column_int(stmt, 0);
	}
	
	if (rc != SQLITE_DONE) {
		std::cout << "Failed to parse query " << sqlite3_errmsg(db) << std::endl;
	}

	sqlite3_finalize(stmt);
	sqlite3_close(db);
}

void SQLQuery(std::string& databasepath, std::string sql, int& rowid){
	sqlite3*											db;
	sqlite3_stmt*										stmt;
	int 												rc;
	int													id;
	std::map<std::string, std::string>					resultlinekeyval;
	
	std::cout << sql << std::endl;
	
	if(sqlite3_open(databasepath.c_str(), &db)){
		std::cout << "Failed to open database " << sqlite3_errmsg(db) << std::endl;
	}
	sqlite3_busy_timeout(db, 5000);
	if(sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, NULL) != SQLITE_OK ){
		std::cout << "Failed to query database " << sqlite3_errmsg(db) << std::endl;
	}

	while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
	}

	if (rc != SQLITE_DONE) {
		std::cout << "Failed to parse query " << sqlite3_errmsg(db) << std::endl;
	}

	if(sqlite3_prepare_v2(db, "SELECT last_insert_rowid();", -1, &stmt, NULL)  != SQLITE_OK ){
		std::cout << "Failed to query database " << sqlite3_errmsg(db) << std::endl;
	}	
	
	while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
		rowid = sqlite3_column_int(stmt, 0);
	}
	
	if (rc != SQLITE_DONE) {
		std::cout << "Failed to parse query " << sqlite3_errmsg(db) << std::endl;
	}

	sqlite3_finalize(stmt);
	sqlite3_close(db);
}

void SQLQuery(std::string& databasepath, std::string sql) {
	sqlite3*											db;
	sqlite3_stmt*										stmt;
	int 												rc;
	int													id;
	std::map<std::string, std::string>					resultlinekeyval;
	
	std::cout << sql << std::endl;
	
	if(sqlite3_open(databasepath.c_str(), &db)){
		std::cout << "Failed to open database " << sqlite3_errmsg(db) << std::endl;
	}
	sqlite3_busy_timeout(db, 5000);
	if(sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, NULL) != SQLITE_OK ){
		std::cout << "Failed to query database " << sqlite3_errmsg(db) << std::endl;
	}
	
	while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
	}

	if (rc != SQLITE_DONE) {
		std::cout << "Failed to parse query " << sqlite3_errmsg(db) << std::endl;
	}
	
	sqlite3_finalize(stmt);
	sqlite3_close(db);
}

void readUniDoc(UniDoc* unidoc, std::string filename){
	
	std::ifstream			input;
	std::string 			line;
	
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

void writeUniDoc(UniDoc* unidoc, std::string filename){
	
	std::ofstream							output;
	std::vector<std::string>::iterator		datait;
	
	output.open(filename.c_str());
	
	if(output.is_open()){
		for(datait = unidoc->data.begin(); datait != unidoc->data.end(); datait++){
			output << *datait << "\n";
		}
		output.close();
	} else {
		std::cout << "Failed to open unidoc" << std::endl;
	}
}

void readMRCFrame(MRCFile* mrcfile, std::string filename, int frameid) {
	
	std::ifstream 				input;
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
	
	input.close();
	
}

void extractMRCBox(MRCFile* mrcfile, MRCFile* box, int xcoord, int ycoord, int boxsize) {
	
	int			xit;
	int			yit;
	float 		count = 0;
	
	box->data = new float[boxsize * boxsize];
	for(yit = 0; yit < boxsize; yit++){
		for(xit = 0; xit < boxsize; xit++){
			box->data[(yit * boxsize) + xit] = mrcfile->data[((ycoord + yit) * mrcfile->nx) + (xcoord + xit)];
			count += mrcfile->data[((ycoord + yit) * mrcfile->nx) + (xcoord + xit)];
		}
	}
	box->nx = boxsize;
	box->ny = boxsize;
	box->nz = boxsize;
}

bool getUniDocValue(UniDoc* unidoc, int line, std::string key, std::string& value){
	
	std::size_t 			startpos;
	std::size_t 			endpos;
	
	startpos = unidoc->data[line].find(key + "=");
	
	value.clear();

	if(startpos != std::string::npos) {
		endpos = unidoc->data[line].find(" ", startpos + 1);
		value = unidoc->data[line].substr(startpos + key.length() + 1, endpos - startpos - key.length() - 1);
		return(true);
	}else{
		return(false);
	}
	
}

bool deleteUniDocValue(UniDoc* unidoc, int line, std::string key){
	
	std::size_t 			startpos;
	std::size_t 			endpos;
	
	startpos = unidoc->data[line].find(key + "=");
	
	if(startpos != std::string::npos) {
		endpos = unidoc->data[line].find(" ", startpos + 1);
		unidoc->data[line].erase(startpos, endpos);
		return(true);
	}else{
		return(false);
	}
	
}

void addUniDocKeyVal(UniDoc* unidoc, int line, std::string key, std::string value){
	
	std::size_t 			startpos;
	std::size_t 			endpos;
	
	startpos = unidoc->data[line].find(key + "=");
	if(startpos != std::string::npos) {
		endpos = unidoc->data[line].find(" ", startpos + 1);
		unidoc->data[line].erase(startpos, endpos  - startpos);
	}
	unidoc->data[line] += " " + key + "=" + value;

}

void addUnidocParticleIdentifiers(UniDoc* unidoc, std::string& unidocfilename){
	
	std::string					stacktabfilename;
	std::string					rootdirectory;
	char*						dirc;
	UniDoc*						stacktab;
	int							stacktabit;
	std::string					stackframestart;
	std::string					stackframeend;
	std::string					stackfile;
	int							stackframeit;
	int							frameit;
	
	dirc = strdup(unidocfilename.c_str());
	rootdirectory = std::string(dirname(dirc));
	stacktabfilename = rootdirectory + "/stktab_info.txt";

	if(fileExists(stacktabfilename)){
		stacktab = new UniDoc();
		readUniDoc(stacktab, stacktabfilename);
		stacktab->data.erase(stacktab->data.begin());
		for(stacktabit = 0; stacktabit < stacktab->data.size(); stacktabit++){
			getElementFromString(stacktab->data[stacktabit], stackfile, 0);
			getElementFromString(stacktab->data[stacktabit], stackframestart, 2);
			getElementFromString(stacktab->data[stacktabit], stackframeend, 3);
			frameit = 0;
			for(stackframeit = std::stoi(stackframestart) - 1; stackframeit < std::stoi(stackframeend); stackframeit++){
				addUniDocKeyVal(unidoc, stackframeit, "stackfile", stackfile);
				addUniDocKeyVal(unidoc, stackframeit, "frameid", std::to_string(frameit));
				frameit++;
			}
		}
		delete stacktab;
	}
}

void writeRelionStar(UniDoc* unidoc, std::string filename){

	UniDoc * 		relionstar;
	int				selectionit;
	std::string		starline;
	std::string		value;
	std::string		stackfile;
	std::string		state;
	float			statef;
	float			defocus;
	
	relionstar = new UniDoc();
	relionstar->data.push_back("");
	relionstar->data.push_back("data_");
	relionstar->data.push_back("");
	relionstar->data.push_back("loop_");
	relionstar->data.push_back("_rlnImageName");
	relionstar->data.push_back("_rlnVoltage");
	relionstar->data.push_back("_rlnDefocusU");
	relionstar->data.push_back("_rlnDefocusV");
	relionstar->data.push_back("_rlnDefocusAngle");
	relionstar->data.push_back("_rlnSphericalAberration");
	relionstar->data.push_back("_rlnAmplitudeContrast");
	relionstar->data.push_back("_rlnMagnification");
	relionstar->data.push_back("_rlnDetectorPixelSize");
	
	for(selectionit = 0; selectionit < unidoc->data.size(); selectionit++){
		getUniDocValue(unidoc, selectionit, "state", state);
		statef = std::stof(state);
		if (statef != 0){
			getUniDocValue(unidoc, selectionit, "frameid", value);
			starline = value + "@";
			getUniDocValue(unidoc, selectionit, "stackfile", value);
			stackfile = value + "s";
			symlink(value.c_str(), stackfile.c_str());
			starline += value + "s";
			getUniDocValue(unidoc, selectionit, "kv", value);
			starline += " " + value;
			getUniDocValue(unidoc, selectionit, "dfx", value);
			defocus = std::stof(value) * 10000;
			starline += " " + std::to_string(defocus);
			getUniDocValue(unidoc, selectionit, "dfy", value);
			defocus = std::stof(value) * 10000;
			starline += " " + std::to_string(defocus);
			getUniDocValue(unidoc, selectionit, "angast", value);
			starline += " " + value;
			getUniDocValue(unidoc, selectionit, "cs", value);
			starline += " " + value;
			getUniDocValue(unidoc, selectionit, "fraca", value);
			starline += " " + value;
			starline += " 10000";
			getUniDocValue(unidoc, selectionit, "smpd", value);
			starline += " " + value;
			relionstar->data.push_back(starline);
		}
	}
	filename = filename.substr(0, filename.length() -7);
	filename += ".star";
	std::cout << "WRITING STAR FILE " + filename << std::endl;
	writeUniDoc(relionstar, filename);
	delete relionstar;
}

void addJob (std::string& jobname, std::string jobdescription, std::string& jobfolder, std::string& jobtype, std::string& table, std::string& runstring, int &jobid) {
	
	std::string				databasepath;
	std::string::size_type  i;
	
	getDatabaseLocation(databasepath);

	SQLQuery(databasepath, "INSERT INTO " + std::string(table) + "(jobname, jobdescription, jobfolder, jobtype, jobstatus, jobpid, jobrerun) VALUES('" + jobname + "','" + jobdescription + "','"	+ jobfolder + "','"	+ jobtype + "','" + "External" + "','"	+ std::to_string(0) + "','" + runstring + "');", jobid);
}

void deleteJob(JSONResponse* response, struct http_message* message) {
	
	std::string				databasepath;
	std::string				table;
	std::string				jobid;
	
	getDatabaseLocation(databasepath);
	
	if (getRequestVariable(message, "table", table) && getRequestVariable(message, "jobid", jobid)) {
		SQLQuery(databasepath, "DELETE FROM " + table + " WHERE jobid=" + jobid + ";");
	}
}

void updateJobStatus (std::string jobstatus, std::string& table, int &jobid) {
	
	std::string				databasepath;
	
	getDatabaseLocation(databasepath);

	SQLQuery(databasepath, "UPDATE " + table + " SET  jobstatus = \'" + jobstatus + "\' WHERE jobid=" + std::to_string(jobid) +";");
}

void updateJobFolder (std::string folder, std::string& table, int &jobid) {
	
	std::string 	databasepath;
	
	getDatabaseLocation(databasepath);
	
	SQLQuery(databasepath, "UPDATE " + table + " SET  jobfolder = \'" + folder + "\' WHERE jobid=" + std::to_string(jobid) +";");
}

void updateJobPID (pid_t pid, std::string& table, int &jobid) {
	
	std::string 	databasepath;
	
	getDatabaseLocation(databasepath);
	
	SQLQuery(databasepath, "UPDATE " + table + " SET  jobpid = \'" + std::to_string(pid) + "\' WHERE jobid=" + std::to_string(jobid) +";");
}

void updateJobRerun (std::string runstring, std::string& table, int &jobid) {
	
	std::string 	databasepath;
	
	getDatabaseLocation(databasepath);
	
	SQLQuery(databasepath, "UPDATE " + table + " SET  jobrerun = \'" + runstring + "\' WHERE jobid=" + std::to_string(jobid) +";");
}

void ctffindPre (std::string micrographsdirectory, std::string unbluroutput){
	
	UniDoc* 									outputunidoc;
	UniDoc* 									unidoc;
	std::ofstream								filetab;								
	std::vector<std::string> 					files;
	std::vector<std::string>::iterator			filesit;
	int											it;
	int											outit;
	std::string									value;
	std::string									micrographpath;
	std::string									micrographlink;
	char*										dname;
	char*										dirc;
	int											status;
	
	outputunidoc = new UniDoc();
	
	filetab.open("micrographs.txt");
	
	if(filetab.is_open()){
	//	status = mkdir("micrographs", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		
		if(micrographsdirectory.size() > 0 && fileExists(micrographsdirectory)){
			status = symlink(micrographsdirectory.c_str(), "micrographs");
			listFilesInDirectory(files, "micrographs" );
			outit = 0;
			for(it = 0; it < files.size(); it++) {
				if(files[it].find(".mrc") != std::string::npos) {
				//	micrographpath = micrographsdirectory + "/" + files[it];
				//	micrographlink = "micrographs/" + files[it];
				//	status = symlink(micrographpath.c_str(), micrographlink.c_str());
				//	if(status == 0){
						outputunidoc->data.push_back("");
						addUniDocKeyVal(outputunidoc, outit, "intg", "micrographs/" + files[it]);
						filetab << micrographsdirectory + "/" + files[it] + "\n";
				//	}
					outit++;
				}
			}
		} else if(unbluroutput.size() > 0 && fileExists(unbluroutput)){
			unidoc = new UniDoc();
			readUniDoc(unidoc, unbluroutput);
			dirc = strdup(unbluroutput.c_str());
			dname = dirname(dirc);
			for(it = 0; it < unidoc->data.size(); it++){
				outputunidoc->data.push_back("");
			//	getUniDocValue(unidoc, it, "movie", std::string(dname) + "/" + value);
			//	addUniDocKeyVal(outputunidoc, it, std::string(dname) + "/" + "movie", value);
				getUniDocValue(unidoc, it, "intg", value);
				addUniDocKeyVal(outputunidoc, it, "intg", std::string(dname) + "/" + value);
				getUniDocValue(unidoc, it, "forctf", value);
				addUniDocKeyVal(outputunidoc, it, "forctf", std::string(dname) + "/" + value);
				filetab << std::string(dname) + "/" + value + "\n";
				getUniDocValue(unidoc, it, "pspec", value);
				addUniDocKeyVal(outputunidoc, it, "pspec", std::string(dname) + "/" + value);
				getUniDocValue(unidoc, it, "thumb", value);
				addUniDocKeyVal(outputunidoc, it, "thumb", std::string(dname) + "/" + value);
			}
			delete unidoc;
		}
	}

	writeUniDoc(outputunidoc, "ctffind_in.simple");
	filetab.close();
	delete outputunidoc;
}

void ctffindPost (std::string directory){
	
	int											it;
	bool										newmicrograph;
	UniDoc* 									ctffindunidocin;
	UniDoc*										ctfparams;
	std::string									intgfilename;
	std::string									ctffindfit;
	std::string									value;
	
	ctffindunidocin = new UniDoc();
	ctfparams = new UniDoc();
	
	if(fileExists(directory + "/ctffind_in.simple") && fileExists(directory + "/ctffind_output_merged.txt")){
		readUniDoc(ctffindunidocin, directory + "/ctffind_in.simple");
		readUniDoc(ctfparams, directory + "/ctffind_output_merged.txt");
		
		for(it = 0; it < ctffindunidocin->data.size(); it++){
			getUniDocValue(ctfparams, it, "kv", value);
			addUniDocKeyVal(ctffindunidocin, it, "kv", value);
			getUniDocValue(ctfparams, it, "cs", value);
			addUniDocKeyVal(ctffindunidocin, it, "cs", value);
			getUniDocValue(ctfparams, it, "fraca", value);
			addUniDocKeyVal(ctffindunidocin, it, "fraca", value);
			getUniDocValue(ctfparams, it, "dfx", value);
			addUniDocKeyVal(ctffindunidocin, it, "dfx", value);
			getUniDocValue(ctfparams, it, "dfy", value);
			addUniDocKeyVal(ctffindunidocin, it, "dfy", value);
			getUniDocValue(ctfparams, it, "angast", value);
			addUniDocKeyVal(ctffindunidocin, it, "angast", value);
			getUniDocValue(ctfparams, it, "ctfres", value);
			addUniDocKeyVal(ctffindunidocin, it, "ctfres", value);
			if(getUniDocValue(ctffindunidocin, it, "forctf", ctffindfit)){
				ctffindfit.replace(ctffindfit.end() - 4, ctffindfit.end(), "_ctffind_diag.mrc");
			} else if (getUniDocValue(ctffindunidocin, it, "intg", ctffindfit)){
				ctffindfit.replace(ctffindfit.end() - 4, ctffindfit.end(), "_ctffind_diag.mrc");
			}
			if(fileExists(ctffindfit)){
				addUniDocKeyVal(ctffindunidocin, it, "ctffindfit", ctffindfit);
			} 
		} 
	} 
	
	delete ctfparams;
	
	writeUniDoc(ctffindunidocin, directory + "/ctffind_out.simple");
	
	delete ctffindunidocin;
}

void ctffitPre (std::string micrographsdirectory, std::string unbluroutput){
	
	UniDoc* 									outputunidoc;
	UniDoc* 									unidoc;
	std::ofstream								filetab;								
	std::vector<std::string> 					files;
	std::vector<std::string>::iterator			filesit;
	int											it;
	int											outit;
	std::string									value;
	std::string									micrographpath;
	std::string									micrographlink;
	char*										dname;
	char*										dirc;
	int											status;
	
	outputunidoc = new UniDoc();
	
	
	filetab.open("micrographs.txt");
	
	if(filetab.is_open()){
		if(micrographsdirectory.size() > 0 && fileExists(micrographsdirectory)){
			//status = symlink(micrographsdirectory.c_str(), "micrographs");
			status = mkdir("micrographs", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			listFilesInDirectory(files, micrographsdirectory);
			outit = 0;
			for(it = 0; it < files.size(); it++) {
				if(files[it].find(".mrc") != std::string::npos) {
					micrographpath = micrographsdirectory + "/" + files[it];
					micrographlink = "micrographs/" + files[it];
					status = symlink(micrographpath.c_str(), micrographlink.c_str());
					outputunidoc->data.push_back("");
					addUniDocKeyVal(outputunidoc, outit, "intg", "micrographs/" + files[it]);
					//filetab << micrographsdirectory + "/" + files[it] + "\n";
					filetab << "micrographs/" + files[it] + "\n";
					outit++;
				}
			}
		} else if(unbluroutput.size() > 0 && fileExists(unbluroutput)){
			unidoc = new UniDoc();
			readUniDoc(unidoc, unbluroutput);
			dirc = strdup(unbluroutput.c_str());
			dname = dirname(dirc);
			status = mkdir("micrographs", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			for(it = 0; it < unidoc->data.size(); it++){
				outputunidoc->data.push_back("");
				
				getUniDocValue(unidoc, it, "intg", value);
				addUniDocKeyVal(outputunidoc, it, "intg", "micrographs/" + value);
				micrographpath = std::string(dname) + "/" + value;
				micrographlink = "micrographs/" + value;
				status = symlink(micrographpath.c_str(), micrographlink.c_str());
				
				getUniDocValue(unidoc, it, "forctf", value);
				addUniDocKeyVal(outputunidoc, it, "forctf", "micrographs/" + value);
				micrographpath = std::string(dname) + "/" + value;
				micrographlink = "micrographs/" + value;
				status = symlink(micrographpath.c_str(), micrographlink.c_str());
				filetab << "micrographs/" + value + "\n";
				
				getUniDocValue(unidoc, it, "pspec", value);
				addUniDocKeyVal(outputunidoc, it, "pspec", "micrographs/" + value);
				micrographpath = std::string(dname) + "/" + value;
				micrographlink = "micrographs/" + value;
				status = symlink(micrographpath.c_str(), micrographlink.c_str());
				
				getUniDocValue(unidoc, it, "thumb", value);
				addUniDocKeyVal(outputunidoc, it, "thumb", "micrographs/" + value);
				micrographpath = std::string(dname) + "/" + value;
				micrographlink = "micrographs/" + value;
				status = symlink(micrographpath.c_str(), micrographlink.c_str());
			}
			delete unidoc;
		}
	}

	writeUniDoc(outputunidoc, "ctffit_in.simple");
	filetab.close();
	delete outputunidoc;
}

void ctffitPost (std::string directory){
	
	int											it;
	bool										newmicrograph;
	UniDoc* 									ctffindunidocin;
	UniDoc*										ctfparams;
	std::string									intgfilename;
	std::string									ctffindfit;
	std::string									value;
	
	ctffindunidocin = new UniDoc();
	ctfparams = new UniDoc();
	
	if(fileExists(directory + "/ctffit_in.simple") && fileExists(directory + "/ctffit_output_merged.txt")){
		readUniDoc(ctffindunidocin, directory + "/ctffit_in.simple");
		readUniDoc(ctfparams, directory + "/ctffit_output_merged.txt");
		
		for(it = 0; it < ctffindunidocin->data.size(); it++){
			getUniDocValue(ctfparams, it, "kv", value);
			addUniDocKeyVal(ctffindunidocin, it, "kv", value);
			getUniDocValue(ctfparams, it, "cs", value);
			addUniDocKeyVal(ctffindunidocin, it, "cs", value);
			getUniDocValue(ctfparams, it, "fraca", value);
			addUniDocKeyVal(ctffindunidocin, it, "fraca", value);
			getUniDocValue(ctfparams, it, "dfx", value);
			addUniDocKeyVal(ctffindunidocin, it, "dfx", value);
			getUniDocValue(ctfparams, it, "dfy", value);
			addUniDocKeyVal(ctffindunidocin, it, "dfy", value);
			getUniDocValue(ctfparams, it, "angast", value);
			addUniDocKeyVal(ctffindunidocin, it, "angast", value);
			getUniDocValue(ctfparams, it, "ctfres", value);
			addUniDocKeyVal(ctffindunidocin, it, "ctfres", value);
			if(getUniDocValue(ctffindunidocin, it, "forctf", ctffindfit)){
				ctffindfit.replace(ctffindfit.end() - 4, ctffindfit.end(), "_ctffit_diag.mrc");
			} else if (getUniDocValue(ctffindunidocin, it, "intg", ctffindfit)){
				ctffindfit.replace(ctffindfit.end() - 4, ctffindfit.end(), "_ctffit_diag.mrc");
			}
			if(fileExists(ctffindfit)){
				addUniDocKeyVal(ctffindunidocin, it, "ctffindfit", ctffindfit);
			} 
		} 
	} 
	
	delete ctfparams;
	
	writeUniDoc(ctffindunidocin, directory + "/ctffit_out.simple");
	
	delete ctffindunidocin;
}

void ini3DPre (std::string simpleinput, std::string& command){ // NEEDS WORK
	
	UniDoc*										inputunidoc;
	std::ofstream								stktab;
	int											inputit;
	std::string									stackfile;
	std::string									stackfilepre;
	std::string									inputroot;
	char*										dname;
	char*										dirc;
	std::string									stackfilepath;
	std::string									stackfilelink;
	std::string									stackfilerealpath;
	
	std::cout << "inipre" << std::endl;
	if(fileExists(simpleinput)){
		mkdir("particles", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		stktab.open("stkparts.txt");
		dirc = strdup(simpleinput.c_str());
		inputroot = std::string(dirname(dirc));
		if(fileExists(inputroot + "/prime2D_selected_clusters.mrc")){
			command += " stk=" + inputroot + "/prime2D_selected_clusters.mrc";
		}
		std::cout << inputroot << std::endl;
		if(stktab.is_open()){
			inputunidoc = new UniDoc();
			readUniDoc(inputunidoc, simpleinput);
			writeUniDoc(inputunidoc, "ini3d_in.simple");
			stackfilepre = "";
			for(inputit = 0; inputit < inputunidoc->data.size(); inputit++) {
				getUniDocValue(inputunidoc, inputit, "stackfile", stackfile);
				if(stackfile != stackfilepre){
					stackfilepath = inputroot + "/" + stackfile;
					stackfilerealpath = std::string(realpath(stackfilepath.c_str(), NULL));
					dirc = strdup(stackfilepath.c_str());
					stktab << "particles/" + std::string(basename(dirc)) << "\n";
					stackfilelink = "particles/" + std::string(basename(dirc));
					symlink(stackfilerealpath.c_str(), stackfilelink.c_str());
					stackfilepre = stackfile;
				}
				deleteUniDocValue(inputunidoc, inputit, "stackfile");
				deleteUniDocValue(inputunidoc, inputit, "frameid");
			}
			stktab.close();
			writeUniDoc(inputunidoc, "ini3d_in.txt");
			delete inputunidoc;
		}
	}
}

void ini3DPost (std::string directory){
	
	FILE*						stream;
	std::string					command;
	std::string					prime3ddoc;
	std::string					stackfile;
	std::string					frameid;
	std::vector<std::string>	files;
	int							it;
	int							iterationcount;
	int							status;
	UniDoc*						inputunidoc;
	UniDoc*						outputunidoc;
	int							inputit;
	
	status = chdir(directory.c_str());
	
	if(fileExists("ini3d_in.simple") && fileExists("ini3d_in.txt")){
		listFilesInDirectory(files, directory);
		iterationcount = 0;
		for(it = 0; it < files.size(); it++) {
			if(files[it].find("prime3Ddoc") != std::string::npos) {
				iterationcount++;
			}
		}	
		prime3ddoc = "prime3Ddoc_" + zeroPadNumber(iterationcount) + ".txt";
		command = "simple_exec prg=map2ptcls_doc oritab=ini3d_in.txt oritab3D="+ prime3ddoc +" outfile=ini3d_particles_out.txt >> simple_job.log";
		std::cout << command << std::endl;
		stream = popen(command.c_str(), "r");
		pclose(stream);
		inputunidoc = new UniDoc();
		outputunidoc = new UniDoc();
		readUniDoc(inputunidoc, "ini3d_in.simple");
		readUniDoc(outputunidoc, "ini3d_particles_out.txt");
		for(inputit = 0; inputit < inputunidoc->data.size(); inputit++) {
			getUniDocValue(inputunidoc, inputit, "stackfile", stackfile);
			getUniDocValue(inputunidoc, inputit, "frameid", frameid);
			addUniDocKeyVal(outputunidoc, inputit, "stackfile", stackfile);
			addUniDocKeyVal(outputunidoc, inputit, "frameid", frameid);
		}
		
		delete inputunidoc;
		writeUniDoc(outputunidoc, "ini3d_particles_out.simple");
		delete outputunidoc;
	}
}

void extractPre (std::string boxdirectory, std::string micrographsdirectory, std::string unidocname){
	
	UniDoc* 									outputunidoc;
	UniDoc* 									unidoc;	
	UniDoc* 									unidocout;	
	UniDoc* 									boxtab;							
	char*										bname;
	char*										dname;
	char*										dirc;
	char*										dird;
	int											it;
	std::string									intg;
	std::string									boxfile;
	
	unidoc = new UniDoc();
	unidocout = new UniDoc();
	boxtab = new UniDoc();

	if(micrographsdirectory.size() > 0){
		
	} else if(unidocname.size() > 0){
		readUniDoc(unidoc, unidocname);
		for(it = 0; it < unidoc->data.size(); it++){
			getUniDocValue(unidoc, it, "intg", intg);
			dirc = strdup(intg.c_str());
			bname = basename(dirc);
			dird = strdup(unidocname.c_str());
			dname = dirname(dird);
		
			dirc = strdup(intg.c_str());
			bname = basename(dirc);
			boxfile = boxdirectory + "/" + std::string(bname);
			boxfile.replace(boxfile.end() - 3, boxfile.end(), "box");
			
			if(fileExists(boxfile)){
				deleteUniDocValue(unidoc, it, "intg");
				addUniDocKeyVal(unidoc, it, "intg", std::string(dname) + "/" + intg);
				unidocout->data.push_back(unidoc->data[it]);
				boxtab->data.push_back(boxfile);
			}
		}
		writeUniDoc(unidocout, "extract_in.simple");
		writeUniDoc(unidocout, "unidoc_in.txt");
		writeUniDoc(boxtab, "boxtab.txt");
	}

	delete unidoc;
	delete unidocout;
	delete boxtab;
}

void extractPost (std::string directory){
	
	UniDoc*							unidoc;
	UniDoc*							parts;
	UniDoc*							unidocpart;
	std::vector<std::string>		files;
	int								filesit;
	std::string						stack;
	int								stackit;
	
	listFilesInDirectory(files, directory);
	unidoc = new UniDoc();
	parts = new UniDoc();
	
	for(filesit = 0; filesit < files.size(); filesit++){
		if(files[filesit].find("extract_params") != std::string::npos){
			stack = files[filesit];
			stack.replace(stack.begin(), stack.begin() + 14, "ptcls_from");
			stack.replace(stack.end() - 3, stack.end(), "mrc");
			if(fileExists(stack)){
				
				unidocpart = new UniDoc();
				readUniDoc(unidocpart, files[filesit]);
				for(stackit = 0; stackit < unidocpart->data.size(); stackit++){
					addUniDocKeyVal(unidocpart, stackit, "stackfile", stack);
					addUniDocKeyVal(unidocpart, stackit, "frameid", std::to_string(stackit + 1));
					unidoc->data.push_back(unidocpart->data[stackit]);
				}
				parts->data.push_back("stackfile=" + stack + " particlecount=" + std::to_string(unidocpart->data.size()));
				delete unidocpart;
			}
		}
	}
	
	writeUniDoc(unidoc, "extract_out.simple");
	writeUniDoc(parts, "extract_parts.txt");
	delete unidoc;
	delete parts;
}

void pickPre (std::string simpleinput, std::string micrographsdirectory, std::string pickrefs, std::string pickvol, std::string pgrp, std::string pcontrast) {
	
	std::string		command;
	std::string		micrograph;
	std::string		rootdir;
	FILE*			stream;
	UniDoc*			simplein;
	UniDoc*			filetab;
	int 			i;
	char*			basec;

	if(pickrefs != "" && pcontrast != ""){
		command = "simple_exec prg=makepickrefs pgrp=" + pgrp + " stk=" + pickrefs + " nthr=1 pcontrast=white >> simple_job.log";
		std::cout << command << std::endl;
		stream = popen(command.c_str(), "r");
		pclose(stream);
	} else if (pickrefs != "" && pcontrast == ""){
		command = "simple_exec prg=makepickrefs pgrp=" + pgrp + " stk=" + pickrefs + " nthr=1 >> simple_job.log";
		std::cout << command << std::endl;
		stream = popen(command.c_str(), "r");
		pclose(stream);
	} else if (pickvol != "" && pcontrast != ""){		
		command = "simple_exec prg=makepickrefs pgrp=" + pgrp + " vol1=" + pickrefs + " nthr=1 pcontrast=white >> simple_job.log";
		std::cout << command << std::endl;
		stream = popen(command.c_str(), "r");
		pclose(stream);
	} else if (pickvol != "" && pcontrast == ""){
		command = "simple_exec prg=makepickrefs pgrp=" + pgrp + " vol1=" + pickrefs + " nthr=1 >> simple_job.log";
		std::cout << command << std::endl;
		stream = popen(command.c_str(), "r");
		pclose(stream);
	}
	
	if (simpleinput != "" && fileExists	(simpleinput)){
		
		simplein = new UniDoc();
		filetab = new UniDoc();
		
		readUniDoc(simplein, simpleinput);
		
		basec = strdup(simpleinput.c_str());
		rootdir = std::string(dirname(basec));
		
		for(i = 0; i < simplein->data.size(); i++){
			getUniDocValue(simplein, i, "intg", micrograph);
			filetab->data.push_back(rootdir + "/" + micrograph);
		}
		
		writeUniDoc(filetab, "picktab.txt");
		writeUniDoc(simplein, "autopick_in.simple");
		
		delete simplein;
		delete filetab;
		
	}
}

void preprocPost (std::string directory){
	
	std::vector<std::string> 					files;
	std::vector<std::string>::iterator			filesit;
	int											datait;
	bool										newmicrograph;
	UniDoc* 									micrographunidoc;
	UniDoc*										preprocunidoc;
	std::string									unidocname;
	std::string									ctffindfit;
	
	listFilesInDirectory(files, directory + "/pipeline/unidocs");
	
	preprocunidoc = new UniDoc();
	
	if(fileExists(directory + "/preproc_out.simple")){
		readUniDoc(preprocunidoc, directory + "/preproc_out.simple");
	}
	
	for(filesit = files.begin(); filesit != files.end(); filesit++) {
		if(std::string(*filesit).find("unidoc") != std::string::npos) {
			newmicrograph = true;
			for(datait = 0; datait < preprocunidoc->data.size(); datait++) {
				if(getUniDocValue(preprocunidoc, datait, "unidoc", unidocname)) {
					if(unidocname == std::string(*filesit)){
						newmicrograph = false;
					}
				}
			}
			if(newmicrograph){
				micrographunidoc = new UniDoc();
				readUniDoc(micrographunidoc, directory + "/pipeline/unidocs/" + std::string(*filesit));
				ctffindfit = std::string(*filesit);
				ctffindfit.replace(ctffindfit.end() - 4, ctffindfit.end(), "_forctf_ctffind_diag.mrc");
				ctffindfit.replace(ctffindfit.begin(), ctffindfit.begin() + 14, "pipeline/ctf/");
				addUniDocKeyVal(micrographunidoc, 0, "ctffindfit", ctffindfit);
				addUniDocKeyVal(micrographunidoc, 0, "unidoc", std::string(*filesit));
				preprocunidoc->data.push_back(micrographunidoc->data[0]);
				delete micrographunidoc;
			}
		}
	} 
	writeUniDoc(preprocunidoc, directory + "/preproc_out.simple");
	delete preprocunidoc;
}

void prime2DPre (std::string simpleinput){
	
	UniDoc*										inputunidoc;
	std::ofstream								stktab;
	int											inputit;
	std::string									stackfile;
	std::string									stackfilepre;
	std::string									inputroot;
	char*										dname;
	char*										dirc;
	std::string									stackfilepath;
	std::string									stackfilelink;
	std::string									stackfilerealpath;
	
	mkdir("particles", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	stktab.open("stkparts.txt");
	dirc = strdup(simpleinput.c_str());
	inputroot = std::string(dirname(dirc));
	
	if(fileExists(simpleinput) && stktab.is_open()){
		inputunidoc = new UniDoc();
		readUniDoc(inputunidoc, simpleinput);
		//if(simpleinput.find("preproc") == std::string::npos) {
		//	std::cout << "adding params" << std::endl;
		//	addUnidocParticleIdentifiers(inputunidoc, simpleinput);
		//}
		stackfilepre = "";
		for(inputit = 0; inputit < inputunidoc->data.size(); inputit++) {
			getUniDocValue(inputunidoc, inputit, "stackfile", stackfile);
			if(stackfile != stackfilepre){
				stackfilepath = inputroot + "/" + stackfile;
				stackfilerealpath = std::string(realpath(stackfilepath.c_str(), NULL));
				dirc = strdup(stackfilepath.c_str());
				stktab << "particles/" + std::string(basename(dirc)) << "\n";
				std::cout << stackfilepath << " " << stackfilerealpath << std::endl;
				stackfilelink = "particles/" + std::string(basename(dirc));
				symlink(stackfilerealpath.c_str(), stackfilelink.c_str());
				stackfilepre = stackfile;
			}
			deleteUniDocValue(inputunidoc, inputit, "stackfile");
			deleteUniDocValue(inputunidoc, inputit, "frameid");
		}
		stktab.close();
		writeUniDoc(inputunidoc, "prime2D_in.txt");
		delete inputunidoc;
	}
}

void prime2DPost (std::string directory){
	
	std::vector<std::string> 					files;
	std::vector<std::string>::iterator			filesit;
	int											datait;
	bool										newclassdoc;
	UniDoc* 									classdocunidoc;
	UniDoc*										prime2dunidoc;
	std::string									classdocname;
	std::string									cavgsdoc;
	std::string									prime2ddoc;
	
	listFilesInDirectory(files, directory);
	
	prime2dunidoc = new UniDoc();

	
	if(fileExists(directory + "/prime2D_iterations.simple")){
		readUniDoc(prime2dunidoc, directory + "/prime2D_iterations.simple");
	}
	
	for(filesit = files.begin(); filesit != files.end(); filesit++) {
		if(std::string(*filesit).find("classdoc_") != std::string::npos) {
			newclassdoc = true;
			for(datait = 0; datait < prime2dunidoc->data.size(); datait++) {
				if(getUniDocValue(prime2dunidoc, datait, "classdoc", classdocname)) {
					if(classdocname == std::string(*filesit)){
						newclassdoc = false;
					}
				}
			}
			if(newclassdoc){
				classdocunidoc = new UniDoc();
				classdocunidoc->data.push_back("");
				cavgsdoc = std::string(*filesit);
				cavgsdoc.replace(0, 9, "cavgs_iter");
				cavgsdoc.replace(cavgsdoc.end() - 3, cavgsdoc.end(), "mrc");
				addUniDocKeyVal(classdocunidoc, 0, "cavgs", cavgsdoc);
				prime2ddoc = std::string(*filesit);
				prime2ddoc.replace(0, 9, "prime2Ddoc_");
				addUniDocKeyVal(classdocunidoc, 0, "prime2ddoc", prime2ddoc);
				addUniDocKeyVal(classdocunidoc, 0, "classdoc", std::string(*filesit));
				prime2dunidoc->data.push_back(classdocunidoc->data[0]);
				delete classdocunidoc;
			}
		}
	} 
	writeUniDoc(prime2dunidoc, directory + "/prime2D_iterations.simple");
	delete prime2dunidoc;
}

void prime3DPre (std::string simpleinput){
	
	UniDoc*										inputunidoc;
	std::ofstream								stktab;
	int											inputit;
	std::string									stackfile;
	std::string									stackfilepre;
	std::string									inputroot;
	char*										dname;
	char*										dirc;
	std::string									stackfilepath;
	std::string									stackfilelink;
	std::string									stackfilerealpath;
	
	mkdir("particles", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	stktab.open("stkparts.txt");
	dirc = strdup(simpleinput.c_str());
	inputroot = std::string(dirname(dirc));
	
	if(fileExists(simpleinput) && stktab.is_open()){
		inputunidoc = new UniDoc();
		readUniDoc(inputunidoc, simpleinput);
		stackfilepre = "";
		for(inputit = 0; inputit < inputunidoc->data.size(); inputit++) {
			getUniDocValue(inputunidoc, inputit, "stackfile", stackfile);
			if(stackfile != stackfilepre){
				stackfilepath = inputroot + "/" + stackfile;
				stackfilerealpath = std::string(realpath(stackfilepath.c_str(), NULL));
				dirc = strdup(stackfilepath.c_str());
				stktab << "particles/" + std::string(basename(dirc)) << "\n";
				std::cout << stackfilepath << " " << stackfilerealpath << std::endl;
				stackfilelink = "particles/" + std::string(basename(dirc));
				symlink(stackfilerealpath.c_str(), stackfilelink.c_str());
				stackfilepre = stackfile;
			}
			deleteUniDocValue(inputunidoc, inputit, "stackfile");
			deleteUniDocValue(inputunidoc, inputit, "frameid");
		}
		stktab.close();
		writeUniDoc(inputunidoc, "prime3D_in.txt");
		delete inputunidoc;
	}
}

void unblurPre (std::string moviesdirectory){
	
	std::ofstream								movietab;
	std::vector<std::string>					directories;
	std::vector<std::string>					files;
	int											i;
	
	listDirectory(files, directories, moviesdirectory);
	
	movietab.open("movies.txt");
	if(movietab.is_open()){
		for(i = 0; i < files.size(); i++){
			if(files[i].find(".mrc") != std::string::npos){
				movietab << moviesdirectory + "/" + files[i] + "\n";
			}
		}
		movietab.close();
	}
}

void unblurPost (std::string directory){
	
	std::vector<std::string> 					files;
	std::vector<std::string>::iterator			filesit;
	int											datait;
	bool										newmicrograph;
	UniDoc* 									micrographunidoc;
	UniDoc*										preprocunidoc;
	std::string									micrographname;
	std::string									unidocline;
	std::string									attribute;
	char*										dirc;
	
	listFilesInDirectory(files, directory);
	
	preprocunidoc = new UniDoc();
	
	if(fileExists(directory + "/unblur_out.simple")){
		readUniDoc(preprocunidoc, directory + "/unblur_out.simple");
	}
	
	for(filesit = files.begin(); filesit != files.end(); filesit++) {
		if(std::string(*filesit).find("intg.mrc") != std::string::npos) {
			newmicrograph = true;
			for(datait = 0; datait < preprocunidoc->data.size(); datait++) {
				if(getUniDocValue(preprocunidoc, datait, "intg", micrographname)) {
					if(micrographname == std::string(*filesit)){
						newmicrograph = false;
					}
				}
			}
			if(newmicrograph){
				attribute = std::string(*filesit);
				dirc = strdup(attribute.c_str());
				attribute = std::string(basename(dirc));
				unidocline = "intg=" + attribute;
				attribute.replace(attribute.end() - 8, attribute.end(), "thumb.mrc");
				unidocline += " thumb=" + attribute;
				attribute = std::string(*filesit);
				attribute.replace(attribute.end() - 8, attribute.end(), "pspec.mrc");
				unidocline += " pspec=" + attribute;
				attribute = std::string(*filesit);
				attribute.replace(attribute.end() - 8, attribute.end(), "forctf.mrc");
				unidocline += " forctf=" + attribute;
				std::cout << unidocline << std::endl;
				preprocunidoc->data.push_back(unidocline);
			}
		}
	} 
	writeUniDoc(preprocunidoc, directory + "/unblur_out.simple");
	delete preprocunidoc;
}

void externalJob (JSONResponse* response, struct http_message* message) {
	
	pid_t 					pid;
	int						jobid;
	std::string 			jobname;
	std::string 			jobdescription;
	std::string 			jobfolder;
	std::string 			jobtype;
	std::string 			table;
	std::string 			loop;
	std::string 			runstring;
	
	runstring = std::string(message->query_string.p);
	runstring = runstring.substr (0,message->query_string.len);
	
	if(getRequestVariable(message, "jobfolder", jobfolder) && getRequestVariable(message, "jobtype", jobtype) && getRequestVariable(message, "table", table)){
		getRequestVariable(message, "jobname", jobname);
		getRequestVariable(message, "jobdescription", jobdescription);
		getRequestVariable(message, "loop", loop);
		
		addJob(jobname, jobdescription, jobfolder, jobtype, table, runstring, jobid);
	
		pid = fork();
	
		if (pid == 0){
			usleep(1000000);
			updateJobStatus("Running", table, jobid);
			if(jobtype == "preproc"){
				if(loop == "yes"){
					while(true){
						std::cout << "looping" << std::endl;
						preprocPost (jobfolder);
						usleep(60000000);
					}
				}else{
					preprocPost (jobfolder);
				}
			} else if (jobtype == "prime2D"){
				prime2DPost (jobfolder);
			} else if (jobtype == "unblur"){
				unblurPost (jobfolder);
			} else if (jobtype == "extract"){
				chdir(jobfolder.c_str());
				extractPost (jobfolder);
			}
			
			/*if (std::strstr(jobtype, "preproc") ) {
				preprocPost (jobfolder);
			}else if (std::strstr(jobtype, "prime2d") ) {
				prime2DPost (jobfolder);
			}else if (std::strstr(jobtype, "stream") ) {
				streamPost (jobfolder);
		*/	
			updateJobPID(0, table, jobid);
			updateJobStatus ("External", table, jobid);
		} else if (pid > 0) {
			waitpid(pid, 0, WNOHANG);
			updateJobPID(pid, table, jobid);
		} else {
			std::cout << "Fork Failed" << std::endl;
		}
	}
}

void generateSimpleEnv (struct http_message* message){
	
	std::string							element;
	char*								simplepath;
	std::ofstream 						envfile;
	
	simplepath = std::getenv("SIMPLE_PATH");
	
	envfile.open ("simple_distr_config.env");
	
	envfile << "simple_path = " << simplepath << "\n";

	if(getRequestVariable(message, "time_per_image", element)){
		envfile << "time_per_image = " << element << "\n";
	}
	
	if(getRequestVariable(message, "qsys_qos", element)){
		envfile << "qsys_qos = " << element << "\n";
	}

	if(getRequestVariable(message, "queue_submit", element)){
		if(element == "no"){
			envfile << "qsys_name = local" << "\n";
		} else if(getRequestVariable(message, "qsys_name", element)){
			envfile << "qsys_name = " << element << "\n";
		}
	}
	
	if(getRequestVariable(message, "qsys_partition", element)){
		envfile << "qsys_partition = " << element << "\n";
	}
	
	//if(getRequestVariable(message, "nparts", element)){
	envfile << "job_ntasks = 1" << "\n";
	//}
	
	if(getRequestVariable(message, "job_memory_per_task", element)){
		envfile << "job_memory_per_task = " << element << "\n";
	}
	
	if(getRequestVariable(message, "job_name", element)){
		envfile << "job_name = " << element << "\n";
	}

	if(getRequestVariable(message, "job_ntasks_per_socket", element)){
		envfile << "job_ntasks_per_socket = " << element << "\n";
	}
	
	if(getRequestVariable(message, "nthr", element)){
		envfile << "job_cpus_per_task = " << element << "\n";
	}

	envfile.close();
}

void getSimpleJobArgs (std::string& executable, std::string& program, std::vector<std::string>& arguments){
	
	FILE*			stream;
	char 			buffer[1024];
	std::string		command;
	std::string		element;
	int				status;
	
	command = executable + " prg=" + program + " 2>&1";
	
	stream = popen(command.c_str(), "r");
	
	if (stream == NULL){
		//ERROR
	}
	
	while (fgets(buffer, 512, stream) != NULL){
		getElementFromString(std::string(buffer), element, 0);
		if (element != "" && element != "\n" && element != " "){
			arguments.push_back(element);
		}
	}
	status = pclose(stream);
	if (status == -1) {
    //ERROR
	}

}

/*
void simpleSlurmSubmit(std::string& command, std::string& directory, int& jobid, char* table, char* databasepath, char* jobtype){
		std::string 						slurmcommand;
		FILE*								stream;
		
		char *cstr = new char[directory.length() + 1];
		strcpy(cstr, directory.c_str());
		
		slurmcommand = "sbatch -n 1 -c 1 --wrap \"" + command + "\"";
	
		stream = popen(slurmcommand.c_str(), "r");
		pclose(stream);
		
		updatePID (jobid, table, databasepath, "slurm");
		
		while(true){ // TEST SLURM JOB!
			usleep(20000000);
			if(jobtype == "prime2D_stream"){
				streamPost(cstr);
			}else if (jobtype == "prime2D"){
				prime2DPost(cstr);
			}else if (jobtype == "preproc"){
				preprocPost(cstr);
			}
		}
		if(jobtype == "prime2D_stream"){
			streamPost(cstr);
		}else if (jobtype == "prime2D"){
			prime2DPost(cstr);
		}else if (jobtype == "preproc"){
			preprocPost(cstr);
		}
		updatePID (jobid, table, databasepath, "0");
		updateStatus (jobid, table, databasepath, "Finished");
		delete [] cstr;
}
*/

void simpleLocalSubmit(std::string& command, std::string& directory, int& jobid, std::string& table, std::string& jobtype){
	
	pid_t 								pid;
	FILE*								stream;
	
	pid = fork();
	
	if (pid == 0){
		stream = popen(command.c_str(), "r");
		pclose(stream);

	} else if (pid > 0) {
		updateJobPID( pid, table, jobid);
		while(waitpid(pid, 0, WNOHANG) == 0){
			usleep(20000000);
			if(jobtype == "prime2D_stream"){
		//		streamPost(directory);
			} else if (jobtype == "prime2D"){
				prime2DPost(directory);
			}else if (jobtype == "preproc"){
				preprocPost(directory);
			}
		}
		if(jobtype == "prime2D_stream"){
		//	streamPost(directory);
		}else if (jobtype == "prime2D"){
			prime2DPost(directory);
		}else if (jobtype == "preproc"){
			preprocPost(directory);
		}else if (jobtype == "ctffind"){
			ctffindPost(directory);
		}else if (jobtype == "ctffit"){
			ctffitPost(directory);	
		}else if (jobtype == "extract"){
			extractPost(directory);
		}else if (jobtype == "ini3D"){
			ini3DPost(directory);
		}else if (jobtype == "unblur"){
			unblurPost(directory);
		}
		updateJobPID (0, table, jobid);
		updateJobStatus ("Finished", table, jobid);
	} else {
		std::cout << "Fork Failed" << std::endl;
	}	

}

void syncJob (JSONResponse* response, struct http_message* message) {
	
	std::string								program;
	std::string								executable;
	std::string								jobname;
	std::string								jobdescription;
	std::string								jobtype;
	std::string								jobfolder;
	std::string								table;	
	std::string								sourcefolder;
	std::string								destinationfolder;
	std::string								fileidentifier;
	std::string								deletesource;
	std::string								command;
	std::string								outputfolder;
	std::string								runstring;
	FILE*									stream;		
	int										jobid;	
	int										status;
	pid_t 									pid;
	pid_t 									sid;
	
	runstring = std::string(message->query_string.p);
	runstring = runstring.substr (0,message->query_string.len);
	
	if(getRequestVariable(message, "jobtype", jobtype) && getRequestVariable(message, "jobfolder", jobfolder) && getRequestVariable(message, "table", table)){
		
		getRequestVariable(message, "jobname", jobname);
		getRequestVariable(message, "jobdescription", jobdescription);
		
		addJob(jobname, jobdescription, jobfolder, jobtype, table, runstring, jobid);
		getRequestVariable(message, "sourcefolder", sourcefolder);
		getRequestVariable(message, "destinationfolder", destinationfolder);
		getRequestVariable(message, "fileidentifier", fileidentifier);
		getRequestVariable(message, "deletesource", deletesource);
		
		pid = fork();
	
		if (pid == 0){
			umask(0);
			sid = setsid();

			updateJobStatus("Running", table, jobid);
			
			outputfolder = jobfolder + "/" + std::to_string(jobid) + "_" + jobtype;
			
			updateJobFolder(outputfolder, table, jobid);
			
			status = mkdir(outputfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			
			status = chdir(outputfolder.c_str());
			
			command = "microscope_sync.sh -s " + sourcefolder + " -d " + destinationfolder + " -i " + fileidentifier;
			
			if(deletesource == "yes") {
				command += " -r";
			}
			command += " > simple_job.log";
			std::cout << command << std::endl;
			stream = popen(command.c_str(), "r");
			pclose(stream);
			
			updateJobStatus ("Failed", table, jobid);
			exit(0);
			
			
		} else if (pid > 0) {
			waitpid(pid, 0, WNOHANG);
			
		} else {
			std::cout << "Fork Failed" << std::endl;
		}
	}
}

void simpleJob (JSONResponse* response, struct http_message* message) {
	
	std::string								program;
	std::string								executable;
	std::string								jobname;
	std::string								jobdescription;
	std::string								jobtype;
	std::string								jobfolder;
	std::string								table;	
	std::string								outputfolder;
	std::string								command;
	std::string								programarg;
	std::string								argval;
	std::string								argval2;
	std::string								argval3;
	std::string								argval4;
	int										jobid;	
	int										status;	
	pid_t 									pid;
	pid_t 									sid;
	std::vector<std::string> 				arguments;
	std::vector<std::string>::iterator		argit;
	std::string								runstring;
	
	runstring = std::string(message->query_string.p);
	runstring = runstring.substr (0,message->query_string.len);
	
	if(getRequestVariable(message, "jobtype", jobtype) && getRequestVariable(message, "jobfolder", jobfolder) && getRequestVariable(message, "table", table)){
		
		getRequestVariable(message, "jobname", jobname);
		getRequestVariable(message, "jobdescription", jobdescription);
		
		addJob(jobname, jobdescription, jobfolder, jobtype, table, runstring, jobid);
		
		pid = fork();
	
		if (pid == 0){
			umask(0);
			sid = setsid();

			updateJobStatus("Running", table, jobid);
			
			outputfolder = jobfolder + "/" + std::to_string(jobid) + "_" + jobtype;
			
			updateJobFolder(outputfolder, table, jobid);
			
			status = mkdir(outputfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			
			status = chdir(outputfolder.c_str());
			
			if(getRequestVariable(message, "executable", executable) && getRequestVariable(message, "program", program)){

				command = executable + " prg=" + std::string(program);

				getSimpleJobArgs(executable, program, arguments);
			}
			
			for(argit = arguments.begin(); argit != arguments.end(); ++argit) {
				programarg = *argit;
				if (getRequestVariable(message, programarg, argval)){
					if(argval[0] != ' '){
						command += " " + programarg + "=" + argval;
					}
				}
				argval[0] = ' ';
			}
			
			
			if(program == "prime2D" && getRequestVariable(message, "simpleinput", argval)){
					prime2DPre(argval);
			} else if (program == "unblur" && getRequestVariable(message, "moviesdirectory", argval)){
					unblurPre(argval);
					command += " filetab=movies.txt";
			}else if (program == "ctffind" && getRequestVariable(message, "micrographsdirectory", argval)){
					ctffindPre(argval, "");
					command += " filetab=micrographs.txt";
			}else if (program == "ctffind" && getRequestVariable(message, "unbluroutput", argval)){
					ctffindPre("", argval);
					command += " filetab=micrographs.txt";
			}else if (program == "ctffit" && getRequestVariable(message, "micrographsdirectory", argval)){
					ctffitPre(argval, "");
					command += " filetab=micrographs.txt";
			}else if (program == "ctffit" && getRequestVariable(message, "unbluroutput", argval)){
					ctffitPre("", argval);
					command += " filetab=micrographs.txt";
			}else if (program == "extract" && getRequestVariable(message, "boxfilesdirectory", argval) && getRequestVariable(message, "ctffindoutput", argval2)){
					extractPre(argval, "", argval2);
					command += " unidoc=unidoc_in.txt boxtab=boxtab.txt";
			}else if (program == "extract" && getRequestVariable(message, "boxfilesdirectory", argval) && getRequestVariable(message, "micrographsdirectory", argval2)){
					extractPre(argval, argval2, "");
					command += " unidoc=unidoc_in.txt boxtab=boxtab.txt";
			}else if (program == "ini3D_from_cavgs" && getRequestVariable(message, "simpleinput", argval)){
					ini3DPre(argval, command);
			}else if (program == "pick" && getRequestVariable(message, "simpleinput", argval) && getRequestVariable(message, "pickrefs", argval2)){
				getRequestVariable(message, "pcontrast", argval3);
				getRequestVariable(message, "pgrp", argval4);
				pickPre(argval, "", argval2, "", argval4, argval3);
				command += " filetab=picktab.txt refs=pickrefs.mrc";
			}else if (program == "pick" && getRequestVariable(message, "micrographsdirectory", argval) && getRequestVariable(message, "pickrefs", argval2)){
				getRequestVariable(message, "pcontrast", argval3);
				getRequestVariable(message, "pgrp", argval4);
				pickPre("", argval, argval2, "", argval4, argval3);
				command += " filetab=picktab.txt refs=pickrefs.mrc";
			}else if (program == "pick" && getRequestVariable(message, "simpleinput", argval) && getRequestVariable(message, "pickvol", argval2)){
				getRequestVariable(message, "pcontrast", argval3);
				getRequestVariable(message, "pgrp", argval4);
				pickPre(argval, "", "", argval2, argval4, argval3);
				command += " filetab=picktab.txt refs=pickrefs.mrc";
			}else if (program == "pick" && getRequestVariable(message, "micrographsdirectory", argval) && getRequestVariable(message, "pickvol", argval2)){
				getRequestVariable(message, "pcontrast", argval3);
				getRequestVariable(message, "pgrp", argval4);
				pickPre("", argval, "", argval2, argval4, argval3);
				command += " filetab=picktab.txt refs=pickrefs.mrc";
			}else if(program == "prime3D" && getRequestVariable(message, "simpleinput", argval)){
					prime3DPre(argval);
				}
			std::cout << command << std::endl;
			command += " >> simple_job.log";
			
			generateSimpleEnv(message);
			
			//if(getRequestVariable(message, "submitmaster", argval)){
				
			//}else{
				std::cout << command << std::endl;
				simpleLocalSubmit(command, outputfolder, jobid, table, jobtype);
			//}
			exit(0);
			
		} else if (pid > 0) {
			waitpid(pid, 0, WNOHANG);
			
		} else {
			std::cout << "Fork Failed" << std::endl;
		}
	}
}

void getJobs (JSONResponse* response, struct http_message* message) {
	
	std::string 			databasepath;
	std::string 			table;
	
	if(getRequestVariable(message, "table", table)){
		getDatabaseLocation(databasepath);
		SQLQuery(databasepath, response->jobs, "CREATE TABLE IF NOT EXISTS " + table + " (jobid integer PRIMARY KEY, jobname text NOT NULL, jobdescription text NOT NULL, jobfolder text NOT NULL, jobtype text NOT NULL, jobstatus text NOT NULL, jobpid integer NOT NULL, jobrerun text NOT NULL);");									
		SQLQuery(databasepath, response->jobs, "SELECT * FROM " + table);	
	}								
}

void killJob (JSONResponse* response, struct http_message* message) {
	
	std::string 				pid;
	std::string 				table;
	std::string 				jobid;
	std::string 				jobfolder;
	std::vector<std::string>	files;
	int							fileit;
	
	if(getRequestVariable(message, "table", table) && getRequestVariable(message, "pid", pid) && getRequestVariable(message, "jobid", jobid) && getRequestVariable(message, "jobfolder", jobfolder)){
		listFilesInDirectory(files, jobfolder);
		for(fileit = 0; fileit < files.size(); fileit++){
		//	if(files[fileit].find(".pid") != std::string::npos
		}
		
	}
}

void newProject (JSONResponse* response, struct http_message* message) {
	
	std::string 			databasepath;
	std::string				projectname;
	std::string				projectdescription;
	std::string				projectfolder;
	std::string  			projecttable; 
	
	if (getRequestVariable(message, "projectname", projectname) && getRequestVariable(message, "projectfolder", projectfolder)) {
		if(fileExists(projectfolder)){
			getRequestVariable(message, "projectdescription", projectdescription);
			getDatabaseLocation(databasepath);
			srand((int)time(0));
			projecttable = "p" + std::to_string(rand());
			SQLQuery(databasepath, response->projects, "CREATE TABLE IF NOT EXISTS projects (projectid integer PRIMARY KEY, projectname text NOT NULL, projectdescription text NOT NULL, projectfolder text NOT NULL, projecttable text NOT NULL, projectusage bool NOT NULL);");	
			SQLQuery(databasepath, response->projects, "INSERT INTO projects(projectname, projectdescription, projectfolder, projecttable, projectusage) VALUES('"	+ projectname + "','" + projectdescription + "','" + projectfolder + "','" + projecttable + "','" + "true" + "');");
		} else {
			response->error = "Project folder does not exist!";
		}
	}
	
}

void getProjects (JSONResponse* response) {

	std::string  														sqlresult;
	std::string															databasepath;
	
	getDatabaseLocation(databasepath);
	
	SQLQuery(databasepath, response->projects, "SELECT * FROM projects");

}

void deleteProject(JSONResponse* response, struct http_message* message) {
	
	std::string				databasepath;
	std::string				table;
	std::string				projectid;
	
	getDatabaseLocation(databasepath);
	
	if (getRequestVariable(message, "projectid", projectid)) {
		SQLQuery(databasepath, "DELETE FROM projects WHERE projectid=" + projectid + ";");
	}
}

void viewCtffind (JSONResponse* response, struct http_message* message) {
	
	std::string 							directory;
	UniDoc*									unidoc;
	int										datait;
	std::map<std::string, std::string>		micrographmap;
	std::string								value;
	int										randomint;
	
	
	if (getRequestVariable(message, "folder", directory)){
		response->rootdirectory = directory;
		if(fileExists(directory + "/ctffind_out.simple")){
			unidoc = new UniDoc();
			readUniDoc(unidoc, directory + "/ctffind_out.simple");
			response->inputfilename = "ctffind_out.simple.";
			for(datait = 0; datait < unidoc->data.size(); datait++){
				micrographmap.clear();
				micrographmap["id"] = std::to_string(datait);
				getUniDocValue(unidoc, datait, "intg", value);
				micrographmap["intg"] = value;
				getUniDocValue(unidoc, datait, "movie", value);
				micrographmap["movie"] = value;
				getUniDocValue(unidoc, datait, "pspec", value);
				micrographmap["pspec"] = value;
				getUniDocValue(unidoc, datait, "thumb", value);
				micrographmap["thumb"] = value;
				getUniDocValue(unidoc, datait, "boxfile", value);
				micrographmap["state"] = value;
				getUniDocValue(unidoc, datait, "smpd", value);
				micrographmap["smpd"] = value;
				getUniDocValue(unidoc, datait, "kv", value);
				micrographmap["kv"] = value;
				getUniDocValue(unidoc, datait, "cs", value);
				micrographmap["cs"] = value;
				getUniDocValue(unidoc, datait, "fraca", value);
				micrographmap["fraca"] = value;
				getUniDocValue(unidoc, datait, "dfx", value);
				micrographmap["dfx"] = value;
				getUniDocValue(unidoc, datait, "dfy", value);
				micrographmap["dfy"] = value;
				getUniDocValue(unidoc, datait, "angast", value);
				micrographmap["angast"] = value;
				getUniDocValue(unidoc, datait, "ctfres", value);
				micrographmap["ctfres"] = value;
				getUniDocValue(unidoc, datait, "ctffindfit", value);
				micrographmap["ctffindfit"] = value;
				response->snapshots.push_back(micrographmap);
			}
			delete unidoc;
		}
	}
	
}

void viewCtffit (JSONResponse* response, struct http_message* message) {
	
	std::string 							directory;
	UniDoc*									unidoc;
	int										datait;
	std::map<std::string, std::string>		micrographmap;
	std::string								value;
	int										randomint;
	
	
	if (getRequestVariable(message, "folder", directory)){
		response->rootdirectory = directory;
		if(fileExists(directory + "/ctffit_out.simple")){
			unidoc = new UniDoc();
			readUniDoc(unidoc, directory + "/ctffit_out.simple");
			response->inputfilename = "ctffit_out.simple";
			for(datait = 0; datait < unidoc->data.size(); datait++){
				micrographmap.clear();
				micrographmap["id"] = std::to_string(datait);
				getUniDocValue(unidoc, datait, "intg", value);
				micrographmap["intg"] = value;
				getUniDocValue(unidoc, datait, "movie", value);
				micrographmap["movie"] = value;
				getUniDocValue(unidoc, datait, "pspec", value);
				micrographmap["pspec"] = value;
				getUniDocValue(unidoc, datait, "thumb", value);
				micrographmap["thumb"] = value;
				getUniDocValue(unidoc, datait, "boxfile", value);
				micrographmap["state"] = value;
				getUniDocValue(unidoc, datait, "smpd", value);
				micrographmap["smpd"] = value;
				getUniDocValue(unidoc, datait, "kv", value);
				micrographmap["kv"] = value;
				getUniDocValue(unidoc, datait, "cs", value);
				micrographmap["cs"] = value;
				getUniDocValue(unidoc, datait, "fraca", value);
				micrographmap["fraca"] = value;
				getUniDocValue(unidoc, datait, "dfx", value);
				micrographmap["dfx"] = value;
				getUniDocValue(unidoc, datait, "dfy", value);
				micrographmap["dfy"] = value;
				getUniDocValue(unidoc, datait, "angast", value);
				micrographmap["angast"] = value;
				getUniDocValue(unidoc, datait, "ctfres", value);
				micrographmap["ctfres"] = value;
				getUniDocValue(unidoc, datait, "ctffindfit", value);
				micrographmap["ctffindfit"] = value;
				response->snapshots.push_back(micrographmap);
			}
			delete unidoc;
		}
	}
	
}

void viewExtract (JSONResponse* response, struct http_message* message) {
	
	std::string 							directory;
	UniDoc*									iterationsunidoc;
	int										datait;
	std::map<std::string, std::string>		iterationmap;
	std::string								value;
	
	if (getRequestVariable(message, "folder", directory)){
		response->rootdirectory = directory;
		if(fileExists(directory + "/extract_parts.txt")){
			iterationsunidoc = new UniDoc();
			readUniDoc(iterationsunidoc, directory + "/extract_parts.txt");
			for(datait = 0; datait < iterationsunidoc->data.size(); datait++){
				iterationmap.clear();
				getUniDocValue(iterationsunidoc, datait, "stackfile", value);
				iterationmap["stackfile"] = value;
				getUniDocValue(iterationsunidoc, datait, "particlecount", value);
				iterationmap["particlecount"] = value;
				response->iterations.push_back(iterationmap);
			}
			delete iterationsunidoc;
		}
		if(fileExists(directory + "/extract_out.simple")){
			iterationsunidoc = new UniDoc();
			readUniDoc(iterationsunidoc, directory + "/extract_out.simple");
			response->particlecount = std::to_string(iterationsunidoc->data.size());
			delete iterationsunidoc;
		}
	}

}

void viewPrime2D (JSONResponse* response, struct http_message* message) {
	
	std::string 							directory;
	UniDoc*									iterationsunidoc;
	int										datait;
	std::map<std::string, std::string>		iterationmap;
	std::string								value;
	
	if (getRequestVariable(message, "folder", directory)){
		response->rootdirectory = directory;
		if(fileExists(directory + "/prime2D_iterations.simple")){
			iterationsunidoc = new UniDoc();
			readUniDoc(iterationsunidoc, directory + "/prime2D_iterations.simple");
			for(datait = 0; datait < iterationsunidoc->data.size(); datait++){
				iterationmap.clear();
				getUniDocValue(iterationsunidoc, datait, "cavgs", value);
				iterationmap["cavgs"] = value;
				getUniDocValue(iterationsunidoc, datait, "classdoc", value);
				iterationmap["classdoc"] = value;
				getUniDocValue(iterationsunidoc, datait, "prime2ddoc", value);
				iterationmap["prime2ddoc"] = value;
				response->iterations.push_back(iterationmap);
			}
			delete iterationsunidoc;
		}
	}
	
}

void viewPrime2DIteration (JSONResponse* response, struct http_message* message) {
	
	std::string 							classdoc;
	std::string 							prime2ddoc;
	std::string 							cavgs;
	UniDoc*									iterationunidoc;
	int										datait;
	std::map<std::string, std::string>		iterationmap;
	std::string								value;
	
	if (getRequestVariable(message, "classdoc", classdoc) && getRequestVariable(message, "prime2ddoc", prime2ddoc) && getRequestVariable(message, "cavgs", cavgs)){
		response->prime2ddoc = prime2ddoc;
		if(fileExists(classdoc)){
			iterationunidoc = new UniDoc();
			readUniDoc(iterationunidoc, classdoc);
			for(datait = 0; datait < iterationunidoc->data.size(); datait++){
				iterationmap.clear();
				iterationmap["frameid"] = std::to_string(datait);
				iterationmap["cavgs"] = cavgs;
				getUniDocValue(iterationunidoc, datait, "class", value);
				iterationmap["class"] = value;
				getUniDocValue(iterationunidoc, datait, "res", value);
				iterationmap["res"] = value;
				getUniDocValue(iterationunidoc, datait, "pop", value);
				iterationmap["pop"] = value;
				getUniDocValue(iterationunidoc, datait, "corr", value);
				iterationmap["corr"] = value;
				getUniDocValue(iterationunidoc, datait, "w", value);
				iterationmap["w"] = value;
				response->snapshots.push_back(iterationmap);
			}
			delete iterationunidoc;
		}
	}
	
}

void viewPrime2DParticles (JSONResponse* response, struct http_message* message) {
	
	std::string 							prime2ddoc;
	std::string								classnumber;
	UniDoc*									unidoc;
	int										datait;
	std::map<std::string, std::string>		iterationmap;
	std::string								value;
	
	if (getRequestVariable(message, "prime2ddoc", prime2ddoc) && getRequestVariable(message, "class", classnumber)){
		if(fileExists(prime2ddoc)){
			unidoc = new UniDoc();
			readUniDoc(unidoc, prime2ddoc);
			addUnidocParticleIdentifiers(unidoc, prime2ddoc);
			for(datait = 0; datait < unidoc->data.size(); datait++){
				iterationmap.clear();
				getUniDocValue(unidoc, datait, "class", value);
				if(value == classnumber){
					getUniDocValue(unidoc, datait, "stackfile", value);
					iterationmap["stackfile"] = value;
					getUniDocValue(unidoc, datait, "frameid", value);
					iterationmap["frameid"] = value;
					response->snapshots.push_back(iterationmap);
				}
			}
			delete unidoc;
		}
	}
	
}

void savePrime2DSelection(JSONResponse* response, struct http_message* message) { //fix this
	
	std::string					inputfilename;
	std::string					outputfilename;
	std::string					selection;
	std::string					selectedframe;
	MRCFile*					mrcfile;
	int							selectionit;
	std::ofstream				output;
	int							framecount;

	if(getRequestVariable(message, "inputfilename", inputfilename) && getRequestVariable(message, "outputfilename", outputfilename) && getRequestVariable(message, "selection", selection)) {
		if (inputfilename.find(".mrc") != std::string::npos && fileExists(inputfilename)){
			output.open(outputfilename.c_str(), std::ifstream::binary);
			if(output.is_open()){
				mrcfile = new MRCFile();
				readMRCFrame(mrcfile, inputfilename, 0);
				mrcfile->extendedheader = 0;
				memcpy(&mrcfile->header[92], &mrcfile->extendedheader, 4);
				output.write(mrcfile->header, 1024);
				selectionit = selection.find(",");
				framecount = 0;
				while(selectionit != std::string::npos){
					selectedframe = selection.substr(0, selectionit);
					readMRCFrame(mrcfile, inputfilename, std::stoi(selectedframe));
					output.write((char*)mrcfile->data, mrcfile->nx * mrcfile->ny * 4);
					//output.write(reinterpret_cast<const char*>(&mrcfile->data), sizeof(mrcfile->data));
					selection.erase(0,selectionit + 1);
					selectionit = selection.find(",");
					framecount++;
				}
				output.seekp(8, std::fstream::beg);
				
				output.write(reinterpret_cast<const char*>(&framecount),4);
				output.close();
			}	
		}
	}
}

void savePrime2DSelectionParticles(JSONResponse* response, struct http_message* message) {
	
	std::string					inputfilename;
	std::string					outputfilename;
	std::string					inverseselection;
	UniDoc*						unidoc;
	int							selectionit;
	std::string					removeclass;
	std::string					classnumber;
	int							datait;
	
	if(getRequestVariable(message, "inputfilename", inputfilename) && getRequestVariable(message, "outputfilename", outputfilename)) {
		getRequestVariable(message, "inverseselection", inverseselection);
		unidoc = new UniDoc();
		readUniDoc(unidoc, inputfilename);
		selectionit = inverseselection.find(",");
		while(selectionit != std::string::npos){
			removeclass = inverseselection.substr(0,selectionit);
			std::cout << removeclass << std::endl;
			for(datait = 0; datait < unidoc->data.size(); datait++){
				getUniDocValue(unidoc, datait, "class", classnumber);
				if(classnumber == removeclass){
					addUniDocKeyVal(unidoc, datait, "state", "0");
				}
			}
			inverseselection.erase(0,selectionit + 1);
			selectionit = inverseselection.find(",");
		}
		addUnidocParticleIdentifiers(unidoc, inputfilename);
		writeUniDoc(unidoc, outputfilename);
		writeRelionStar(unidoc, outputfilename);
		delete unidoc;
	}
}

void viewPreproc (JSONResponse* response, struct http_message* message) {
	
	std::string 							directory;
	UniDoc*									unidoc;
	int										datait;
	std::map<std::string, std::string>		micrographmap;
	std::string								value;
	int										randomint;
	
	
	if (getRequestVariable(message, "folder", directory)){
		response->rootdirectory = directory;
		if(fileExists(directory + "/preproc_out.simple")){
			srand((int)time(0));
			randomint = rand();
			unidoc = new UniDoc();
			readUniDoc(unidoc, directory + "/preproc_out.simple");
			writeUniDoc(unidoc, directory + "/preproc_out.simple." + std::to_string(randomint) + "~");
			response->inputfilename = "preproc_out.simple." + std::to_string(randomint) + "~";
			for(datait = 0; datait < unidoc->data.size(); datait++){
				micrographmap.clear();
				micrographmap["id"] = std::to_string(datait);
				getUniDocValue(unidoc, datait, "intg", value);
				micrographmap["intg"] = value;
				getUniDocValue(unidoc, datait, "movie", value);
				micrographmap["movie"] = value;
				getUniDocValue(unidoc, datait, "pspec", value);
				micrographmap["pspec"] = value;
				getUniDocValue(unidoc, datait, "thumb", value);
				micrographmap["thumb"] = value;
				getUniDocValue(unidoc, datait, "boxfile", value);
				micrographmap["boxfile"] = value;
				getUniDocValue(unidoc, datait, "state", value);
				micrographmap["state"] = value;
				getUniDocValue(unidoc, datait, "smpd", value);
				micrographmap["smpd"] = value;
				getUniDocValue(unidoc, datait, "kv", value);
				micrographmap["kv"] = value;
				getUniDocValue(unidoc, datait, "cs", value);
				micrographmap["cs"] = value;
				getUniDocValue(unidoc, datait, "fraca", value);
				micrographmap["fraca"] = value;
				getUniDocValue(unidoc, datait, "dfx", value);
				micrographmap["dfx"] = value;
				getUniDocValue(unidoc, datait, "dfy", value);
				micrographmap["dfy"] = value;
				getUniDocValue(unidoc, datait, "angast", value);
				micrographmap["angast"] = value;
				getUniDocValue(unidoc, datait, "ctfres", value);
				micrographmap["ctfres"] = value;
				getUniDocValue(unidoc, datait, "nptcls", value);
				micrographmap["nptcls"] = value;
				getUniDocValue(unidoc, datait, "ctffindfit", value);
				micrographmap["ctffindfit"] = value;
				response->snapshots.push_back(micrographmap);
			}
			delete unidoc;
		}
		
		if(fileExists(directory + "/preproc_out_selected.simple")){
			
			unidoc = new UniDoc();
			readUniDoc(unidoc, directory + "/preproc_out_selected.simple");
			for(datait = 0; datait < unidoc->data.size(); datait++){
				micrographmap.clear();
				getUniDocValue(unidoc, datait, "state", value);
				response->snapshots[datait]["state"] = value;
				response->snapshots[datait]["viewed"] = "1";
			}
			delete unidoc;
		}
		
	}
	
}

void savePreprocSelection(JSONResponse* response, struct http_message* message) {
	
	std::string					inputfilename;
	std::string					outputfilename;
	std::string					inverseselection;
	UniDoc*						unidoc;
	int							selectionit;
	std::string					removeline;
	int							datacount;
	
	if(getRequestVariable(message, "inputfilename", inputfilename) && getRequestVariable(message, "outputfilename", outputfilename)) {
		getRequestVariable(message, "inverseselection", inverseselection);
		std::cout << inverseselection << std::endl;
		unidoc = new UniDoc();
		if(fileExists(inputfilename)){
			readUniDoc(unidoc, inputfilename);
			for(datacount = 0; datacount < unidoc->data.size(); datacount++){
				addUniDocKeyVal(unidoc, datacount, "state", "1");
			}
			
			selectionit = inverseselection.find(",");
			while(selectionit != std::string::npos){
				removeline = inverseselection.substr(0,selectionit);
				addUniDocKeyVal(unidoc, std::stoi(removeline), "state", "0");
				inverseselection.erase(0,selectionit + 1);
				selectionit = inverseselection.find(",");
			}
			
			writeUniDoc(unidoc, outputfilename);
			delete unidoc;
		}else{
			std::cout << inputfilename + "~ Does not exist" << std::endl;
		}
	}
}

void savePreprocSelectionParticles(JSONResponse* response, struct http_message* message) {
	
	std::string					inputfilename;
	std::string					outputfilename;
	std::string					inverseselection;
	UniDoc*						unidoc;
	UniDoc*						partunidoc;
	UniDoc*						outputunidoc;
	UniDoc*						relionstar;
	int							selectionit;
	std::string					removeline;
	int							particleit;
	std::string					state;
	std::string					intg;
	std::string					stackfile;
	std::string					extractparams;
	std::string					value;
	char*						dirc;
	std::string					directory;
	float						defocusU;
	float						defocusV;

	if(getRequestVariable(message, "inputfilename", inputfilename) && getRequestVariable(message, "outputfilename", outputfilename)) {
		getRequestVariable(message, "inverseselection", inverseselection);
		unidoc = new UniDoc();
		if(fileExists(inputfilename)){
			readUniDoc(unidoc, inputfilename);
			selectionit = inverseselection.find(",");
			while(selectionit != std::string::npos){
				removeline = inverseselection.substr(0,selectionit);
				addUniDocKeyVal(unidoc, std::stoi(removeline), "state", "0");
				inverseselection.erase(0,selectionit + 1);
				selectionit = inverseselection.find(",");
			}	
			
			dirc = strdup(inputfilename.c_str());
			directory = std::string(dirname(dirc));
			
			outputunidoc = new UniDoc();
			
			for(selectionit = 0; selectionit < unidoc->data.size(); selectionit++){
				getUniDocValue(unidoc, selectionit, "state", state);
				if(state != "0"){
					getUniDocValue(unidoc, selectionit, "intg", intg);
					stackfile = intg;
					//stackfile.replace(0, 21,  directory + "/pipeline/particles/ptcls_from_");
					stackfile.replace(0, 21, "pipeline/particles/ptcls_from_");
					stackfile.replace(stackfile.end() - 9, stackfile.end(), ".mrc");
					extractparams = intg;
					extractparams.replace(0, 21, directory + "/pipeline/particles/extract_params_");
					extractparams.replace(extractparams.end() - 9, extractparams.end(), ".txt");

					std::cout << stackfile << " " << extractparams << std::endl;
					if(fileExists(stackfile) && fileExists(extractparams)){
						partunidoc = new UniDoc();
						readUniDoc(partunidoc, extractparams);

						for(particleit = 0; particleit <  partunidoc->data.size(); particleit++){
							addUniDocKeyVal(partunidoc, particleit, "stackfile", stackfile);
							addUniDocKeyVal(partunidoc, particleit, "frameid", std::to_string(particleit + 1));
							outputunidoc->data.push_back(partunidoc->data[particleit]);
						}
						delete partunidoc;
					}
				}
			}
			std::cout << outputfilename << std::endl;
			writeUniDoc(outputunidoc, outputfilename);
			writeRelionStar(outputunidoc, outputfilename);
			delete outputunidoc;
		}
		delete unidoc;
	}
}

void viewManualPick (JSONResponse* response, struct http_message* message) {
	
	std::string 							unidocfilename;
	UniDoc*									unidoc;
	int										datait;
	std::map<std::string, std::string>		micrographmap;
	std::string								value;
	char*									dname;
	char*									dirc;
	std::string 							jobtype;
	std::string 							jobfolder;
	std::string 							table;
	std::string 							jobname;
	std::string 							jobdescription;
	std::string 							boxfile;
	std::string								runstring;
	int										jobid;
	int										status;
	std::string 							outputfolder;	
			
	runstring = std::string(message->query_string.p);
	runstring = runstring.substr (0,message->query_string.len);
							
	if(getRequestVariable(message, "jobtype", jobtype) && getRequestVariable(message, "jobfolder", jobfolder) && getRequestVariable(message, "table", table)){
		
		getRequestVariable(message, "jobname", jobname);
		getRequestVariable(message, "jobdescription", jobdescription);
		
		addJob(jobname, jobdescription, jobfolder, jobtype, table, runstring, jobid);
		
		outputfolder = jobfolder + "/" + std::to_string(jobid) + "_" + jobtype;
			
		updateJobFolder(outputfolder, table, jobid);
			
		status = mkdir(outputfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		
		response->jobfolder = outputfolder;
		
		if (getRequestVariable(message, "unbluroutput", unidocfilename)){
			dirc = strdup(unidocfilename.c_str());
			dname = dirname(dirc);
			response->rootdirectory = std::string(dname);
			if(fileExists(unidocfilename)){
				unidoc = new UniDoc();
				readUniDoc(unidoc, unidocfilename);
				for(datait = 0; datait < unidoc->data.size(); datait++){
					micrographmap.clear();
					micrographmap["id"] = std::to_string(datait);
					getUniDocValue(unidoc, datait, "intg", value);
					micrographmap["micrograph"] = value;
					value.replace(value.end() - 3, value.end(), "box");
					dirc = strdup(value.c_str());
					boxfile = outputfolder + "/" + std::string(basename(dirc));
					//getUniDocValue(unidoc, datait, "boxfile", value);
					//if(value.size() > 0){
					//	micrographmap["boxfile"] = value;
					//} else {
						micrographmap["boxfile"] = boxfile;
					//}
					response->snapshots.push_back(micrographmap);
				}
				delete unidoc;
			}
			
		}
	}
}

void viewUnblur (JSONResponse* response, struct http_message* message) {
	
	std::string 							directory;
	UniDoc*									unidoc;
	int										datait;
	std::map<std::string, std::string>		micrographmap;
	std::string								value;
	int										randomint;
	
	
	if (getRequestVariable(message, "folder", directory)){
		response->rootdirectory = directory;
		if(fileExists(directory + "/unblur_out.simple")){
			unidoc = new UniDoc();
			readUniDoc(unidoc, directory + "/unblur_out.simple");
			response->inputfilename = directory + "/unblur_out.simple";
			for(datait = 0; datait < unidoc->data.size(); datait++){
				micrographmap.clear();
				micrographmap["id"] = std::to_string(datait);
				getUniDocValue(unidoc, datait, "intg", value);
				micrographmap["intg"] = value;
				//getUniDocValue(unidoc, datait, "movie", value);
				//micrographmap["movie"] = value;
				getUniDocValue(unidoc, datait, "pspec", value);
				micrographmap["pspec"] = value;
				getUniDocValue(unidoc, datait, "thumb", value);
				micrographmap["thumb"] = value;
				response->snapshots.push_back(micrographmap);
			}
			delete unidoc;
		}
		
	}
	
}

void saveUnblurSelection(JSONResponse* response, struct http_message* message) {
	
	std::string					inputfilename;
	std::string					outputfilename;
	std::string					inverseselection;
	UniDoc*						unidoc;
	int							selectionit;
	std::string					removeline;
	int							datacount;
	
	if(getRequestVariable(message, "inputfilename", inputfilename) && getRequestVariable(message, "outputfilename", outputfilename) && getRequestVariable(message, "inverseselection", inverseselection)) {
		unidoc = new UniDoc();
		if(fileExists(inputfilename)){
			readUniDoc(unidoc, inputfilename);
			for(datacount = 0; datacount < unidoc->data.size(); datacount++){
				addUniDocKeyVal(unidoc, datacount, "state", "1");
			}
			
			selectionit = inverseselection.find(",");
			while(selectionit != std::string::npos){
				removeline = inverseselection.substr(0,selectionit);
				addUniDocKeyVal(unidoc, std::stoi(removeline), "state", "0");
				inverseselection.erase(0,selectionit + 1);
				selectionit = inverseselection.find(",");
			}
			
			writeUniDoc(unidoc, outputfilename);
			delete unidoc;
		}else{
			std::cout << inputfilename + "~ Does not exist" << std::endl;
		}
	}
}

void getPixelsFromMRC(JPEGResponse* response, std::string filename, struct http_message* message){
	
	std::string			contrast;
	std::string			brightness;
	std::string			frameid;
	MRCFile*			mrcfile;
	float				datamean;
	int					datalength;
	int					nx;
	int					ny;
	int					datacount;
	float 				datadiff2;
	float 				datasd;
	float 				datanorm;
	
	if (getRequestVariable(message, "contrast", contrast) && getRequestVariable(message, "brightness", brightness) && getRequestVariable(message, "frameid", frameid)){
		mrcfile = new MRCFile();
		readMRCFrame(mrcfile, filename, std::stoi(frameid));

		datamean = 0;

		memcpy(&nx, &mrcfile->header[0], 4);
		memcpy(&ny, &mrcfile->header[4], 4);
		
		datalength = nx * ny;
		
		for(datacount = 0; datacount < datalength; datacount++){
			datamean += mrcfile->data[datacount];
		}	

		datamean = datamean / datalength;
		datadiff2 = 0;

		for(datacount = 0; datacount < datalength; datacount++){
			datadiff2 += pow((mrcfile->data[datacount] - datamean), 2);
		}
		
		datasd = sqrt(datadiff2 /datalength);

		response->pixels = new unsigned char[datalength];

		for(datacount = 0; datacount < datalength; datacount++){
			datanorm = (mrcfile->data[datacount] - datamean) / datasd;
			datanorm *= 128 / (10.5 - std::stof(contrast));
			datanorm += std::stoi(brightness);
			if(datanorm > 254){
				response->pixels[datacount] = (int) 254;
			} else if(datanorm < 0){
				response->pixels[datacount] = (int) 0;
			} else {
				response->pixels[datacount] = round(datanorm);
			}
		}
		
		response->xdim = nx;
		response->ydim = ny;
		
		delete [] mrcfile->data;
		delete mrcfile;
	}
	
}

void clearBoxes (JSONResponse* response, struct http_message* message) {

	std::string								filename;

	if (getRequestVariable(message, "filename", filename)){
		if(fileExists(filename)){
			std::remove(filename.c_str());
		}
	}
}

void getBoxes (JSONResponse* response, struct http_message* message) {

	std::ifstream 							input;
	std::string								filename;
	std::string 							line;
	std::string 							element;
	std::map<std::string, std::string>		box;
	
	if (getRequestVariable(message, "filename", filename)){
		input.open(filename.c_str());
		if(input.is_open()){
			while(getline(input, line)){
				box.clear();
				getElementFromString(line, element, 0);
				box["x"] = element;
				getElementFromString(line, element, 1);
				box["y"] = element;
				getElementFromString(line, element, 2);
				box["size"] = element;
				response->boxes.push_back(box);
			}
			input.close();
		}
	}	
}

void getBoxes (JSONResponse* response, std::string filename) {

	std::ifstream 							input;
	std::string 							line;
	std::string 							element;
	std::map<std::string, std::string>		box;
	
	input.open(filename.c_str());
	if(input.is_open()){
		while(getline(input, line)){
			box.clear();
			getElementFromString(line, element, 0);
			box["x"] = element;
			getElementFromString(line, element, 1);
			box["y"] = element;
			getElementFromString(line, element, 2);
			box["size"] = element;
			response->boxes.push_back(box);
		}
		input.close();
	}	
}

void getLogFile (JSONResponse* response, struct http_message* message) {
	
	std::string								folder;
	std::string								logfilename;
	std::ifstream 							input;
	std::string 							line;
//	char									finalchar;
//	int										finalascii;
	
	if (getRequestVariable(message, "folder", folder)){
		logfilename = folder + "/simple_job.log";
		if(fileExists(logfilename)){
			input.open(logfilename.c_str(), std::ios::binary);
			if(input.is_open()){
				while(getline(input, line)){
					//line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
					std::replace(line.begin(), line.end(), '\r', '@');
					response->logfile += line + "@";
				}
				input.close();
				std::replace( response->logfile.begin(), response->logfile.end(), '\\', '^');
			}
		
		}
	}
	
}

void autoPick (JSONResponse* response, struct http_message* message){
	
	std::string		pickrefs;
	std::string		pickvol;
	std::string		pgrp;
	std::string		boxfile;
	std::string		folder;
	std::string		micrograph;
	std::string		command;
	std::string		thres;
	std::string		ndev;
	std::string		picktab;
	std::string		smpd;
	std::string		pcontrast;
	
	FILE*								stream;
	int				status;
	std::ofstream					filetab;
	char*			basec;
	
	if (getRequestVariable(message, "folder", folder) && getRequestVariable(message, "micrograph", micrograph) && getRequestVariable(message, "pgrp", pgrp)){
		if(getRequestVariable(message, "pickrefs", pickrefs) && getRequestVariable(message, "pcontrast", pcontrast) && fileExists(folder)){
			command = "cd " + folder + " && simple_exec prg=makepickrefs pgrp=" + pgrp + " stk=" + pickrefs + " nthr=1 pcontrast=" + pcontrast + " >> simple_job.log";
			std::cout << command << std::endl;
			stream = popen(command.c_str(), "r");
			pclose(stream);
		} else if (getRequestVariable(message, "pickvol", pickvol) && fileExists(folder)){
			command = "cd " + folder + " && simple_exec prg=makepickrefs pgrp=" + pgrp + " vol1=" + pickrefs + " nthr=1 pcontrast=" + pcontrast + " >> simple_job.log";
			stream = popen(command.c_str(), "r");
			pclose(stream);
		}
		
		pickrefs = folder + "/pickrefs.mrc";
		picktab = folder + "/picktab.txt";

		if(fileExists(pickrefs)){
			filetab.open(picktab.c_str());
			if(filetab.is_open()){
				basec = strdup(micrograph.c_str());
				boxfile = folder + "/";
				boxfile += std::string(basename(basec));
				boxfile.replace(boxfile.end() - 3, boxfile.end(), "box");
				if(fileExists(boxfile)){
					remove(boxfile.c_str());
				}
				filetab << micrograph << "\n";
				filetab.close();
				command = "cd " + folder + " && simple_exec prg=pick refs=pickrefs.mrc filetab=picktab.txt nthr=4";
				if(getRequestVariable(message, "thres", thres)){
					command += " thres="+thres;
				}
				if(getRequestVariable(message, "ndev", ndev)){
					command += " ndev="+ndev;
				}
				if(getRequestVariable(message, "smpd", smpd)){
					command += " smpd="+smpd;
				}
				command += " >> simple_job.log";
				std::cout << command << std::endl;
				std::cout << boxfile << std::endl;
				stream = popen(command.c_str(), "r");
				pclose(stream);
				
				getBoxes(response, boxfile);
			}
		}
	}
}

void circularScore(MRCFile* box, int particlediameter, int boxsize, float& score){
	
	float 			r;
	int				xit;
	int				yit;
	
	score = 0;
	for(yit = 0; yit < boxsize; yit++){
		for(xit = 0; xit < boxsize; xit++){
			r = sqrt(pow((xit - (boxsize/2)), 2) + pow((yit - (boxsize/2)), 2));
			if(r <= (particlediameter/2)){
				score += box->data[(yit * boxsize) + xit];
			}
		}
	}
}

void refinePick (JSONResponse* response, struct http_message* message){
	
	std::string		micrographname;
	std::string		xcoord;
	std::string		ycoord;
	std::string		boxsize;
	std::string		particlediameter;
	MRCFile*		micrograph;
	MRCFile*		box;
	int 			searchsize = 30;
	int				xit;
	int				yit;
	float 			score;
	int				i;
	std::vector<std::map<std::string, float> >		scores;
	std::map<std::string, float>					scoreelement;
	std::map<std::string, std::string>				returnbox;
	
	getRequestVariable(message, "micrograph", micrographname);
	getRequestVariable(message, "xcoord", xcoord);
	getRequestVariable(message, "ycoord", ycoord);
	getRequestVariable(message, "boxsize", boxsize);
	getRequestVariable(message, "particlediameter", particlediameter);
	
	micrograph = new MRCFile();
	box = new MRCFile();
	
	readMRCFrame(micrograph, micrographname, 0);
	
	std::cout << xcoord << " " << ycoord << std::endl;
	
	for(yit = -searchsize; yit <=  searchsize; yit++){
		for(xit = -searchsize; xit <=  searchsize; xit++){
			extractMRCBox(micrograph, box, std::atoi(xcoord.c_str()) - (std::atoi(boxsize.c_str())/2) + xit, std::atoi(ycoord.c_str()) - (std::atoi(boxsize.c_str())/2) + yit, std::atoi(boxsize.c_str()));
			circularScore(box, std::atoi(particlediameter.c_str()), std::atoi(boxsize.c_str()), score);
			scoreelement.clear();
			scoreelement["x"] = std::atoi(xcoord.c_str()) - (std::atoi(boxsize.c_str())/2) + xit;
			scoreelement["y"] = std::atoi(ycoord.c_str()) - (std::atoi(boxsize.c_str())/2) + yit;
			scoreelement["score"] = score;
			scores.push_back(scoreelement);
			delete [] box->data;
			
		}
	}
	
	returnbox["x"] = "0";
	returnbox["y"] = "0";
	returnbox["size"] = boxsize;
	score = 1000000;
	
	for(i = 0; i <scores.size(); i++){
		std::cout << scores[i]["score"] << std::endl;
		if(scores[i]["score"] < score){
			returnbox["x"] = std::to_string((int) floor(scores[i]["x"]));
			returnbox["y"] = std::to_string((int) floor(scores[i]["y"]));
			score = scores[i]["score"];
		}
	}

	std::cout << returnbox["x"] << " " << returnbox["y"] << " " << returnbox["size"] << std::endl;
	response->boxes.push_back(returnbox);
	
	delete [] micrograph->data;
	delete micrograph;
	delete box;
}

void swarmLikePick(JSONResponse* response, struct http_message* message){
	
	std::string		pickrefs;
	std::string		boxfile;
	std::string		folder;
	std::string		micrographname;
	std::string		xcoord;
	std::string		ycoord;
	std::string		boxsize;
	std::string		particlediameter;
	std::string		command;
	std::string		thres;
	std::string		ndev;
	FILE*								stream;
	int				status;
	std::ofstream					filetab;
	std::ofstream					refs;
	char*			basec;
	MRCFile*		micrograph;
	MRCFile*		box;
	MRCFile*		oldbox;
	int				datalength;
	int				datacount;
	float			datamean;
	float			datadiff2;
	float			datasd;
	float			datanorm;
	
	
	getRequestVariable(message, "micrograph", micrographname);
	getRequestVariable(message, "xcoord", xcoord);
	getRequestVariable(message, "ycoord", ycoord);
	getRequestVariable(message, "boxsize", boxsize);
	getRequestVariable(message, "particlediameter", particlediameter);
	
	
	micrograph = new MRCFile();
	box = new MRCFile();
	oldbox = new MRCFile();
	
	filetab.open("/tmp/picktab.txt");
	if(filetab.is_open()){
			basec = strdup(micrographname.c_str());
				//boxfile = folder + "/";
				boxfile += std::string(basename(basec));
				boxfile.replace(boxfile.end() - 3, boxfile.end(), "box");
				if(fileExists(boxfile)){
					remove(boxfile.c_str());
				}
				filetab << micrographname << "\n";
				filetab.close();
	}
			
	readMRCFrame(micrograph, micrographname, 0);
	extractMRCBox(micrograph, box, std::atoi(xcoord.c_str()) - (std::atoi(boxsize.c_str())/2), std::atoi(ycoord.c_str()) - (std::atoi(boxsize.c_str())/2), std::atoi(boxsize.c_str()));
	
	delete [] micrograph->data;
	delete micrograph;
	
	if(!fileExists("/tmp/refs.mrc")){
		box->nx = std::stoi(boxsize);
		box->ny = std::stoi(boxsize);
		box->nz = 1;
		box->mode = 2;		
	
		memcpy(&box->header[0], &box->nx, 4);
		memcpy(&box->header[4], &box->ny, 4);
		memcpy(&box->header[8], &box->nz, 4);
		memcpy(&box->header[12], &box->mode, 4);
		
		refs.open("/tmp/refs.mrc");
		refs.write(box->header, 1024);
		refs.write((char*)box->data, box->nx * box->ny * box->nz * 4);
		refs.close();
		delete [] box->data;
		delete box;
		
	} else {
		readMRCFrame(oldbox, "/tmp/refs.mrc", 0);
		for(int i = 0; i < oldbox->nx * oldbox->ny; i++){
			oldbox->data[i] += box->data[i];
		}
	
		refs.open("/tmp/refs.mrc");
		refs.write(oldbox->header, 1024);
		refs.write((char*)oldbox->data, oldbox->nx * oldbox->ny * oldbox->nz * 4);
		refs.close();
		delete [] box->data;
		delete [] oldbox->data;
		delete oldbox;
		delete box;
	}
	
	micrograph = new MRCFile();
	box = new MRCFile();
	
	readMRCFrame(micrograph, "/tmp/refs.mrc", 0);
	readMRCFrame(box, "/tmp/refs.mrc", 0);
	
	datamean = 0;
		
	datalength = micrograph->nx * micrograph->ny;
		
	for(datacount = 0; datacount < datalength; datacount++){
		datamean += micrograph->data[datacount];
	}	

	datamean = datamean / datalength;
		
	datadiff2 = 0;

	for(datacount = 0; datacount < datalength; datacount++){
		datadiff2 += pow((micrograph->data[datacount] - datamean), 2);
	}
		
	datasd = sqrt(datadiff2 /datalength);

	for(datacount = 0; datacount < datalength; datacount++){
		datanorm = (micrograph->data[datacount] - datamean) / datasd;
		box->data[datacount] = datanorm;
	}
	
	refs.open("/tmp/refs_norm.mrc");
	refs.write(box->header, 1024);
	refs.write((char*)box->data, box->nx * box->ny * box->nz * 4);
	refs.close();
	
	delete [] box->data;
	delete [] micrograph->data;
	delete micrograph;
	delete box;
	
	command = "simple_exec prg=pick smpd=6.8 refs=/tmp/refs_norm.mrc filetab=/tmp/picktab.txt nthr=1";
		if(getRequestVariable(message, "thres", thres)){
			command += " thres="+thres;
		}
				if(getRequestVariable(message, "ndev", ndev)){
					command += " ndev="+ndev;
				}
				command += " > simple_pick.log";
				std::cout << command << std::endl;
				std::cout << boxfile << std::endl;
				stream = popen(command.c_str(), "r");
				pclose(stream);
				
				std::cout << boxfile << std::endl;
				getBoxes(response, boxfile);
			
		
	
}

void massCentre(JSONResponse* response, struct http_message* message){
	
	std::string		micrographname;
	std::string		xcoord;
	std::string		ycoord;
	std::string		boxsize;
	std::string		particlediameter;
	std::string		boxfile;
	int 			searchsize = 10;
	MRCFile*		micrograph;
	MRCFile*		box;
	float			datamean;
	float			datadiff2;
	float			datasd;
	float			datanorm;
	int				datacount;
	int				datalength;
	float 			backaverage;
	int				backgroundcount;
	float 			datasub;
	int				xit;
	int				yit;
	int 			searchxit;
	int 			searchyit;
	float 			r;
	float 			xcen;
	float			ycen;
	std::map<std::string, std::string>				returnbox;
	std::ofstream 							output;
	
	getRequestVariable(message, "micrograph", micrographname);
	getRequestVariable(message, "xcoord", xcoord);
	getRequestVariable(message, "ycoord", ycoord);
	getRequestVariable(message, "boxsize", boxsize);
	getRequestVariable(message, "particlediameter", particlediameter);
	getRequestVariable(message, "boxfile", boxfile);
	
	micrograph = new MRCFile();
	readMRCFrame(micrograph, micrographname, 0);
	
	if ((std::atoi(xcoord.c_str()) - std::atoi(boxsize.c_str())) > 0 &&
		(std::atoi(xcoord.c_str()) + std::atoi(boxsize.c_str())) < micrograph->nx &&
		(std::atoi(ycoord.c_str()) - std::atoi(boxsize.c_str())) > 0 &&
		(std::atoi(ycoord.c_str()) + std::atoi(boxsize.c_str())) < micrograph->ny){
	
		box = new MRCFile();
		
		extractMRCBox(micrograph, box, std::atoi(xcoord.c_str()) - (std::atoi(boxsize.c_str())/2), std::atoi(ycoord.c_str()) - (std::atoi(boxsize.c_str())/2), std::atoi(boxsize.c_str()));
		
		datamean = 0;
			
		datalength = box->nx * box->ny;
			
		for(datacount = 0; datacount < datalength; datacount++){
			datamean += box->data[datacount];
		}	

		datamean = datamean / datalength;
			
		datadiff2 = 0;

		for(datacount = 0; datacount < datalength; datacount++){
			datadiff2 += pow((box->data[datacount] - datamean), 2);
		}
			
		datasd = sqrt(datadiff2 /datalength);

		for(datacount = 0; datacount < datalength; datacount++){
			datanorm = (box->data[datacount] - datamean) / datasd;
			box->data[datacount] = datanorm;
		}
		
		backgroundcount = 0;
		
		for(yit = 0; yit < std::atoi(boxsize.c_str()); yit++){
			for(xit = 0; xit < std::atoi(boxsize.c_str()); xit++){
				r = sqrt(pow((xit - (std::atoi(boxsize.c_str())/2)), 2) + pow((yit - std::atoi(boxsize.c_str())/2), 2));
				if(r > (std::atoi(particlediameter.c_str()) / 2)){
					backaverage += box->data[(yit * box->nx) + xit];
					backgroundcount++;
				}
			}
		}
		
		backaverage = backaverage / backgroundcount;
		
		for(datacount = 0; datacount < datalength; datacount++){
			datasub = (box->data[datacount] - backaverage);
			box->data[datacount] = datasub;
		}
		
		backgroundcount = 0;
		xcen = 0;
		ycen = 0;
		
		for(yit = 0; yit < std::atoi(boxsize.c_str()); yit++){
			for(xit = 0; xit < std::atoi(boxsize.c_str()); xit++){
				r = sqrt(pow((xit - (std::atoi(boxsize.c_str())/2)), 2) + pow((yit - std::atoi(boxsize.c_str())/2), 2));
				if(r < (std::atoi(particlediameter.c_str()) / 2)){
					xcen += (box->data[(yit * box->nx) + xit] * (xit - (std::atoi(boxsize.c_str())/2)));
					ycen += (box->data[(yit * box->nx) + xit] * (yit - (std::atoi(boxsize.c_str())/2)));
					backgroundcount++;
				}
			}
		}
		xcen = xcen / backgroundcount;
		ycen = ycen / backgroundcount;
		
		delete [] box->data;
		delete box;
		std::cout << "centred" << std::endl;
		
	} else {
		xcen = 0;
		ycen = 0;
		std::cout << "non-centred" << std::endl;
	}
	
	delete [] micrograph->data;
	delete micrograph;
	
	output.open(boxfile.c_str(), std::ios_base::app);
	if(output.is_open()){
		output << std::to_string((int) round(std::atoi(xcoord.c_str()) - (std::atoi(boxsize.c_str())/2) - xcen)) + " " + std::to_string((int) round(std::atoi(ycoord.c_str()) - (std::atoi(boxsize.c_str())/2) - ycen)) + " " + boxsize + " " + boxsize + "\n";
		output.close();
	}
	
	getBoxes(response, boxfile);
}

void deleteBox(JSONResponse* response, struct http_message* message){
	
	std::string		xcoord;
	std::string		ycoord;
	std::string		boxsize;
	std::string		boxfilename;
	UniDoc*			boxfile;
	UniDoc*			newboxfile;
	int				i;
	int 			xmin;
	int				xmax;
	int				ymin;
	int				ymax;
	int				x;
	int				y;
	
	getRequestVariable(message, "xcoord", xcoord);
	getRequestVariable(message, "ycoord", ycoord);
	getRequestVariable(message, "boxsize", boxsize);
	getRequestVariable(message, "boxfile", boxfilename);
	
	boxfile = new UniDoc();
	newboxfile = new UniDoc();
	
	readUniDoc(boxfile, boxfilename);
	
	x = std::atoi(xcoord.c_str());
	y = std::atoi(ycoord.c_str());
	
	std::cout << x << " " << y << std::endl;
	
	for(i = 0; i < boxfile->data.size(); i++){
		getElementFromString(boxfile->data[i], xcoord, 0);
		getElementFromString(boxfile->data[i], ycoord, 1);
		xmin = std::atoi(xcoord.c_str());
		ymin = std::atoi(ycoord.c_str());
		xmax = xmin + std::atoi(boxsize.c_str());
		ymax = ymin + std::atoi(boxsize.c_str());
		std::cout << xmin << " " << ymin << std::endl;
		if(x > xmin &&
		   x < xmax &&
		   y > ymin &&
		   y < ymax){
			   std::cout << "ptcl" << std::endl;
		   } else {
			   newboxfile->data.push_back(boxfile->data[i]);
		   }
	}
	
	writeUniDoc(newboxfile, boxfilename);
	
	delete boxfile;
	delete newboxfile;
	
	getBoxes(response, boxfilename);
}

void generateGaussianRef (JPEGResponse* response, struct http_message* message) {
	
	std::string						boxsize;
	std::string						sigma;	
	std::string						directory;			
	int								xit;
	int								yit;
	float							r;
	float							gauss;
	MRCFile*						mrcfile;
	std::ofstream					output;

	if (getRequestVariable(message, "boxsize", boxsize) && getRequestVariable(message, "sigma", sigma) && getRequestVariable(message, "directory", directory)){
		
		mrcfile = new MRCFile();
		
		mrcfile->nx = std::stoi(boxsize);
		mrcfile->ny = std::stoi(boxsize);
		mrcfile->nz = 1;
		mrcfile->mode = 2;
		mrcfile->data = new float[mrcfile->nx * mrcfile->ny];
		
		for(yit = 0; yit < mrcfile->ny; yit++){
			for(xit = 0; xit < mrcfile->nx; xit++){ 
				r = sqrt(pow((xit - (mrcfile->nx/2)), 2) + pow((yit - (mrcfile->ny/2)), 2));
				gauss = exp(-((pow(r, 2))/(2 * pow((std::stoi(sigma) / 2), 2))));
				//mrcfile->data[(yit * mrcfile->nx) + xit] = (0.5 - gauss); // This Works
				mrcfile->data[(yit * mrcfile->nx) + xit] = (-0.1 - gauss);// No zeros ;-)
			}
		}
		
		memcpy(&mrcfile->header[0], &mrcfile->nx, 4);
		memcpy(&mrcfile->header[4], &mrcfile->ny, 4);
		memcpy(&mrcfile->header[8], &mrcfile->nz, 4);
		memcpy(&mrcfile->header[12], &mrcfile->mode, 4);

		output.open(directory + "/gausspick.mrc",std::ifstream::binary);
		output.write(mrcfile->header, 1024);
		output.write((char*)mrcfile->data, mrcfile->nx * mrcfile->ny * mrcfile->nz * 4);
		output.close();
		
		getPixelsFromMRC(response, directory + "/gausspick.mrc", message);
			
		delete [] mrcfile->data;
		delete mrcfile;
	}
}

void getDirectoryContents (JSONResponse* response, struct http_message* message) {
	
	std::string						directory;
	std::string						filter;
	int								fileit;
	
	if(!getRequestVariable(message, "directoryname", directory)) {
		directory = "/";
	}
	
	
	if(fileExists(directory)){
		listDirectory(response->files, response->directories, directory);
		response->rootdirectory = directory;
	}
	
	if(getRequestVariable(message, "filefilter", filter)) {
		for(fileit = response->files.size() - 1; fileit >= 0; fileit--){
			if(response->files[fileit].find(filter) == std::string::npos){
				response->files.erase(response->files.begin() + fileit);
			}
		}

	}

}

void encodeJSON(JSONResponse* response){
	
	int 											i;
	int												j;
	std::map<std::string,std::string>::iterator 	mapit;
	
	response->JSONstring = "{ ";
	
	if (response->projects.size() > 0) {
		response->JSONstring += "\"projects\" : [ ";
		for (i = 0; i < response->projects.size(); i++){
			response->JSONstring += "{ ";
			for (mapit=response->projects[i].begin(); mapit!=response->projects[i].end(); ++mapit){
				response->JSONstring += "\"" + mapit->first + "\" : \"" + mapit->second + "\",";
			}
			response->JSONstring.pop_back();
			response->JSONstring += "},";
		}
		response->JSONstring.pop_back();
		response->JSONstring += "],";
	}
	
	if (response->snapshots.size() > 0) {
		response->JSONstring += "\"snapshots\" : [ ";
		for (i = 0; i < response->snapshots.size(); i++){
			response->JSONstring += "{ ";
			for (mapit=response->snapshots[i].begin(); mapit!=response->snapshots[i].end(); ++mapit){
				response->JSONstring += "\"" + mapit->first + "\" : \"" + mapit->second + "\",";
			}
			response->JSONstring.pop_back();
			response->JSONstring += "},";
		}
		response->JSONstring.pop_back();
		response->JSONstring += "],";
	}
		
	if (response->jobs.size() > 0) {
		response->JSONstring += "\"jobs\" : [ ";
		for (i = 0; i < response->jobs.size(); i++){
			response->JSONstring += "{ ";
			for (mapit=response->jobs[i].begin(); mapit!=response->jobs[i].end(); ++mapit){
				response->JSONstring += "\"" + mapit->first + "\" : \"" + mapit->second + "\",";
			}
			response->JSONstring.pop_back();
			response->JSONstring += "},";
		}
		response->JSONstring.pop_back();
		response->JSONstring += "],";
	}
	
	if (response->iterations.size() > 0) {
		response->JSONstring += "\"iterations\" : [ ";
		for (i = 0; i < response->iterations.size(); i++){
			response->JSONstring += "{ ";
			for (mapit=response->iterations[i].begin(); mapit!=response->iterations[i].end(); ++mapit){
				response->JSONstring += "\"" + mapit->first + "\" : \"" + mapit->second + "\",";
			}
			response->JSONstring.pop_back();
			response->JSONstring += "},";
		}
		response->JSONstring.pop_back();
		response->JSONstring += "],";
	}
	
	if (response->boxes.size() > 0) {
		response->JSONstring += "\"boxes\" : [ ";
		for (i = 0; i < response->boxes.size(); i++){
			response->JSONstring += "[ ";
			response->JSONstring += "\"" + response->boxes[i]["x"] + "\", \"" + response->boxes[i]["y"] + "\", \"" + response->boxes[i]["size"] + "\"";
			response->JSONstring += "],";
		}
		response->JSONstring.pop_back();
		response->JSONstring += "],";
	}
	
	if (response->files.size() > 0) {
		response->JSONstring += "\"files\" : [ ";
		for (i = 0; i < response->files.size(); i++){
			response->JSONstring += "\"" + response->files[i] + "\",";
		}
		response->JSONstring.pop_back();
		response->JSONstring += "],";
	}
	
	if (response->directories.size() > 0) {
		response->JSONstring += "\"directories\" : [ ";
		for (i = 0; i < response->directories.size(); i++){
			response->JSONstring += "\"" + response->directories[i] + "\",";
		}
		response->JSONstring.pop_back();
		response->JSONstring += "],";
	}
	
	if (response->rootdirectory.size() > 0) {
		response->JSONstring += "\"rootdirectory\" : \"" + response->rootdirectory + "\",";
	}
	
	if (response->prime2ddoc.size() > 0) {
		response->JSONstring += "\"prime2ddoc\" : \"" + response->prime2ddoc + "\",";
	}
	
	if (response->inputfilename.size() > 0) {
		response->JSONstring += "\"inputfilename\" : \"" + response->inputfilename + "\",";
	}
	
	if (response->logfile.size() > 0) {
		response->JSONstring += "\"logfile\" : \"" + response->logfile + "\",";
	}
	
	if (response->error.size() > 0) {
		response->JSONstring += "\"error\" : \"" + response->error + "\",";
	}
	
	if (response->particlecount.size() > 0) {
		response->JSONstring += "\"particlecount\" : \"" + response->particlecount + "\",";
	}
	
	if (response->jobfolder.size() > 0) {
		response->JSONstring += "\"jobfolder\" : \"" + response->jobfolder + "\",";
	}
	
	response->JSONstring.pop_back();
	response->JSONstring += "}";
}

void JSONHandler (struct mg_connection* http_connection, struct http_message* message) {
	
	std::string				function;
	JSONResponse*			response;
	
	response = new JSONResponse();
	
	if (getRequestVariable(message, "function", function)) {
		std::cout << function << std::endl;
		if (function == "getprojects") {
			getProjects(response);
		} else if (function == "listdir") {
			getDirectoryContents(response, message);
		} else if (function == "newproject") {
			newProject(response, message);
		} else if (function == "getjobs") {
			getJobs(response, message);
		} else if (function == "externaljob") {
			externalJob(response, message);
		} else if (function == "pipelineview") {
			viewPreproc(response, message);
		} else if (function == "boxfiledata") {
			getBoxes(response, message);
		} else if (function == "savesel") {
			savePreprocSelection(response, message);
		} else if (function == "preprocsaveselectionparticles") {
			savePreprocSelectionParticles(response, message);
		} else if (function == "simplejob") {
			simpleJob(response, message);
		} else if (function == "2dview") {
			viewPrime2D(response, message);
		} else if (function == "2dviewiteration") {
			viewPrime2DIteration(response, message);
		} else if (function == "2dviewparticles") {
			viewPrime2DParticles(response, message);
		} else if (function == "2dsaveselection") {
			savePrime2DSelection(response, message);
		} else if (function == "2dsaveselectionparticles") {
			savePrime2DSelectionParticles(response, message);
		} else if (function == "viewlogfile") {
			getLogFile(response, message);
		} else if (function == "refinebox") {
			massCentre(response, message);
		} else if (function == "deletebox") {
			deleteBox(response, message);
		} else if (function == "unblurview") {
			viewUnblur(response, message);
		} else if (function == "saveunblurselection") {
			saveUnblurSelection(response, message);
		} else if (function == "manualpick") {
			viewManualPick(response, message);
		} else if (function == "deletejob") {
			deleteJob(response, message);
		//} else if (function == "ctffindview") {
		//	viewCtffind(response, message);
		} else if (function == "deleteproject") {
			deleteProject(response, message);
		} else if (function == "ctffindview") {
			viewCtffit(response, message);
		} else if (function == "extractview") {
			viewExtract(response, message);
		} else if (function == "syncjob") {
			syncJob(response, message);
		} else if (function == "clearboxes") {
			clearBoxes(response, message);
			getBoxes(response, message);
		} else if (function == "autopick") {
			autoPick(response, message);
		} else if (function == "teststring") {
			std::string inverseselection;
			getRequestVariable(message, "inverseselection", inverseselection);
			std::cout << inverseselection << std::endl;
		}else if (function == "ini3dpost") {
			std::string dir;
			getRequestVariable(message, "directory", dir);
			std::cout << dir << std::endl;
			ini3DPost(dir);
		}
	}
	
	encodeJSON(response);

	mg_send_head(http_connection, 200, response->JSONstring.length(), "Content-Type: application/json");
	mg_send(http_connection, response->JSONstring.c_str(), response->JSONstring.length());
	
	delete response;
	
}

void encodeJPEG(JPEGResponse* response, int size) {
	
	struct jpeg_compress_struct 		cinfo;
	struct jpeg_error_mgr 				jerr;
	int									numerator;
	int 								newsize;
	
	for(numerator = 1; numerator <= 8; numerator++){
		newsize = response->xdim * numerator / 8;
		if(newsize > size){
			std::cout << size << " " << newsize << " " << numerator << std::endl;
			break;
		}
	}
	
	if(response->xdim > 0){
		cinfo.err = jpeg_std_error(&jerr);
		jpeg_create_compress(&cinfo);
		jpeg_mem_dest(&cinfo, &response->jpeg, &response->jpegsize);
		cinfo.image_width = response->xdim;
		cinfo.image_height = response->ydim;
		cinfo.input_components = 1;
		cinfo.in_color_space = JCS_GRAYSCALE;
		jpeg_set_defaults(&cinfo);
		cinfo.scale_num = numerator;
		cinfo.scale_denom = 8;
		jpeg_set_quality(&cinfo, 50, TRUE);
		jpeg_start_compress(&cinfo, TRUE);
		JSAMPROW row_pointer;
		for(int i = 0; i < response->ydim; i++){
			row_pointer = (JSAMPROW) &response->pixels[i * response->xdim];
			jpeg_write_scanlines(&cinfo, &row_pointer, 1);
		}
		jpeg_finish_compress(&cinfo);
		jpeg_destroy_compress(&cinfo);
	}
}

void encodeJPEG(JPEGResponse* response) {
	
	struct jpeg_compress_struct 		cinfo;
	struct jpeg_error_mgr 				jerr;
	std::cout << sizeof(response->pixels) << std::endl;
	if(response->xdim > 0){
		cinfo.err = jpeg_std_error(&jerr);
		jpeg_create_compress(&cinfo);
		jpeg_mem_dest(&cinfo, &response->jpeg, &response->jpegsize);
		cinfo.image_width = response->xdim;
		cinfo.image_height = response->ydim;
		cinfo.input_components = 1;
		cinfo.in_color_space = JCS_GRAYSCALE;
		jpeg_set_defaults(&cinfo);
		jpeg_set_quality(&cinfo, 50, TRUE);
		jpeg_start_compress(&cinfo, TRUE);
		JSAMPROW row_pointer;
		for(int i = 0; i < response->ydim; i++){
			row_pointer = (JSAMPROW) &response->pixels[i * response->xdim];
			jpeg_write_scanlines(&cinfo, &row_pointer, 1);
		}
		jpeg_finish_compress(&cinfo);
		jpeg_destroy_compress(&cinfo);
	}
}

void JPEGHandler (struct mg_connection* http_connection, struct http_message* message) {
 
	std::string			filename;
	std::string			size;
	std::string			function;
	JPEGResponse*		response;
	
	response = new JPEGResponse();
	
	if (getRequestVariable(message, "function", function)) {
		if (function == "createpickref") {
			generateGaussianRef(response, message);
		}
	}
	
	if (getRequestVariable(message, "filename", filename)) {
		if (filename.find(".mrc") != std::string::npos && fileExists(filename)){
			getPixelsFromMRC(response, filename, message); 
		}
	}

	if(getRequestVariable(message, "size", size)) {
		encodeJPEG(response, std::stoi(size));
	}else{
		encodeJPEG(response);
	}
	
	mg_send_head(http_connection, 200, response->jpegsize, "Content-Type: image/jpeg");
	mg_send(http_connection, response->jpeg, response->jpegsize);
	delete [] response->pixels;
	delete [] response->jpeg;
	delete response;
}

void webEventHandler (struct mg_connection* http_connection, int event, void *eventdata) {
	
	struct http_message*   					message;
	
	message = (struct http_message *) eventdata;
	
	switch (event) {
		case MG_EV_HTTP_REQUEST: {
			if (mg_vcmp(&message->uri, "/JSONhandler") == 0) {
				JSONHandler(http_connection, message);
			} else if (mg_vcmp(&message->uri, "/JPEGhandler") == 0) {
				JPEGHandler(http_connection, message);
			} else {
				mg_serve_http(http_connection, message, http_server_opts);
			}
		}
		default: {
			break;
		}
	}
}

void webServer(){
	
	struct mg_mgr				http_manager;
	struct mg_connection*		http_connection;
	sig_atomic_t				http_signal_received;
	
	mg_mgr_init(&http_manager, NULL);
	
	std::cout << "Starting webserver on port 8088" << " ... ";
	
	http_connection = mg_bind(&http_manager, "8088", webEventHandler);
	
	if (http_connection == NULL) {
		std::cout << "Failed" << std::endl;
		exit(1);
	}else{
		std::cout << "Success" << std::endl;
	}
	
	mg_set_protocol_http_websocket(http_connection);	

	http_server_opts.document_root = ".";
	
	http_server_opts.enable_directory_listing = "no";
	
	http_signal_received = 0;
	
	while (http_signal_received == 0) {
		mg_mgr_poll(&http_manager, 200);
	}

	mg_mgr_free(&http_manager);
}

int main(void) {
    
    std::string     wwwdirectory;
    char*           simplepath;
    
    simplepath = getenv("SIMPLE_PATH");
    
    if(simplepath == NULL){
        std::cout << "SIMPLE_PATH is not set" << std::endl;
        return 1;
    }
    
    wwwdirectory = std::string(simplepath) + "/www";
    
    if(fileExists(wwwdirectory)){
        if(chdir(wwwdirectory.c_str()) == 0){
            webServer();
        }else{
            std::cout << "Error: Could not change directory to www directory" << std::endl;
            return 1;
        }
    }else {
        std::cout << "Error: Directory set in SIMPLE_PATH does not exist" << std::endl;
        
    }
	return 0;
}
