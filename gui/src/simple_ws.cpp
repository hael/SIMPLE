#define KEY(X) #X
#define VALUE(X) KEY(X)
#define SIMPLE_DIR_VALUE VALUE(SIMPLE_DIR)

#include "base64.h"
#include "lodepng.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <pwd.h>
#include <unistd.h>
#include <sys/stat.h>
#include <dirent.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstring>
#include <cmath>
#include <signal.h>

void missingArgumentError(std::string arg){ 											// Print error and exit
	std::cout << "status=error message=" << arg << " not declared" << std::endl;
	exit(1);
}

void fatalError(std::string message){ 													// Print error and exit
	std::cout << "status=error message=" << message << std::endl;
	exit(1);
}

void completionSuccess(){ 																// Print success and exit
	std::cout << "status=complete" << std::endl;
	exit(0);
}	
	
void readMRC(std::string mrcfilename, int *nx, int *ny, int *nz, int *mode, std::vector<float> *data){
	std::ifstream mrcfile;
	int pixcount;
	float pixvalue;
	
	mrcfile.open(mrcfilename.c_str(), std::ifstream::binary);
	
	if (mrcfile.is_open()){
		mrcfile.read ((char*)nx, sizeof(int));
		mrcfile.read ((char*)ny, sizeof(int));
		mrcfile.read ((char*)nz, sizeof(int));
		mrcfile.read ((char*)mode, sizeof(int));
		mrcfile.seekg(1024, std::fstream::beg);
		if(*mode == 2){
			pixcount = 0;
			while(pixcount < (*nx * *ny)){
				mrcfile.read((char*)&pixvalue, sizeof(float));
				data->push_back(pixvalue);
				pixcount++;
			}
		}
	}
	mrcfile.close();
}

void readMRCHeader(std::string mrcfilename, int *nx, int *ny, int *nz, int *mode){
	std::ifstream mrcfile;
	
	mrcfile.open(mrcfilename.c_str(), std::ifstream::binary);
	
	if (mrcfile.is_open()){
		mrcfile.read ((char*)nx, sizeof(int));
		mrcfile.read ((char*)ny, sizeof(int));
		mrcfile.read ((char*)nz, sizeof(int));
		mrcfile.read ((char*)mode, sizeof(int));
	}
	mrcfile.close();
}

std::string getKeyValFromString(std::string argstring, std::string key){ 				// Return string containing value from space separated keyval pairs in input string
	std::size_t location;
	char character;
	std::string returnstring;
	
	returnstring = "";
	key.append("=");
	location = argstring.find(key) + key.length();
	if(location != key.length() - 1){
		while ( location < argstring.length() ){
			character = argstring.c_str()[location];
			if ( character != ' '){
				returnstring.push_back(character);
			} else {
				break;
			}
			location = location + 1;
		}
	}
	return returnstring;
}

std::string getElementFromString(std::string argstring, int element){ 					// Return string containing value from space separated keyval pairs in input string
	std::size_t location;
	char character;
	std::string returnstring;
	int currentelement;
	
	returnstring = "";
	location = 0;
	currentelement = 0;
	
	while(argstring.c_str()[location] == ' ' && location < argstring.length()){
		location = location + 1;
	}
	
	while ( location < argstring.length() ){
		character = argstring.c_str()[location];
		if ( character != ' '){
			if(currentelement == element){
				returnstring.push_back(character);
			}
		} else {
			currentelement = currentelement + 1;
		}
		location = location + 1;
	}
	return returnstring;
}

std::string getUser(std::string argstring){ 											// Return user or daemon user if local
	std::string user;
	struct passwd* pw;
	register uid_t uid;
	
	user = getKeyValFromString(argstring, "usr");
	
	if(user == "" || user == "local"){
		uid = geteuid();
		pw = getpwuid(uid);
		user = pw->pw_name;
	}
	
	return user;
}

std::string getUserHome(std::string user){ 												// Return home directory of user
	std::string userhome;
	struct passwd* pw;
	
	pw = getpwnam(user.c_str());
	userhome = pw->pw_dir;
	
	return userhome;
}

std::vector<std::string> readProjectConf(std::string argstring){						// Return array with contents of project configuration file
	std::string user, userproject, line;
	std::ifstream projectfile;
	std::vector<std::string> returnvector;
	struct stat buffer;
	
	user = getUser(argstring);
	userproject = getUserHome(user);
	userproject.append("/.simple/projects.simple");
	
	if(stat(userproject.c_str(), &buffer) == 0){
		projectfile.open(userproject.c_str());
		while (std::getline(projectfile, line)){
			returnvector.push_back(line);
		}
		projectfile.close();
	}else{
		//CREATE PROJECT FILE
	}
	
	return returnvector;
}

std::vector<std::string> readBoxes(std::string boxfile){
	struct stat buffer;
	std::string line, returnstring;
	std::vector<std::string> returnarray;
	std::ifstream openboxfile;
	bool nextspace = false;
	
	if(stat(boxfile.c_str(), &buffer) == 0){
		openboxfile.open(boxfile.c_str());
		if (openboxfile.is_open()){
			while(getline(openboxfile,line)){
				nextspace = false;
				returnstring = "";
				for(unsigned count=0; count < line.length(); count++) {
					if(line[count] != ' ') {
						nextspace = true;
						returnstring.append(1, line[count]);
					} else {
						if(nextspace){
							returnstring.append(",");
						}
						nextspace = false;
					}
				}
				returnarray.push_back(returnstring);
				
			}
			openboxfile.close();
		}
	}	
	return returnarray;
}

std::string getProjectFolder(std::string argstring){									// Return project directory for project
	std::string projectname, projecttest, returnpath;
	std::vector<std::string> projectconf;
	
	returnpath = "";
	projectname = getKeyValFromString(argstring, "projectname");
	
	if(projectname == ""){
		missingArgumentError("projectname");
	}
	
	projectconf = readProjectConf(argstring);
	
	for (std::vector<std::string>::iterator it=projectconf.begin(); it!=projectconf.end(); ++it){
		projecttest = getKeyValFromString(*it, "projectname");
		if(projectname == projecttest){
			returnpath = getKeyValFromString(*it, "projectdir");
			break;
		}
	}
	return returnpath;
}

std::vector<std::string> readJobConf(std::string argstring){							// Return array with contents of job configuration file
	std::string projectjob, line;
	std::ifstream jobfile;
	std::vector<std::string> returnvector;
	struct stat buffer;
	
	projectjob = getProjectFolder(argstring);
	projectjob.append("/.simple/jobs.simple");
	
	if(stat(projectjob.c_str(), &buffer) == 0){
		jobfile.open(projectjob.c_str());
		while (std::getline(jobfile, line)){
			returnvector.push_back(line);
		}
		jobfile.close();
	}else{
		//CREATE JOB FILE
	}
	
	return returnvector;
}

std::string getJobFolder(std::string argstring, int jobid){									// Return project directory for project
	std::string jobtest, returnpath;
	std::vector<std::string> jobconf;
	std::stringstream ss;
	
	returnpath = "";
	ss << jobid;
	
	jobconf = readJobConf(argstring);
	
	for (std::vector<std::string>::iterator it=jobconf.begin(); it!=jobconf.end(); ++it){
		jobtest = getKeyValFromString(*it, "jobid");
		if(ss.str() == jobtest){
			returnpath = getKeyValFromString(*it, "jobdir");
			break;
		}
	}
	return returnpath;
}

std::string readCTFParams(std::string ctfparamsfile){									// Return contents of CTF parameters file
	struct stat buffer;
	std::string returnstring;
	std::ifstream ctffile;
	
	returnstring = "";
	
	if(stat(ctfparamsfile.c_str(), &buffer) == 0){
		ctffile.open(ctfparamsfile.c_str());
		if (ctffile.is_open()){
			getline(ctffile,returnstring);
			ctffile.close();
		}
	}
	return returnstring;
}

std::string data2PNGString(std::vector<float>* data, int nx, int ny, float brightness , float contrast){
	float datanorm;
	float datatotal;
	float datamean;
	float datadiff2;
	float datasd;
	std::vector<unsigned char> image;
	std::vector<unsigned char> png;
	std::string pngstring, returnstring;
	datatotal = 0.0;

	for (std::vector<float>::iterator it=data->begin(); it!=data->end(); ++it){
		datatotal += *it;
	}
	datamean = datatotal / data->size();
	
	datadiff2 = 0.0;
	for (std::vector<float>::iterator it=data->begin(); it!=data->end(); ++it){
		datadiff2 += pow((*it - datamean),2);
	}
	datasd = sqrt(datadiff2 / data->size());
	
	for (std::vector<float>::iterator it=data->begin(); it!=data->end(); ++it){
		datanorm = (*it - datamean) / datasd;
		datanorm = datanorm * (255 / contrast);
		datanorm = datanorm + brightness;
		if(datanorm > 255){
			datanorm = 255;
		} else if(datanorm < 0){
			datanorm = 0;
		}
		image.push_back(datanorm);
		image.push_back(datanorm);
		image.push_back(datanorm);
		image.push_back(255);
	}
	
	lodepng::encode(png, image, nx, ny);
	
	for (std::vector<unsigned char>::iterator it=png.begin(); it!=png.end(); ++it){
		pngstring.push_back(*it);
	}
	
	returnstring = base64_encode(reinterpret_cast<const unsigned char*>(pngstring.c_str()), pngstring.length()); 
	
	return returnstring;
}

std::string pngFromMRC(std::string mrcfile, int brightness, int contrast){											// Return base64 encoded png string  
	int nx, ny, nz, mode;
	std::vector<float> data;
	std::string returnstring;
	struct stat buffer;
	
	returnstring = "";
	if(stat(mrcfile.c_str(), &buffer) == 0){
		readMRC(mrcfile, &nx, &ny, &nz, &mode, &data);
		returnstring = data2PNGString(&data, nx, ny, brightness, contrast);
	}
	return returnstring;
}

void pngsFromMRCStack(std::string mrcfilename, int brightness, int contrast){											// Return base64 encoded png string  
	int nx, ny, nz, mode, pixcount, zcount;
	std::vector<float> data;
	std::string returnstring;
	struct stat buffer;
	float pixvalue;
	std::ifstream mrcfile;
	
	if(stat(mrcfilename.c_str(), &buffer) == 0){
		mrcfile.open(mrcfilename.c_str(), std::ifstream::binary);
		if (mrcfile.is_open()){
			mrcfile.read((char*)&nx, sizeof(int));
			mrcfile.read((char*)&ny, sizeof(int));
			mrcfile.read((char*)&nz, sizeof(int));
			mrcfile.read((char*)&mode, sizeof(int));
			mrcfile.seekg(1024, std::fstream::beg);
			zcount = 0;
			if(mode == 2){
				while(zcount < nz){
					pixcount = 0;
					data.clear();
					while(pixcount < (nx * ny)){
						mrcfile.read((char*)&pixvalue, sizeof(float));
						data.push_back(pixvalue);
						pixcount++;
					}
					returnstring = data2PNGString(&data, nx, ny, brightness, contrast);
					std::cout << "image=" << returnstring << " nz=" <<  zcount << " nx=" << nx << " ny=" << ny << std::endl;
					zcount++;
				}
			}
		}
		mrcfile.close();
	}

}

int addJob(std::string argstring){												// Add job to job configuration file
	std::string jobdir, jobtype, jobname, projectdir, line;
	std::fstream jobfile;
	std::ostringstream jobid;
	int jobcount;
	struct stat buffer;
	
	jobdir = getKeyValFromString(argstring, "jobdir");

	projectdir = getProjectFolder(argstring);
	projectdir.append("/.simple");

	jobcount = 1;
	
	if(stat(projectdir.c_str(), &buffer) != 0){
		mkdir(projectdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}

	projectdir.append("/jobs.simple");
	jobfile.open(projectdir.c_str(), std::ios::in);
	
	if (jobfile.is_open()){
		while (std::getline(jobfile, line)){
			jobcount++;
		}
		jobfile.close();
	}
	
	if(jobdir == ""){
		jobid << jobcount;
		jobtype = getKeyValFromString(argstring, "jobtype");
		jobdir = getProjectFolder(argstring);
		jobdir.append("/");
		jobdir.append(jobid.str());
		jobdir.append("_");
		jobdir.append(jobtype);
	}	
	
	jobfile.open(projectdir.c_str(), std::ios::app);
	
	if (jobfile.is_open()){
		jobfile << "jobid=" << jobcount << " jobdir=" << jobdir << " jobstatus=Running" << " jobpid=" << " " << argstring << "\n";
		jobfile.close();
	}else{
		fatalError("Failed to open jobfile");
	}
	
	return jobcount;
}

void editJob(std::string argstring, int jobid, std::string key, std::string value){
	std::string projectjob, projectjobtmp, jobline;
	std::vector<std::string> jobconf;
	int random;
	std::stringstream ssrandom, ssjobid;
	size_t pos, pos2;
	char character;
	std::ofstream tmpfile;
	
	srand (time(NULL));
	random = rand() % 100000 + 1000;

	ssrandom << random;
	ssjobid << jobid;
	
	jobconf = readJobConf(argstring);
	
	projectjob = getProjectFolder(argstring);
	projectjob.append("/.simple/jobs.simple");
	
	projectjobtmp = getProjectFolder(argstring);
	projectjobtmp.append("/.simple/jobs.simple." + ssrandom.str() + ".tmp");
	
	tmpfile.open(projectjobtmp.c_str());
	
	for (std::vector<std::string>::iterator it=jobconf.begin(); it!=jobconf.end(); ++it){
		jobline = *it;
		if(getKeyValFromString(jobline, "jobid") == ssjobid.str()){
			key.append("=");
			pos = jobline.find(key);
			if(pos != key.length() - 1 && pos != std::string::npos){
				pos2 = pos;
				character = jobline[pos];
				while( character != ' ' && pos2 != std::string::npos){
					pos2++;
					character = jobline[pos2];
				}
				jobline.erase(pos, (pos2 - pos + 1));
			}
			jobline.append(" ");
			jobline.append(key);
			jobline.append(value);
		}
		tmpfile << jobline << std::endl;
	}
	tmpfile.close();
	rename(projectjobtmp.c_str(), projectjob.c_str());
}
	
void killJob(std::string argstring){
	std::vector<std::string> jobconf;
	std::stringstream ssjobid;
	std::string pid, jobidstr;
	int jobid;
	
	jobconf = readJobConf(argstring);
	jobidstr = getKeyValFromString(argstring, "jobid");
	jobid = std::atoi(jobidstr.c_str());
	
	ssjobid << jobid;
	
	for (std::vector<std::string>::iterator it=jobconf.begin(); it!=jobconf.end(); ++it){
		if(getKeyValFromString(*it, "jobid") == ssjobid.str()){
			pid = getKeyValFromString(*it, "jobpid");
			kill(std::atoi(pid.c_str()), 9);
			editJob(argstring, jobid, "jobstatus", "Cancelled");
		}
	}	
}

void createSimpleEnv(std::string jobdir, std::string argstring){
	std::ofstream envfile;
	
	jobdir.append("/simple_distr_config.env");
	envfile.open(jobdir.c_str());
	envfile << "simple_path = " << SIMPLE_DIR_VALUE << "\n";
	
	if(getKeyValFromString(argstring, "time_per_image") != ""){
		envfile << "time_per_image = " << getKeyValFromString(argstring, "time_per_image") << "\n";
	}
	
	if(getKeyValFromString(argstring, "qsys_qos") != ""){
		envfile << "qsys_qos = " << getKeyValFromString(argstring, "qsys_qos") << "\n";
	}
	
	if(getKeyValFromString(argstring, "qsys_name") != ""){
		envfile << "qsys_name = " << getKeyValFromString(argstring, "qsys_name") << "\n";
	}
	
	if(getKeyValFromString(argstring, "qsys_partition") != ""){
		envfile << "qsys_partition = " << getKeyValFromString(argstring, "qsys_partition") << "\n";
	}
	
	if(getKeyValFromString(argstring, "job_ntasks") != ""){
		envfile << "job_ntasks = " << getKeyValFromString(argstring, "job_ntasks") << "\n";
	}
	
	if(getKeyValFromString(argstring, "job_memory_per_task") != ""){
		envfile << "job_memory_per_task = " << getKeyValFromString(argstring, "job_memory_per_task") << "\n";
	}
	
	if(getKeyValFromString(argstring, "job_ntasks_per_socket") != ""){
		envfile << "job_ntasks_per_socket = " << getKeyValFromString(argstring, "job_ntasks_per_socket") << "\n";
	}
	
	if(getKeyValFromString(argstring, "job_cpus_per_task") != ""){
		envfile << "job_cpus_per_task = " << getKeyValFromString(argstring, "job_cpus_per_task") << "\n";
	}
	envfile.close();
}

void deleteJob(std::string argstring){
	int jobid;
	std::string jobidstr;
	
	jobidstr = getKeyValFromString(argstring, "jobid");
	jobid = std::atoi(jobidstr.c_str());
	
	editJob(argstring, jobid, "deleted", "yes");
	// DELETE FOLDER HERE!!!!!!
}

void listDirectory(std::string argstring){ 												// List contents of directory if given, else users home
	std::string user, directory;
	struct stat buffer;
	DIR *directorycontents;
	struct dirent *entry;
	std::vector<std::string> directories, files;
	
	directory = getKeyValFromString(argstring, "dir");
	
	if(directory == ""){
		user = getUser(argstring);
		directory = getUserHome(user);
	}
	
	if(stat(directory.c_str(), &buffer) == 0){
		std::cout << "parentdirectory=" << directory << std::endl;
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
		std::sort (directories.begin(), directories.end());
		for (std::vector<std::string>::iterator it=directories.begin(); it!=directories.end(); ++it){
			std::cout << "directory=" << directory << "/" << *it << std::endl;
		}
		std::sort (files.begin(), files.end());
		for (std::vector<std::string>::iterator it=files.begin(); it!=files.end(); ++it){
			std::cout << "file=" << directory << "/" << *it << std::endl;
		}
	}else{
		fatalError("Directory doesn't exist");
	}
}

void listProjects(std::string argstring){												// List projects from users project configuration file
	std::vector<std::string> projectconf;
	
	projectconf = readProjectConf(argstring);
	
	for (std::vector<std::string>::iterator it=projectconf.begin(); it!=projectconf.end(); ++it){
		std::cout << *it << std::endl;
	}
}

void addProject(std::string argstring){													// Add project to users project configuration file
	std::string projectname, projectdir, projecttest, user, userproject;
	std::vector<std::string> projectconf;
	std::ofstream projectfile;
	
	projectname = getKeyValFromString(argstring, "projectname");
	projectdir = getKeyValFromString(argstring, "projectdir");
	
	if(projectname == ""){
		missingArgumentError("projectname");
	}
	
	if(projectdir == ""){
		missingArgumentError("projectdir");
	}
	
	user = getUser(argstring);
	userproject = getUserHome(user);
	userproject.append("/.simple/projects.simple");

	projectconf = readProjectConf(argstring);
	
	for (std::vector<std::string>::iterator it=projectconf.begin(); it!=projectconf.end(); ++it){
		projecttest = getKeyValFromString(*it, "projectname");
		if(projecttest == projectname){
			fatalError("A project with this name already exists");
		}
	}
	
	projectfile.open(userproject.c_str(), std::ios::app);
	projectfile << "projectname=" << projectname << " " << "projectdir=" << projectdir << " " << "\n";
	projectfile.close();
	
	projectconf = readProjectConf(argstring);
	
	for (std::vector<std::string>::iterator it=projectconf.begin(); it!=projectconf.end(); ++it){
		projecttest = getKeyValFromString(*it, "projectname");
		if(projecttest == projectname){
			completionSuccess();
		}
	}
	fatalError("Failed to create project");
}

void listJobs(std::string argstring){													// List jobs from belonging to project
	std::vector<std::string> jobconf;
	
	jobconf = readJobConf(argstring);
	
	for (std::vector<std::string>::iterator it=jobconf.begin(); it!=jobconf.end(); ++it){
		if(getKeyValFromString(*it, "deleted") != "yes"){
			std::cout << *it << std::endl;
		}
	}
}

void listJob(std::string argstring){													// List jobs from belonging to project
	std::vector<std::string> jobconf;
	std::string jobid;
	
	jobconf = readJobConf(argstring);
	jobid = getKeyValFromString(argstring, "jobid");
	
	for (std::vector<std::string>::iterator it=jobconf.begin(); it!=jobconf.end(); ++it){
		if(getKeyValFromString(*it, "jobid") == jobid){
			std::cout << *it << std::endl;
			break;
		}
	}
}

void addExternalJob(std::string argstring){
	std::string jobdir;
	int jobid;
	
	jobid = addJob(argstring);
	editJob(argstring, jobid, "jobstatus", "External");
}

void simpleDistrExec(std::string argstring){
	FILE* pipe;
	int pid, jobid, returnstatus;
	char buffer[512];
	std::stringstream ss, sspid;
	std::string argscommand, line, command, arg, jobdir, path, simplecommand;
	std::vector<std::string> commandargs;
	struct stat dirbuffer;
	
	jobid = addJob(argstring);
	jobdir = getJobFolder(argstring, jobid);
	
	if(stat(jobdir.c_str(), &dirbuffer) != 0){
		mkdir(jobdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}

	pid = fork();
	
	if(pid == 0){
		path = std::getenv("PATH");
		path.append(":");
		path.append(SIMPLE_DIR_VALUE);
		path.append("/bin");
		setenv("PATH", path.c_str(), 1); 
		
		argscommand = "simple_distr_exec prg=";
		argscommand.append(getKeyValFromString(argstring, "prg"));
	
		pipe = popen(argscommand.c_str(), "r");
		while (!feof(pipe)) {
			if (fgets(buffer, sizeof(buffer), pipe) != NULL){
				ss << buffer;
			}
		}
	
		pclose(pipe);
	
		while(std::getline(ss,line,'\n')){
			commandargs.push_back(getElementFromString(line, 0));
		}
	
		createSimpleEnv(jobdir, argstring);
	
		simplecommand = "simple_distr_exec prg=";
		simplecommand.append(getKeyValFromString(argstring, "prg"));
	
		for (std::vector<std::string>::iterator it=commandargs.begin(); it!=commandargs.end(); ++it){
			arg = getKeyValFromString(argstring, *it);
			if(arg != "" && *it != ""){
				simplecommand.append(" ");
				simplecommand.append(*it);
				simplecommand.append("=");
				simplecommand.append(arg);
			}
		}
		
		simplecommand.append(" 2>&1 > job.log");
		
		chdir(jobdir.c_str());
		pipe = popen(simplecommand.c_str(), "r");
		returnstatus = pclose(pipe);
		if(WEXITSTATUS(returnstatus) == 0){
			editJob(argstring, jobid, "jobstatus", "Complete");
		}else{
			editJob(argstring, jobid, "jobstatus", "Failed");
		}
		
	}else if(pid > 0){
		sspid << pid;
		editJob(argstring, jobid, "jobpid", sspid.str());
	}
}

void simpleExec(std::string argstring){
	FILE* pipe;
	int pid, jobid, returnstatus;
	char buffer[512];
	std::stringstream ss, sspid;
	std::string argscommand, line, command, arg, jobdir, path, simplecommand, boxtabloc, intg;
	std::vector<std::string> commandargs;
	struct stat dirbuffer, filebuffer;
	std::ifstream unidoc;
	std::ofstream boxtab;
	
	jobid = addJob(argstring);
	jobdir = getJobFolder(argstring, jobid);
	
	if(stat(jobdir.c_str(), &dirbuffer) != 0){
		mkdir(jobdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
	
	if(getKeyValFromString(argstring, "prg") == "extract"){
		if(getKeyValFromString(argstring, "boxdir") == ""){
			missingArgumentError("boxdir");
		}
		if(stat(getKeyValFromString(argstring, "unidoc").c_str(), &filebuffer) != 0){
			fatalError("Unidoc file specified does not exist");
		}
		
		boxtabloc = jobdir;
		boxtabloc.append("/boxtab.txt");
		
		unidoc.open(getKeyValFromString(argstring, "unidoc").c_str());
		boxtab.open(boxtabloc.c_str());
		
		while(getline(unidoc,line)){
			intg = getKeyValFromString(line, "intg").c_str();
			intg.replace(intg.find(".mrc"), 4, ".box");
			boxtab << getKeyValFromString(argstring, "boxdir") << "/" << intg << "\n";
		}
		
		boxtab.close();
		unidoc.close();
		
		argstring.append(" boxtab=");
		argstring.append(boxtabloc);
	}

	pid = fork();
	
	if(pid == 0){
		path = std::getenv("PATH");
		path.append(":");
		path.append(SIMPLE_DIR_VALUE);
		path.append("/bin");
		setenv("PATH", path.c_str(), 1); 
		
		argscommand = "simple_exec prg=";
		argscommand.append(getKeyValFromString(argstring, "prg"));
	
		pipe = popen(argscommand.c_str(), "r");
		while (!feof(pipe)) {
			if (fgets(buffer, sizeof(buffer), pipe) != NULL){
				ss << buffer;
			}
		}
	
		pclose(pipe);
	
		while(std::getline(ss,line,'\n')){
			commandargs.push_back(getElementFromString(line, 0));
		}
	
		simplecommand = "simple_exec prg=";
		simplecommand.append(getKeyValFromString(argstring, "prg"));
	
		for (std::vector<std::string>::iterator it=commandargs.begin(); it!=commandargs.end(); ++it){
			arg = getKeyValFromString(argstring, *it);
			if(arg != "" && *it != ""){
				simplecommand.append(" ");
				simplecommand.append(*it);
				simplecommand.append("=");
				simplecommand.append(arg);
			}
		}
		
		simplecommand.append(" 2>&1 > job.log");
		
		chdir(jobdir.c_str());
		pipe = popen(simplecommand.c_str(), "r");
		returnstatus = pclose(pipe);
		if(WEXITSTATUS(returnstatus) == 0){
			editJob(argstring, jobid, "jobstatus", "Complete");
		}else{
			editJob(argstring, jobid, "jobstatus", "Failed");
		}
		
	}else if(pid > 0){
		sspid << pid;
		editJob(argstring, jobid, "jobpid", sspid.str());
	}
}

void simpleExecOLD(std::string argstring){
	FILE* pipe;
	int pid;
	char buffer[128];
	std::stringstream ss;
	std::string argscommand, line, command, arg, jobdir, path;
	std::vector<std::string> commandargs;
	struct stat dirbuffer;
	
	path = std::getenv("PATH");
	path.append(":/beegfs/software/electron_microscopy/simple/simple_git/bin"); // SET AT COMPILE
	setenv("PATH", path.c_str(), 1); // SET AT COMPILE
	argscommand = "simple_exec prg=";
	argscommand.append(getKeyValFromString(argstring, "prg"));
	
	pipe = popen(argscommand.c_str(), "r");
	while (!feof(pipe)) {
		if (fgets(buffer, 128, pipe) != NULL){
			ss << buffer;
		}
	}
	
	pclose(pipe);
	
	while(std::getline(ss,line,'\n')){
		commandargs.push_back(getElementFromString(line, 0));
	}
	
	jobdir = addJob(argstring);
	if(stat(jobdir.c_str(), &dirbuffer) != 0){
		mkdir(jobdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
	
	createSimpleEnv(jobdir, argstring);
	
	command = "cd ";
	command.append(jobdir);
	//command.append(" && nohup sh -c \"");
	command.append(" && ");
	command.append(argscommand);
	
	for (std::vector<std::string>::iterator it=commandargs.begin(); it!=commandargs.end(); ++it){
		arg = getKeyValFromString(argstring, *it);
		if(arg != "" && *it != ""){
			command.append(" ");
			command.append(*it);
			command.append("=");
			command.append(arg);
		}
	}
	//command.append(" && echo $? > out\"&");
	//system(command.c_str());
	command.append(" > out");
	pid = fork();
	system(command.c_str());
}

void pipelineSaveSelectionHistory(std::string argstring){
	std::string dir, selections, filename, selection;
	std::ofstream selectionfile;
	size_t pos = 0;
	
	dir = getKeyValFromString(argstring, "dir");
	selections = getKeyValFromString(argstring, "selections");
	dir.append("/.selection");
	
	selectionfile.open(dir.c_str());
	
	while ((pos = selections.find(";")) != std::string::npos) {
		selection = selections.substr(0, pos);
		if (selection != ""){
			if (selectionfile.is_open()){
				selectionfile << selection << "\n";
			}
		}
		selections.erase(0, pos + 1);
	}
	selectionfile.close();
}

void pipelineViewMaster(std::string argstring){
	std::string unidoc, returnstring, line, select;
	struct stat buffer;
	std::ifstream selectionfile, unidocstream;
	
	unidoc = getKeyValFromString(argstring, "dir");

	if(unidoc == ""){
		missingArgumentError("dir");
	}
	
	unidoc.append("/simple_unidoc.txt");
	
	returnstring = "micrographs=";
	
	if(stat(unidoc.c_str(), &buffer) == 0){
		unidocstream.open(unidoc.c_str());
		while(getline(unidocstream, line)){
			returnstring.append(getKeyValFromString(line, "intg"));
			returnstring.append(",");
			returnstring.append(getKeyValFromString(line, "pspec"));
			returnstring.append(",");
			returnstring.append(getKeyValFromString(line, "thumb"));
			returnstring.append(",");
			returnstring.append(getKeyValFromString(line, "dfx"));
			returnstring.append(",");
			returnstring.append(getKeyValFromString(line, "dfy"));
			returnstring.append(",");
			returnstring.append(getKeyValFromString(line, "angast"));
			returnstring.append(";");	
		}
		unidocstream.close();
		
		select = getKeyValFromString(argstring, "dir");
		select.append("/.selection");
	
		if(stat(select.c_str(), &buffer) == 0){
			selectionfile.open(select.c_str());
			if (selectionfile.is_open()){
				returnstring.append(" selected=");
				while(getline(selectionfile,line)){
					returnstring.append(line);
					returnstring.append(";");
				}
				selectionfile.close();
			}
		}
		std::cout << returnstring << std::endl;
	}
}

void pipelineViewMasterOLD(std::string argstring){
	std::string dir, returnstring, line, select, ctf, ctfparams;
	DIR *directory;
	struct dirent *entry;
	struct stat buffer, buffer2;
	std::vector<std::string> pipelineviewlist;
	std::ifstream selectionfile;
	
	dir = getKeyValFromString(argstring, "dir");

	if(dir == ""){
		missingArgumentError("dir");
	}
	
	returnstring = "micrographs=";
	
	if(stat(dir.c_str(), &buffer) == 0){
		directory = opendir(dir.c_str());
		entry = readdir(directory);
		while (entry != NULL){
			if(strstr(entry->d_name, ".mrc") && strstr(entry->d_name, "intg")){
				//pipelineviewlist.push_back(entry->d_name);
				returnstring.append(entry->d_name);
				ctf = entry->d_name;
				ctf.replace(ctf.find(".mrc"), 4, ".txt");
				ctf.replace(0,ctf.find("intg"), "");
				ctf.replace(ctf.find("intg"), 4, "ctffind_output_part");
				ctf.insert(0, "/");
				ctf.insert(0, dir);
				ctfparams = readCTFParams(ctf);
				returnstring.append(",");
				returnstring.append(getKeyValFromString(ctfparams, "dfx"));
				returnstring.append(",");
				returnstring.append(getKeyValFromString(ctfparams, "dfy"));
				returnstring.append(",");
				returnstring.append(getKeyValFromString(ctfparams, "angast"));
				returnstring.append(";");
			}
			entry = readdir(directory);
		}
		closedir(directory);
    //	sort(pipelineviewlist.begin(), pipelineviewlist.end());

    //	returnstring = "micrographs=";
    //	for (std::vector<std::string>::iterator it=pipelineviewlist.begin(); it!=pipelineviewlist.end(); ++it){
    //		returnstring.append(*it);
    //		returnstring.append(";");
    //	}
		
		select = dir;
		select.append("/.selection");
	
		if(stat(select.c_str(), &buffer2) == 0){
			selectionfile.open(select.c_str());
			if (selectionfile.is_open()){
				returnstring.append(" selected=");
				while(getline(selectionfile,line)){
					returnstring.append(line);
					returnstring.append(";");
				}
				selectionfile.close();
			}
		}
		std::cout << returnstring << std::endl; // MAYBE GET CTF INFO TOO SO CAN SORT
	}
}

void pipelineViewSlave(std::string argstring){
	std::string dir, micrographs, micrograph, thumb, ctf, pspec, box, pspecbrightness, pspeccontrast, thumbbrightness, thumbcontrast;
	std::vector<std::string> boxes;
	size_t pos = 0;
	
	dir = getKeyValFromString(argstring, "dir");
	micrographs = getKeyValFromString(argstring, "micrographs");
	pspecbrightness = getKeyValFromString(argstring, "pspecbrightness");
	pspeccontrast = getKeyValFromString(argstring, "pspeccontrast");
	thumbbrightness = getKeyValFromString(argstring, "thumbbrightness");
	thumbcontrast = getKeyValFromString(argstring, "thumbcontrast");
	
	pipelineSaveSelectionHistory(argstring);

	if(dir == ""){
		missingArgumentError("dir");
	}
	
  if(micrographs == ""){
		missingArgumentError("micrographs");
	}
	
	if(pspecbrightness == ""){
		pspecbrightness = "64";
	}
	
	if(pspeccontrast == ""){
		pspeccontrast = "4";
	}
	
	if(thumbbrightness == ""){
		thumbbrightness = "128";
	}
	
	if(thumbcontrast == ""){
		thumbcontrast = "2";
	}
	
	dir.append("/");
	
	while ((pos = micrographs.find(";")) != std::string::npos) {
		micrograph = micrographs.substr(0, pos);
	
		thumb = micrograph;
		thumb.replace(thumb.find("intg"), 4, "thumb");
		thumb.insert(0, dir);
	
		pspec = micrograph;
		pspec.replace(pspec.find("intg"), 4, "pspec");
		pspec.insert(0, dir);
		
    //	ctf = micrograph;
    //	ctf.replace(ctf.find(".mrc"), 4, ".txt");
    //	ctf.replace(0,ctf.find("intg"), "");
    //	ctf.replace(ctf.find("intg"), 4, "ctffind_output_part");
    //	ctf.insert(0, dir);
		
		box = micrograph;
		box.replace(box.find(".mrc"), 4, ".box");
		box.insert(0, dir);
		
    //	std::cout << "ctf=" << micrograph << " " << readCTFParams(ctf) << std::endl;
		std::cout << "pspec=" << micrograph << " image=" << pngFromMRC(pspec, std::strtof(pspecbrightness.c_str(), NULL), std::strtof(pspeccontrast.c_str(), NULL)) << std::endl;
		std::cout << "thumb=" << micrograph << " image=" << pngFromMRC(thumb, std::strtof(thumbbrightness.c_str(), NULL), std::strtof(thumbcontrast.c_str(), NULL)) << std::endl;
		boxes = readBoxes(box);
		std::cout << "box="<< micrograph << " count=" << boxes.size() << std::endl;
		micrographs.erase(0, pos + 1);
	}
}

void pipelineSaveSelection(std::string argstring){
	std::string dir, selections, filename, selection;
	std::ofstream selectionfile;
	
	size_t pos = 0;
	
	pipelineSaveSelectionHistory(argstring);
	
	dir = getKeyValFromString(argstring, "dir");
	filename = getKeyValFromString(argstring, "filename");
	selections = getKeyValFromString(argstring, "selections");
	
	selectionfile.open(filename.c_str());
	
	while ((pos = selections.find(";")) != std::string::npos) {
		selection = selections.substr(0, pos);
		if (selection != ""){
			if (selectionfile.is_open()){
				selectionfile << dir << "/" << selection << "\n";
			}
		}
		selections.erase(0, pos + 1);
	}
	selectionfile.close();
}

void pipelineViewBoxes(std::string argstring){
	std::string micrograph, boxfile, thumbbrightness, thumbcontrast, returnstring, line;
	std::vector<std::string> boxes;
	int nx, ny, nz, mode;
	struct stat buffer;
	bool nextspace;
	std::ostringstream nxss, nyss;
	
	micrograph = getKeyValFromString(argstring, "dir");
	micrograph.append("/");
	micrograph.append(getKeyValFromString(argstring, "micrograph"));
	
	thumbbrightness = getKeyValFromString(argstring, "thumbbrightness");
	thumbcontrast = getKeyValFromString(argstring, "thumbcontrast");
	
	boxfile = micrograph;
	boxfile.replace(boxfile.find(".mrc"), 4, ".box");
	boxes = readBoxes(boxfile);
	
	returnstring = "";
	
	if(stat(micrograph.c_str(), &buffer) == 0){
		readMRCHeader(micrograph, &nx, &ny, &nz, &mode);
		returnstring.append("nx=");
		nxss << nx;
		returnstring.append(nxss.str());
		returnstring.append(" ny=");
		nyss << ny;
		returnstring.append(nyss.str());
		returnstring.append(" image=");
		returnstring.append(pngFromMRC(micrograph, std::strtof(thumbbrightness.c_str(), NULL), std::strtof(thumbcontrast.c_str(), NULL)));
		
	}
	if(stat(boxfile.c_str(), &buffer) == 0){
		returnstring.append(" boxes=");
		for (std::vector<std::string>::iterator it=boxes.begin(); it!=boxes.end(); ++it){
			nextspace = false;
			line = *it;
			for(unsigned count=0; count < line.length(); count++) {
				if(line[count] != ' ') {
					nextspace = true;
					returnstring.push_back(line[count]);
				} else {
					if(nextspace){
						returnstring.append(",");
					}
					nextspace = false;
				}
			}
			returnstring.append(";");
		}
	}
	std::cout << returnstring << std::endl;
}

void import(std::string argstring){
	std::string directory, suffix, outfile, jobdir;
	struct stat buffer;
	DIR *directorycontents;
	struct dirent *entry;
	std::vector<std::string> files;
	std::ofstream outputfile;
	int jobid;
	
	directory = getKeyValFromString(argstring, "dir");
	suffix = getKeyValFromString(argstring, "suffix");
	outfile = getKeyValFromString(argstring, "outfile");
	
	if(outfile == ""){
		outfile = "import_" + suffix + ".txt";
	}
	
	jobid = addJob(argstring);
	jobdir = getJobFolder(argstring, jobid);
	
	if(stat(jobdir.c_str(), &buffer) != 0){
		mkdir(jobdir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
	
	jobdir.append("/");
	jobdir.append(outfile);
	
	if(stat(directory.c_str(), &buffer) == 0){
		directorycontents = opendir(directory.c_str());
		entry = readdir(directorycontents);
		while (entry != NULL){
			if(entry->d_name[0] != '.'){
				if(std::strstr(entry->d_name, suffix.c_str())){
					files.push_back(entry->d_name);
				}
			}
			entry = readdir(directorycontents);
		}
		
		closedir(directorycontents);
		
		std::sort (files.begin(), files.end());
		
		outputfile.open(jobdir.c_str());
		
		for (std::vector<std::string>::iterator it=files.begin(); it!=files.end(); ++it){
			outputfile << directory << "/" << *it << "\n";
		}
		
		outputfile.close();
		editJob(argstring, jobid, "jobstatus", "Complete");
	}else{
		editJob(argstring, jobid, "jobstatus", "Error");
		fatalError("Directory doesn't exist");
	}
}

void showFile(std::string argstring){
	std::string file, contrast, brightness, returnstring, line;
	struct stat buffer;
	int nx, ny, nz, mode;
	std::ifstream logfile, txtfile;
	
	file = getKeyValFromString(argstring, "file");
	contrast = getKeyValFromString(argstring, "contrast");
	brightness = getKeyValFromString(argstring, "brightness");
	
	if(stat(file.c_str(), &buffer) == 0){
		if(std::strstr(file.c_str(), ".mrc")){
			readMRCHeader(file, &nx, &ny, &nz, &mode);
			if(nz == 1 && mode == 2){
				returnstring = "";
				returnstring.append("image=");
				returnstring.append(pngFromMRC(file, std::strtof(brightness.c_str(), NULL), std::strtof(contrast.c_str(), NULL)));
				std::cout << returnstring << " nz=" << nz  << " nx=" << nx << " ny=" << ny << std::endl;
			}else if(nz > 1 && mode == 2){
				returnstring = "";
				pngsFromMRCStack(file, std::strtof(brightness.c_str(), NULL), std::strtof(contrast.c_str(), NULL));
			}
		}else if(std::strstr(file.c_str(), ".log")){
			logfile.open(file.c_str());
			returnstring.append("log=");
			while(getline(logfile, line)){
				std::replace(line.begin(), line.end(), ' ', '^');
				returnstring.append(line);
			}
			logfile.close();
			std::cout << returnstring << std::endl;
		}else if(std::strstr(file.c_str(), ".txt")){
			txtfile.open(file.c_str());
			returnstring.append("txt=");
			while(getline(txtfile, line)){
				std::replace(line.begin(), line.end(), ' ', '^');
				std::replace(line.begin(), line.end(), ' ', '|');
				returnstring.append(line);
			}
			txtfile.close();
			std::cout << returnstring << std::endl;
		}
	}
}

void filesViewSaveSelection(std::string argstring){
	std::string outfilename, selections, filename, selection;
	std::vector<int> selectionarray;
	std::ifstream infile; 
	std::ofstream outfile;
	char *buffer;
	int *bufferint, nx, ny, nz, zcount, framesize;
	
	size_t pos = 0;
	
	filename = getKeyValFromString(argstring, "filename");
	outfilename = getKeyValFromString(argstring, "outfile");
	selections = getKeyValFromString(argstring, "selections");
	
	while ((pos = selections.find(";")) != std::string::npos) {
		selection = selections.substr(0, pos);
		selectionarray.push_back(std::atoi(selection.c_str()));
		selections.erase(0, pos + 1);
	}
	
	infile.open(filename.c_str(), std::ifstream::binary);
	outfile.open(outfilename.c_str(), std::ofstream::binary);
	
	if (infile.is_open()){
		buffer = new char [1024];
		
		infile.read(buffer, 1024);
		
		bufferint = (int*)buffer;
		nx = bufferint[0];
		ny = bufferint[1];
		nz = bufferint[2];
		bufferint[2] = selectionarray.size();
		outfile.write(buffer, 1024);
		
		delete[] buffer;
		
		infile.seekg(1024, std::fstream::beg);
		zcount = 0;
		framesize = nx * ny * sizeof(double);
		
		while(zcount < nz){
			buffer = new char[framesize];
			infile.read(buffer, framesize);
			for (std::vector<int>::iterator it=selectionarray.begin(); it!=selectionarray.end(); ++it){
				if(zcount == *it){
					outfile.write(buffer, framesize);
					break;
				}
			}
			delete[] buffer;
			zcount++;
		}
	}
	infile.close();
	outfile.close();

}

int main(int argc, char* argv[]){ 														// Main
	std::string argstring, cmd;
	
	if(argc > 1){
		for(int i = 0; i < argc; i++){
			argstring.append(argv[i]);
			argstring.append(" ");
		}
	}else{
		std::getline(std::cin, argstring);
	}
	
	cmd = getKeyValFromString(argstring, "cmd");
	
	if(cmd != ""){
		if(cmd == "dirlist"){
			listDirectory(argstring);
		}else if(cmd == "projectlist"){
			listProjects(argstring);
		}else if(cmd == "addproject"){
			addProject(argstring);
		}else if(cmd == "jobslist"){
			listJobs(argstring);
		}else if(cmd == "addjob"){
			addExternalJob(argstring);
		}else if(cmd == "simple_distr_exec"){
			simpleDistrExec(argstring);
		}else if(cmd == "simple_exec"){
			simpleExec(argstring);	
		}else if(cmd == "pipelineviewmaster"){
			pipelineViewMaster(argstring);
		}else if(cmd == "pipelineviewslave"){
			pipelineViewSlave(argstring);
		}else if(cmd == "pipelinesaveselection"){
			pipelineSaveSelection(argstring);
		}else if(cmd == "pipelineviewboxes"){
			pipelineViewBoxes(argstring);
		}else if(cmd == "import"){
			import(argstring);
		}else if(cmd == "showfile"){
			showFile(argstring);
		}else if(cmd == "killjob"){
			killJob(argstring);
		}else if(cmd == "deletejob"){
			deleteJob(argstring);
		}else if(cmd == "jobparameters"){
			listJob(argstring);
		}else if(cmd == "savestackselection"){
			filesViewSaveSelection(argstring);
		}
	}else{
		missingArgumentError("cmd");
	}
	completionSuccess();
	return 0;
}
