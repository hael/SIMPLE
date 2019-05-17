// open include guard

#ifndef __ENVIRONMENT_H_INCLUDED__
#define __ENVIRONMENT_H_INCLUDED__

//=================================


// included external dependencies

#include <string> 

//=================================


// struct declaration

extern struct Environment {
	std::string			simplepath;
	std::string			databasepath;
	std::string			logpath;
	std::string			webpath;
	std::string			user;
	int					port = 8088;
	bool				multiuser;
	bool				verbose;
} environment;

// included internal dependencies

#include "simple_gui_logging.h"
#include "simple_gui_filesystem.h"

//=================================


//  check environment

void environmentCheck() {
	
	if(getenv("SIMPLE_PATH") != NULL){
		environment.simplepath = getenv("SIMPLE_PATH");
		log("SIMPLE_PATH is set to " + environment.simplepath);
    } else {
		error("SIMPLE_PATH is not set");
	}
	
	if(pathExists(environment.simplepath)){
		log("Simple path is set to " + environment.simplepath);
	} else {
		error("Simple path (" + environment.simplepath + ") does not exist");
	}
	
	environment.logpath = environment.simplepath + "/log";
	
	if(pathExists(environment.logpath)){
		log("Log path is set to " + environment.logpath);
	} else {
		error("Log path (" + environment.logpath + ") does not exist");
	}
	
	environment.webpath = environment.simplepath + "/www";
	
	if(pathExists(environment.webpath)){
		log("Web path is set to " + environment.webpath);
	} else {
		error("Web path (" + environment.webpath + ") does not exist");
	}
	
	getUserHome(environment.databasepath);
	
	if(pathExists(environment.databasepath)){
		log("Database path is set to " + environment.databasepath);
	} else {
		error("Database path (" + environment.databasepath + ") does not exist");
	}
	
	if(environment.multiuser){
		log("Multi user mode is enabled");
	} else {
		log("Single user mode is enabled");
	}
	
	log("Port is set to " +  std::to_string(environment.port));
	
	
}

//=================================


//  setup environment

void environmentSetup(int& argc, char* argv[]) {
	int				option;
	
	while( ( option = getopt (argc, argv, "mvp:") ) != -1 ) {
        switch(option){
            case 'm':
                environment.multiuser = true;
                break;
            case 'v':
                environment.verbose = true;
                break;
			case 'p':
                environment.port = std::atoi(optarg);
                break;      
        }
    }

}

//=================================


// close include guard

#endif

//=================================


