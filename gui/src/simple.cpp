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

std::string getUser(){ 											// Return user or daemon user if local
	std::string user;
	struct passwd* pw;
	register uid_t uid;
	
	uid = geteuid();
	pw = getpwuid(uid);
	user = pw->pw_name;
	
	return user;
}

std::string getUserHome(std::string user){ 												// Return home directory of user
	std::string userhome;
	struct passwd* pw;
	
	pw = getpwnam(user.c_str());
	userhome = pw->pw_dir;
	
	return userhome;
}

std::string getConfOption(std::string key){
	std::string conffile, conffileline, returnstring;
	std::ifstream conffilehandle;
	
	conffile = SIMPLE_DIR_VALUE;
	conffile.append("/gui/etc/simple_gui.conf");

	returnstring = "";
	
	conffilehandle.open(conffile.c_str());
	
	while (std::getline(conffilehandle, conffileline)){
		if(getKeyValFromString(conffileline, key) != ""){
			returnstring = getKeyValFromString(conffileline, key);
			break;
		}
	}
	
	conffilehandle.close();
	
	return returnstring;
}

void checkRunning(std::string port){
	std::string pidfile, pid;
	std::ifstream pidfilehandle;
	
	pidfile = SIMPLE_DIR_VALUE;
	pidfile.append("/gui/.pid");
	
	pidfilehandle.open(pidfile.c_str());
	
	if(pidfilehandle.is_open()){
		std::getline(pidfilehandle, pid);
		if (kill(std::atoi(pid.c_str()), 0) == 0){
			std::cout << "A Simple GUI instance is already running on port " << port << "!" << std::endl;
			exit(1);
		}
	}
	
	pidfilehandle.close();
}

void checkNotRunning(std::string port){
	std::string pidfile, pid;
	std::ifstream pidfilehandle;
	
	pidfile = SIMPLE_DIR_VALUE;
	pidfile.append("/gui/.pid");
	
	pidfilehandle.open(pidfile.c_str());
	
	if(pidfilehandle.is_open()){
		std::getline(pidfilehandle, pid);
		if (kill(std::atoi(pid.c_str()), 0) != 0){
			std::cout << "The Simple GUI server is not running. Please start it using \"simple start\" "<< std::endl;
			exit(1);
		}
	}
	
	pidfilehandle.close();
}

void killRunning(){
	std::string pidfile, pid;
	std::ifstream pidfilehandle;
	
	pidfile = SIMPLE_DIR_VALUE;
	pidfile.append("/gui/.pid");
	
	pidfilehandle.open(pidfile.c_str());
	
	if(pidfilehandle.is_open()){
		std::getline(pidfilehandle, pid);
		if (kill(std::atoi(pid.c_str()), 9) == 0){
			std::cout << "Simple GUI server killed" << std::endl;
			exit(0);
		}
	}
	
	pidfilehandle.close();
}

void openBrowser(std::string port){
	std::string browsercommand;
	
	browsercommand = getConfOption("browser_command");
	browsercommand.append(" http://localhost:");
	browsercommand.append(port);
	system(browsercommand.c_str());
}

void startDaemon(std::string port){

	std::string websocketcommand;
	
	websocketcommand = SIMPLE_DIR_VALUE;
	websocketcommand.append("/gui/bin/websocketd --staticdir=");
	websocketcommand.append(SIMPLE_DIR_VALUE);
	websocketcommand.append("/gui/www --port=");
	websocketcommand.append(port);
	websocketcommand.append(" ");
	websocketcommand.append(SIMPLE_DIR_VALUE);
	websocketcommand.append("/gui/bin/simple_ws & echo $! > ");
	websocketcommand.append(SIMPLE_DIR_VALUE);
	websocketcommand.append("/gui/.pid");
	
	system(websocketcommand.c_str());

	openBrowser(port);

}

int main(int argc, char* argv[]){ 														// Main
	std::string user, userconfdir, browser, port;
	
	port  = getConfOption("port");

	if(argc > 1){
		for(int i = 0; i < argc; i++){
			if(strcmp(argv[1], "start") == 0){
				checkRunning(port);
				startDaemon(port);
			}
			if(strcmp(argv[1], "stop") == 0){
				killRunning();
			}
		}
	} else {
		checkNotRunning(port);
		openBrowser(port);
	}
	
	
	
	
	return 0;
}
