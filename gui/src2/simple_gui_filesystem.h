// open include guard

#ifndef __FILESYSTEM_H_INCLUDED__
#define __FILESYSTEM_H_INCLUDED__

//=================================


// included dependencies

#include <string>
#include <sys/stat.h>
#include <pwd.h>

//=================================


//  check file or folder exists

bool pathExists(std::string path) {
	struct stat 		buf;

	if(stat(path.c_str(), &buf) == 0) {
		return(true);
	} else {
		return(false);
	}
}

//=================================


//  get home directory of user running the server

void getUserHome(std::string& databasepath) {
	
	struct passwd* 		pw;
	uid_t 				uid;
	
	uid = geteuid();
	pw = getpwuid(uid);
	databasepath = std::string(pw->pw_dir);
	
}
//=================================


// close include guard

#endif

//=================================
