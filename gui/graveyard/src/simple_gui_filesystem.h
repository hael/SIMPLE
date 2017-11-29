#ifndef SIMPLE_GUI_FILESYSTEM_H_INCLUDED
#define SIMPLE_GUI_FILESYSTEM_H_INCLUDED

#include <pwd.h>

void listDirectory (struct mg_connection* http_connection, struct http_message* message) {
	
 std::string 	JSONstring;
 std::string 	directorystring;
 std::string 	filestring;
 char         	directoryname[1024];
 char         	filefilter[128];
 const char*	JSONchar;
 struct stat 	statbuffer;
 DIR*			directorycontents;
 struct dirent*	entry;
 
 if (mg_get_http_var(&message->query_string, "directoryname", directoryname, sizeof(directoryname)) <= 0) {
	std::strcpy(directoryname, "/");
 }

 if (mg_get_http_var(&message->query_string, "filefilter", filefilter, sizeof(filefilter)) <= 0) {
	std::strcpy(filefilter, "*");
 }
	JSONstring = "{";
	if(stat(directoryname, &statbuffer) != 0){
		JSONstring += "\"error\":\"Directory does not exist\"}";
		mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
		mg_send(http_connection, JSONchar, JSONstring.length());
		return;
	}
	JSONstring += "\"directory\":\"";
	JSONstring += directoryname;
	JSONstring += "\",";
	directorycontents = opendir(directoryname);
	entry = readdir(directorycontents);
	while (entry != NULL){
		if (entry->d_type == DT_DIR){
			if(entry->d_name[0] != '.'){
				directorystring += "\"";
				directorystring += entry->d_name;
				directorystring += "\",";
			} 
		}else{
			if(entry->d_name[0] != '.'){
				std::cout << filefilter << std::endl;
				if (std::strcmp(filefilter, "*") == 0) {
					filestring += "\"";
					filestring += entry->d_name;
					filestring += "\",";
				} else if (std::strstr(entry->d_name, filefilter) != NULL) {
					filestring += "\"";
					filestring += entry->d_name;
					filestring += "\",";
				}
			}
		}
		entry = readdir(directorycontents);
	}
	closedir(directorycontents);
	if (directorystring.length() > 0){
		directorystring.pop_back();
	}
	if (filestring.length() > 0){
		filestring.pop_back();
	}
	JSONstring += "\"directories\":[" + directorystring + "],";
	JSONstring += "\"files\":[" + filestring + "]}";
	JSONchar = JSONstring.c_str();
	mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
	mg_send(http_connection, JSONchar, JSONstring.length());
	return;

}

void databaseLocation (struct http_message* message, char* databasepath) {
	
	struct 		passwd* 	pw;
	uid_t 		uid;
	
	uid = geteuid();
	pw = getpwuid(uid);
	strcpy(databasepath, pw->pw_dir);
	strcat(databasepath, "/.simple.sqlite");
}

#endif

