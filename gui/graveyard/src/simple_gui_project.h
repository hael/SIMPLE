#ifndef SIMPLE_GUI_PROJECT_H_INCLUDED
#define SIMPLE_GUI_PROJECT_H_INCLUDED

#include "simple_gui_filesystem.h"
#include "simple_gui_sql.h"

#include <unistd.h>
#include <cstring>
#include <string>
#include <stdlib.h> 

extern "C" {
    #include "ext/mongoose/mongoose.h"
}

void newProject (struct mg_connection* http_connection, struct http_message* message) {
	char         projectname[128];
	char         projectdescription[2048];
	char         projectfolder[2048];
	std::string  JSONstring;
	std::string  projecttable = "p" + std::to_string(rand());
	const char*	JSONchar;
	char		databasepath[1024];
	int			rowid;
	 
	if (mg_get_http_var(&message->query_string, "projectname", projectname, sizeof(projectname)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "projectdescription", projectdescription, sizeof(projectdescription)) <= 0){
		//error
	}
	if (mg_get_http_var(&message->query_string, "projectfolder", projectfolder, sizeof(projectfolder)) <= 0){
		//error
	}
	
	databaseLocation(message, databasepath);
	
	SQLQuery(databasepath, JSONstring, "CREATE TABLE IF NOT EXISTS projects (\
										projectid integer PRIMARY KEY,\
										projectname text NOT NULL,\
										projectdescription text NOT NULL,\
										projectfolder text NOT NULL,\
										projecttable text NOT NULL,\
										projectusage bool NOT NULL);", rowid);
	
				
	SQLQuery(databasepath, JSONstring, "INSERT INTO projects(\
										projectname,\
										projectdescription,\
										projectfolder,\
										projecttable,\
										projectusage)\
										VALUES('"
										+ std::string(projectname) + "','"
										+ std::string(projectdescription) + "','"
										+ std::string(projectfolder) + "','"
										+ projecttable + "','"
										+ "true" + "');", rowid);
	
	JSONchar = JSONstring.c_str();
	mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
	mg_send(http_connection, JSONchar, JSONstring.length());

}


void getProjects (struct mg_connection* http_connection, struct http_message* message) {
	std::string  JSONstring;
	std::string  sqlresult;
	const char*	 JSONchar;
	char		databasepath[1024];
	int			rowid;
	
	databaseLocation(message, databasepath);
	
	SQLQuery(databasepath, sqlresult, "SELECT * FROM projects", rowid);
	
	JSONstring = "{\"projects\": [";
	JSONstring += sqlresult;
	JSONstring += "]}";
	JSONchar = JSONstring.c_str();

	mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
	mg_send(http_connection, JSONchar, JSONstring.length());	
}
#endif
