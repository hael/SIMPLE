// open include guard

#ifndef __PROJECT_H_INCLUDED__
#define __PROJECT_H_INCLUDED__

//=================================


// included external dependencies

#include <string>  
#include <vector>  
#include <sqlite3.h> 
 
//=================================


// included internal dependencies

#include "simple_gui_sqlite.h" 

//=================================


//  List projects

void listProjects(struct WebResponse* webresponse) {
	
	std::vector<std::map<std::string, std::string> >	result;
	std::string											projects;
	int													rowid;
	int													resultcount;
	
	SQLQuery (result, "SELECT * FROM projects", rowid);
	
	projects = "[";
	
	for(resultcount = 0; resultcount < result.size(); resultcount++){
		projects += "[";
		projects += "\"" + result[resultcount]["projectname"] + "\",";
		projects += "\"" + result[resultcount]["projectdescription"] + "\",";
		projects += "\"" + result[resultcount]["projectfolder"] + "\"";
		projects += "],";
	}
	projects.pop_back();
	projects += "]";
	webresponse->responsemap["projects"] = projects;
}

//=================================


//  Create project

void creatroject(struct WebResponse* webresponse) {


}

//=================================


// close include guard

#endif

//=================================
