// open include guard

#ifndef __SQLITE_H_INCLUDED__
#define __SQLITE_H_INCLUDED__

//=================================


// included external dependencies

#include <string>
#include <vector>

#include <sqlite3.h>

//=================================


// included internal dependencies

#include "simple_gui_environment.h"
#include "simple_gui_filesystem.h"
#include "simple_gui_logging.h"

//=================================


//  Get location of the database

void getDatabaseLocation (std::string& databasepath) {
	
	struct passwd* 		pw;
	uid_t 				uid;
	
	uid = geteuid();
	pw = getpwuid(uid);
	
	if(environment.multiuser){
		databasepath = std::string(pw->pw_dir) + "/.simple." + environment.user + ".sqlite";
	} else {
		databasepath = std::string(pw->pw_dir) + "/.simple.sqlite";
	}
}

//=================================


//  SQL query

void SQLQuery(std::vector<std::map<std::string, std::string> > &result, std::string sql, int& rowid){
	
	std::string											databasepath;
	sqlite3*											db;
	sqlite3_stmt*										stmt;
	int 												rc;
	int													id;
	std::map<std::string, std::string>					resultlinekeyval;
	
	getDatabaseLocation(databasepath);
	
	if(sqlite3_open(databasepath.c_str(), &db) != SQLITE_OK){
		log("Failed to open sqlite database");
		log("Database path " + databasepath);
		log(sqlite3_errmsg(db));
	}
	
	sqlite3_busy_timeout(db, 5000);
	
	if(sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, NULL) != SQLITE_OK ){
		log("Failed to query sqlite database");
		log(sqlite3_errmsg(db));
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
		log("Failed to parse sqlite query result");
		log(sqlite3_errmsg(db));
	}

	if(sqlite3_prepare_v2(db, "SELECT last_insert_rowid();", -1, &stmt, NULL)  != SQLITE_OK ){
		log("Failed to obtain last insert row id");
		log(sqlite3_errmsg(db));
	}	
	
	while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
		rowid = sqlite3_column_int(stmt, 0);
	}
	
	if (rc != SQLITE_DONE) {
		log("Failed to query sqlite database");
		log(sqlite3_errmsg(db));
	}

	sqlite3_finalize(stmt);
	sqlite3_close(db);
}

//=================================


// close include guard

#endif

//=================================
