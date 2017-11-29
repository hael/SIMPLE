#ifndef SIMPLE_GUI_SQL_H_INCLUDED
#define SIMPLE_GUI_SQL_H_INCLUDED

#include <sqlite3.h> 
#include <string>

void SQLQuery(char* databasepath, std::string& returnstring, std::string sql, int &rowid){
	sqlite3*		db;
	sqlite3_stmt*	stmt;
	int 			rc;
	int				id;
	
	
	std::cout << sql << std::endl;
	
	rc = sqlite3_open(databasepath, &db);
	if(rc){
		std::cout << "Failed to open database " << sqlite3_errmsg(db) << std::endl;
	}
	
	rc = sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, NULL);
	if( rc != SQLITE_OK ){
		std::cout << "Failed to query database " << sqlite3_errmsg(db) << std::endl;
	}
	returnstring = " ";
	while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
		returnstring += "{";
		int id = sqlite3_column_int(stmt, 0);
		returnstring += "\"id\":\"";
		returnstring += std::to_string(id);
		returnstring += "\",";
		for(int i = 1; i < sqlite3_column_count(stmt); i++){
			const char *columnname = sqlite3_column_name(stmt, i);
			const unsigned char *columnvalue = sqlite3_column_text(stmt, i);
			returnstring += "\"";
			returnstring += std::string(columnname);
			returnstring += "\":\"";
			returnstring += std::string((char*)columnvalue);
			returnstring += "\",";
		}
		returnstring.pop_back();
		returnstring += "},";
	}
	returnstring.pop_back();
	
	if (rc != SQLITE_DONE) {
		std::cout << "Failed to parse query " << sqlite3_errmsg(db) << std::endl;
	}

	rc = sqlite3_prepare_v2(db, "SELECT last_insert_rowid();", -1, &stmt, NULL);
	if( rc != SQLITE_OK ){
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

#endif
