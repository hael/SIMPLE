// open include guard

#ifndef __LOGGING_H_INCLUDED__
#define __LOGGING_H_INCLUDED__

//=================================


// included dependencies

#include <string>
#include <ctime>
#include <iostream>
#include "simple_gui_environment.h"
//=================================

//  check environment

void log(std::string logtext) {
	std::string			logstring;
	time_t 				date;
	
	date = time(0);
	logstring = ctime(&date);
	logstring.pop_back();												// Removetrailing newline
	logstring += " : " + logtext;
		
	if(environment.verbose){
		std::cout << logstring << std::endl;
	}
}

void error(std::string errortext) {
	std::string			errorstring;
	time_t 				date;
	
	date = time(0);
	errorstring = ctime(&date);
	errorstring.pop_back();												// Removetrailing newline
	errorstring += " : ERROR :" + errortext;
		
	if(environment.verbose){
		std::cout << errorstring << std::endl;
	}
	exit(1);
}
//=================================


// close include guard

#endif

//=================================
