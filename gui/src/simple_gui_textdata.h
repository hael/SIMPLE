#ifndef SIMPLE_GUI_TEXTDATA_H_INCLUDED
#define SIMPLE_GUI_TEXTDATA_H_INCLUDED

#include <string>

void getKeyValFromString(std::string *output, std::string argstring, std::string key){ 				
	
	std::size_t		location;
	char 			character;
	
	output->clear();
	key.append("=");
	location = argstring.find(key) + key.length();
	if(location != key.length() - 1){
		while ( location < argstring.length() ){
			character = argstring.c_str()[location];
			if (character == '^'){
				*output += " ";
			} else if (character != ' '){
				*output += character;
			} else {
				break;
			}
			location = location + 1;
		}
	}
}

std::string getElementFromString(std::string argstring, int element){ 
	std::size_t 	location;
	char 			character;
	int 			currentelement;
	int 			nextelement;
	std::string		returnstring;

	location = argstring.find("  ");
	while (location != std::string::npos) {
		argstring.replace (location, 2, " ");
		location = argstring.find("  ");
	}

	location = 0;
	
	while(argstring.c_str()[0] == ' '){
		argstring.erase (0, 1);
	}
	
	currentelement = 0;
	
	while (currentelement  < element){
		argstring.erase(0, argstring.find(" ") + 1);
		currentelement += 1;
	}
	
	returnstring = argstring.substr(0, argstring.find(" "));
	
	
	/*
	
	while (currentelement  < element){
		location = argstring.find(currentelement + 1, " ");
		currentelement += 1;
	}
	
	nextelement = argstring.find(location + 1, " ");
	
	/*std::cout << argstring << std::endl;
	
	location = 0;
	
	while(location < argstring.length()){
		if(argstring.c_str()[location] == '\t'){
			argstring.replace (location, 1, " ");
		}
		location += 1;
	}
	
	location = 0;
	
	while(location < argstring.length()){
		if(argstring.c_str()[location] == ' ' && argstring.c_str()[location + 1] == ' '){
			argstring.erase (location,1);
		}else {
			location += 1;
		}
	}
	
	//std::cout << "cut " <<argstring << std::endl;*/

/*	location = 0;
	currentelement = 0;
	
	while(argstring.c_str()[location] == ' ' && location < argstring.length()){
		location = location + 1;
	}
	
	while ( location < argstring.length() ){
		character = argstring.c_str()[location];
		if ( character == ' ' ){
			while (character == ' '&& location < argstring.length()){
				location += 1;
			}
			currentelement += 1;
		}
		if(currentelement == element){
			returnstring.push_back(character);
		}
		location += 1;
		
		
		if ( character != ' ' && character != '\t'){
			if(currentelement == element){
				returnstring.push_back(character);
			}
		} else {
			currentelement += 1;
		}
		location += 1;
	}*/
	
	return returnstring;
}

#endif

