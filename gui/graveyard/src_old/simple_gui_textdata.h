#ifndef SIMPLE_GUI_TEXTDATA_H_INCLUDED
#define SIMPLE_GUI_TEXTDATA_H_INCLUDED

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

void getElementFromString(std::string *output, std::string argstring, int element){ 
	std::size_t location;
	char character;
	int currentelement;
	
	output->clear();
	location = 0;
	currentelement = 0;
	
	while(argstring.c_str()[location] == ' ' && location < argstring.length()){
		location = location + 1;
	}
	
	while ( location < argstring.length() ){
		character = argstring.c_str()[location];
		if ( character != ' '){
			if(currentelement == element){
				output->push_back(character);
			}
		} else {
			currentelement = currentelement + 1;
		}
		location = location + 1;
	}
}

#endif

