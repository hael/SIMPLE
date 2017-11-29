#ifndef SIMPLE_GUI_STARDATA_H_INCLUDED
#define SIMPLE_GUI_STARDATA_H_INCLUDED

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

struct StarFile {
	std::vector<std::string> columns;
	std::vector<std::string> data;
};

void readStarFile(StarFile* starfile, char* filename){
	
	std::ifstream 	input;
	std::string 	line;
	std::size_t 	found;

	input.open(filename);
	if(input.is_open()){
		while(getline(input, line)){
			if(line.length() > 0){
				if(line.at(0) == '_'){
					line.erase(0,1);
					found = line.find("#");
					if(found != std::string::npos){
						line.erase(line.begin() + found, line.end());
					}
					found = line.find(" ");
					if(found != std::string::npos){
						line.erase(line.begin() + found, line.end());
					}
					starfile->columns.push_back(line);
				} else if(line.at(0) == '#' || line.at(0) == ' '){
				} else if(line.find("data_") != std::string::npos || line.find("loop_") != std::string::npos) {
				} else {
					starfile->data.push_back(line);
				}
			}
		}
		input.close();
	} else {
		std::cout << "no starfil" << std::endl;
	}
}

static void writeStarFile(StarFile *starfile, const char* filename){
	std::ofstream output;
	std::vector<std::string>::iterator it;
	
	output.open(filename);
	if(output.is_open()){
		output << "data_" << "\n";
		output << "loop_" << "\n";
		for (it = starfile->columns.begin() ; it != starfile->columns.end(); ++it){
			output << "_" << *it << "\n";
		}
		for (it = starfile->data.begin() ; it != starfile->data.end(); ++it){
			output << *it << "\n";
		}
		output.close();
	} else {
		
	}
}

void getStarColumn(StarFile *starfile, std::string columnname, int *columnnumber){
	std::vector<std::string>::iterator 	it;
	int 								columncount = 0;
	
	for(it = starfile->columns.begin() ; it != starfile->columns.end(); ++it){
		if(*it == columnname){
			*columnnumber = columncount;
			return;
		}
		columncount++;
	}
}

void replaceStarColumn(StarFile *starfile, int row, int column, std::string value){
	std::size_t 	location;
	std::size_t 	endlocation;
	int 			currentcolumn;
	std::string		argstring;
	
	argstring = starfile->data[row];

	location = argstring.find("  ");
	while (location != std::string::npos) {
		argstring.replace (location, 2, " ");
		location = argstring.find("  ");
	}

	location = 0;
	
	while(argstring.c_str()[0] == ' '){
		argstring.erase (0, 1);
	}
	
	currentcolumn = 0;
	location = 0;
	
	while (currentcolumn < column){
		location = argstring.find(" ", location + 1);
		currentcolumn += 1;
	}
	
	endlocation = argstring.find(" ", location + 1);
	argstring.replace(location + 1, (endlocation - location) - 1, value);
	starfile->data[row] = argstring;
}
#endif
