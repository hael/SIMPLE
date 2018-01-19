#ifndef SIMPLE_GUI_STARDATA_H_INCLUDED
#define SIMPLE_GUI_STARDATA_H_INCLUDED

struct StarFile {
	std::vector<std::string> columns;
	std::vector<std::string> data;
};

void readStarFile(StarFile *starfile, std::string filename){
	
	std::ifstream 	input;
	std::string 	line;

	input.open(filename.c_str());
	
	if(input.is_open()){
		while(getline(input, line)){
			if(line.at(0) == '_'){
				line.erase(0,1);
				starfile->columns.push_back(line);
			} else if(line.at(0) == '#' || line.at(0) == ' '){
			} else if(line.find("data_") != std::string::npos || line.find("loop_") != std::string::npos) {
			} else {
				starfile->data.push_back(line);
			}
		}
		input.close();
	} else {
		std::cout << "no starfil" << std::endl;
	}
}

static void writeStarFile(StarFile *starfile, std::string filename){
	std::ofstream output;
	std::vector<std::string>::iterator it;
	
	output.open(filename.c_str());
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

#endif
