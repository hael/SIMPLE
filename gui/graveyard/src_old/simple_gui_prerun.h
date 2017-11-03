#ifndef SIMPLE_GUI_PRERUN_H_INCLUDED
#define SIMPLE_GUI_PRERUN_H_INCLUDED

void preRunUnblur (struct mg_connection *nc, std::string &messagedata, std::string &argscommand) {
	
	int 								columnid;
	std::string							starfilename;
	std::string							element;
	std::ofstream 						outfile;
	std::vector<std::string>::iterator	it;

	StarFile* starfile = new StarFile();

	getKeyValFromString(&starfilename, messagedata, "starfile");
	readStarFile(starfile, starfilename);
	outfile.open("movies.txt");
	if(outfile.is_open()){
		getStarColumn(starfile, "splMovie", &columnid);
		for(it = starfile->data.begin() ; it != starfile->data.end(); ++it){
			getElementFromString(&element, *it, columnid);
			outfile << element << "\n";
		}
		outfile.close();
		argscommand.append(" filetab=movies.txt");
	} else {
		return;
	}
	
	delete starfile;
	
}

void preRunUnblurCtffind (struct mg_connection *nc, std::string &messagedata, std::string &argscommand) {
	
	int 								columnid;
	std::string							starfilename;
	std::string							element;
	std::ofstream 						outfile;
	std::vector<std::string>::iterator	it;

	StarFile* starfile = new StarFile();
	
	getKeyValFromString(&starfilename, messagedata, "starfile");
	readStarFile(starfile, starfilename);
	outfile.open("movies.txt");
	if(outfile.is_open()){	
		getStarColumn(starfile, "splMovie", &columnid);	
		for(it = starfile->data.begin() ; it != starfile->data.end(); ++it){
			getElementFromString(&element, *it, columnid);
			outfile << element << "\n";
		}
		outfile.close();
		argscommand.append(" filetab=movies.txt");
	} else {
		return;
	}
	
	delete starfile;
	
}

void preRunCtffind (struct mg_connection *nc, std::string &messagedata, std::string &argscommand) {
	
	int 								columnid;
	std::string							starfilename;
	std::string							element;
	std::ofstream 						outfile;
	std::vector<std::string>::iterator	it;

	StarFile* starfile = new StarFile();

	getKeyValFromString(&starfilename, messagedata, "starfile");
	readStarFile(starfile, starfilename);
	outfile.open("micrographsforctf.txt");
	if(outfile.is_open()){
		getStarColumn(starfile, "splForCTFName", &columnid);
		for(it = starfile->data.begin() ; it != starfile->data.end(); ++it){
			getElementFromString(&element, *it, columnid);
			outfile << element << "\n";
		}
		outfile.close();
		argscommand.append(" filetab=micrographsforctf.txt");
	} else {
		return;
	}
	
	delete starfile;
	
}

void preRunPreproc (struct mg_connection *nc, std::string &messagedata, std::string &argscommand) {
	
	int				watcherpid;
	std::string		movieloc;
	
	getKeyValFromString(&movieloc, messagedata, "movieloc");
	symlink(movieloc.c_str(), "movies");
	mkdir("pipeline", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	argscommand.append(" dir_movies=movies dir_target=pipeline");
	watcherpid = fork();
	if (watcherpid == 0) {
		watchRunPreproc(nc, messagedata, 0);
	}
	
}

void jobPreRun (struct mg_connection *nc, std::string &messagedata, std::string &argscommand) {
	
	std::string		prg;
	
	getKeyValFromString(&prg, messagedata, "prg");
	
	if(prg == "unblur"){
		preRunUnblur(nc, messagedata, argscommand);
	}else if(prg == "unblur_ctffind"){
		preRunUnblurCtffind(nc, messagedata, argscommand);
	}else if(prg == "ctffind"){
		preRunCtffind(nc, messagedata, argscommand);
	}else if(prg == "preproc"){
		preRunPreproc(nc, messagedata, argscommand);
	}
	
}



#endif
