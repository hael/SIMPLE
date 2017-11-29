#ifndef SIMPLE_GUI_HTTP_H_INCLUDED
#define SIMPLE_GUI_HTTP_H_INCLUDED

#include "simple_gui_view.h"
#include "simple_gui_filesystem.h"
#include "simple_gui_job.h"
#include "simple_gui_project.h"

#include <string>
#include <cstring>
#include <iostream>
#include <cstdlib>

extern "C" {
    #include "ext/mongoose/mongoose.h"
}

const char* 					          http_port = "8084";
struct mg_serve_http_opts 	http_server_opts;
sig_atomic_t 				          http_signal_received = 0;
struct mg_mgr 				         http_manager;

void JSONHandler (struct mg_connection* http_connection, struct http_message* message) {
 char function[32];
 mg_get_http_var(&message->query_string, "function", function, sizeof(function));
 
 if (std::strcmp(function,"viewfile") == 0) {
  fileViewerInit (http_connection, message);
 } else if (std::strcmp(function,"listdir") == 0) {
	listDirectory (http_connection, message);
 } else if (std::strcmp(function,"savesel") == 0) {
	saveSelection (http_connection, message);
 } else if (std::strcmp(function,"externaljob") == 0) {
	externalJob (http_connection, message);
 } else if (std::strcmp(function,"newproject") == 0) {
	newProject (http_connection, message);
 } else if (std::strcmp(function,"2dview") == 0) {
	view2DInit (http_connection, message);
 } else if (std::strcmp(function,"pipelineview") == 0) {
	pipelineviewInit (http_connection, message);
 } else if (std::strcmp(function,"boxfiledata") == 0) {
	boxFileData (http_connection, message);
 }	else if (std::strcmp(function,"getprojects") == 0) {
	getProjects (http_connection, message);
 }	else if (std::strcmp(function,"getjobs") == 0) {
	getJobs (http_connection, message);
 } /*else if (std::strcmp(function,"addjob") == 0) {
	addJob (http_connection, message);
 }*/
  else if (std::strcmp(function,"particleviewer") == 0) {
	particleViewInit (http_connection, message);
 }else if (std::strcmp(function,"save2dselection") == 0) {
	save2DSelection (http_connection, message);
 }else if (std::strcmp(function,"save2dselectionparticles") == 0) {
	save2DSelectionParticles (http_connection, message);
 }
}

void JPEGHandler (struct mg_connection* http_connection, struct http_message* message) {
 char filename[1024];
 mg_get_http_var(&message->query_string, "filename", filename, sizeof(filename));
 
 if (std::strstr(filename,".mrc") != NULL) {
  MRCToJPEG (http_connection, message);
 }
}

void webEventHandler (struct mg_connection* http_connection, int event, void *eventdata) {
 struct http_message*   message;
 
 message = (struct http_message *) eventdata;

	switch (event) {
		case MG_EV_HTTP_REQUEST: {
   if (mg_vcmp(&message->uri, "/JSONhandler") == 0) {
    std::cout << "JSON HTTP Request Served" << std::endl;
    JSONHandler(http_connection, message);
   }else if (mg_vcmp(&message->uri, "/JPEGhandler") == 0) {
    std::cout << "JPEG HTTP Request Served" << std::endl;
    JPEGHandler(http_connection, message);
   } else {
    std::cout << "Static HTTP Request Served" << std::endl;
    mg_serve_http(http_connection, message, http_server_opts);
   }
  }
  default: {
   break;
  }
 }
}

void webServer(int &port){
 struct mg_mgr 				         http_manager;
 struct mg_connection*      http_connection;
 
	mg_mgr_init(&http_manager, NULL);
	std::cout << "Starting webserver on port " << port << " ... ";
	http_connection = mg_bind(&http_manager, http_port, webEventHandler);
	if (http_connection == NULL) {
		std::cout << "Failed" << std::endl;
		exit(1);
	}else{
		std::cout << "Success" << std::endl;
	}
	mg_set_protocol_http_websocket(http_connection);	
	http_server_opts.document_root = "www";
	http_server_opts.enable_directory_listing = "no";
	while (http_signal_received == 0) {
		mg_mgr_poll(&http_manager, 200);
	}
	mg_mgr_free(&http_manager);
}

void returnWSMessage(struct mg_connection *nc, std::string errorstring){
	std::string returnstring = "message=";
	for(unsigned int i = 0; i < errorstring.length(); i++){
		if(errorstring.c_str()[i] == ' '){
			returnstring += "^";
		}else{
			returnstring += errorstring.c_str()[i];
		}
	}
	mg_printf_websocket_frame(nc, WEBSOCKET_OP_TEXT, returnstring.c_str(), returnstring.length());
}


#endif
