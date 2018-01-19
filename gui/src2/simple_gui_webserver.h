// open include guard

#ifndef __WEBSERVER_H_INCLUDED__
#define __WEBSERVER_H_INCLUDED__

//=================================


// included dependencies

#include <string>
#include <map>
extern "C" {
    #include "../ext/mongoose/mongoose.h"
}

//=================================


// struct declaration

struct mg_serve_http_opts 	http_server_opts;

struct WebResponse {
	std::string							user;
	std::map<std::string, std::string>	requestmap;
	std::map<std::string, std::string>	responsemap;
	unsigned char* 						jpeg;
	unsigned long						jpegsize;
	int									jpegx;
	int									jpegy;
};


//=================================


// included internal dependencies

#include "simple_gui_environment.h"

//=================================


// send response

void processRequest(struct http_message* message, struct WebResponse* webresponse) {
	
	std::string		query_string;
	std::string		key;
	std::string		value;
	int 	equals;
	int 	ampersand;
	
	query_string.resize(message->query_string.len);
	memcpy(&query_string[0], message->query_string.p, message->query_string.len);
	
	ampersand = -1;
	equals = -1;

	while (true){
		equals = query_string.find("=", ampersand + 1);
		key = query_string.substr(ampersand + 1, equals - ampersand - 1);
		ampersand = query_string.find("&", equals + 1);
		value = query_string.substr(equals + 1, ampersand - equals - 1);
		if(key.size() > 0 && value.size() > 0){
			webresponse->requestmap[key] = value;
		}
		if(ampersand == std::string::npos){
			break;
		}
	}
	query_string.clear();
}
 
//=================================


// send response

void handleRequest(struct WebResponse* webresponse) {
	
	if (webresponse->requestmap.find("function") != webresponse->requestmap.end()) {
		log("Request function is " + webresponse->requestmap["function"]);
		if(webresponse->requestmap["function"] == "viewfile"){
			
		}else if (webresponse->requestmap["function"] == "newproject"){
			
		}else{
			log("Request function is invalid");
			webresponse->responsemap["error"] = "Invalid function in request";
		}
	} else {
		log("Request is missing a function");
		webresponse->responsemap["error"] = "Function missing from request";
	}
}
 
//=================================


// send response

void webResponse(struct mg_connection* http_connection, struct WebResponse* webresponse) {
	
	std::string JSONstring;
	
	if(sizeof(webresponse->jpeg) > 0){
		mg_send_head(http_connection, 200, webresponse->jpegsize, "Content-Type: image/jpeg");
		mg_send(http_connection, webresponse->jpeg, webresponse->jpegsize);
		delete [] webresponse->jpeg;
	} else {
		JSONstring = "joe";
		mg_send_head(http_connection, 200, JSONstring.length(), "Content-Type: application/json");
		mg_send(http_connection, JSONstring.c_str(), JSONstring.length());
	}

}
 
//=================================


//  event handler

void webEventHandler (struct mg_connection* http_connection, int event, void *eventdata) {
	
	struct http_message*   		message;
	struct WebResponse*			webresponse;
	
	message = (struct http_message *) eventdata;
	
	switch (event) {
		case MG_EV_HTTP_REQUEST: {
			if (mg_vcmp(&message->uri, "/handler") == 0) {
				log("Received HTTP request for handler");
				webresponse = new WebResponse();
				processRequest(message, webresponse);
				handleRequest(webresponse);
				webResponse(http_connection, webresponse);
				delete webresponse;
			} else {
				log("Received HTTP request for file from document root");
				mg_serve_http(http_connection, message, http_server_opts);
			}
		}
		default: {
			break;
		}
	}
}

//=================================


//  check file or folder exists

void webServer(){
	
	struct mg_mgr				http_manager;
	struct mg_connection*		http_connection;
	struct mg_serve_http_opts 	http_server_opts;
	sig_atomic_t				http_signal_received;
	
	if(chdir(environment.webpath.c_str()) == 0){
		log("Moved to web directory");
	}else{
		error("Failed to move to web directory");
    }
        
	mg_mgr_init(&http_manager, NULL);
		
	http_connection = mg_bind(&http_manager, std::to_string(environment.port).c_str(), webEventHandler); // cumbersome
	
	if (http_connection == NULL) {
		error("Failed to start web server");
	}else{
		log("Web server started");
	}
	
	mg_set_protocol_http_websocket(http_connection);	

	http_server_opts.document_root = ".";
	
	http_server_opts.enable_directory_listing = "no";
	
	if (environment.multiuser) {
		http_server_opts.per_directory_auth_file = ".htpasswd";
		http_server_opts.auth_domain = "simple";
	}
	
	http_signal_received = 0;
	
	while (http_signal_received == 0) {
		mg_mgr_poll(&http_manager, 200);
	}

	mg_mgr_free(&http_manager);
}

//=================================


// close include guard

#endif

//=================================
