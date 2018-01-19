// open include guard

#ifndef __JSONREQUEST_H_INCLUDED__
#define __JSONREQUEST_H_INCLUDED__

//=================================


// included dependencies

#include <map>

extern "C" {
    #include "../ext/mongoose/mongoose.h"
}

//=================================



// class declaration

class JSONRequest {
	
	std::string							user;
	std::map<std::string, std::string>	requestmap;

	public:
		
		JSONRequest(struct http_message* message) {
			
		}
};

// class declaration

class JSONResponse {
	
	std::map<std::string, std::string>	responsemap;

	public:
		
		JSONResponse() {
			
		}
};

//=================================



// close include guard

#endif

//=================================
