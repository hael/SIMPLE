#ifndef SIMPLE_GUI_HTTP_H_INCLUDED
#define SIMPLE_GUI_HTTP_H_INCLUDED

static void returnWSMessage(struct mg_connection *nc, std::string errorstring){
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
