#include <cstdlib>
#include <unistd.h>


#include "simple_gui_environment.h"
#include "simple_gui_logging.h"
#include "simple_gui_webserver.h"

Environment environment;

int main(int argc, char* argv[]) {
    environmentSetup(argc, argv);
    environmentCheck();
    webServer();
	return 0;
}
