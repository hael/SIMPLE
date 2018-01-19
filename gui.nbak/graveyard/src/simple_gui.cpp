#define KEY(X) #X
#define VALUE(X) KEY(X)
#define SIMPLE_DIR_VALUE VALUE(SIMPLE_DIR)

/*#include <stdio.h>
#include <jpeglib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <sqlite3.h> 
#include <sys/types.h>
#include <pwd.h>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <ctime>
#include <vector>
#include <signal.h>*/


#include "simple_gui_http.h"
/*#include "simple_gui_postrun.h"

#include "simple_gui_textdata.h"
#include "simple_gui_stardata.h"
#include "simple_gui_mrcdata.h"
#include "simple_gui_unidocdata.h"
#include "simple_gui_midrun.h"
#include "simple_gui_prerun.h"
#include "simple_gui_view.h"*/


int main(void) {
 int port = 8084;
 webServer(port);
 return 0;
}
