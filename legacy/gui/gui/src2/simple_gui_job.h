// open include guard

#ifndef __JOB_H_INCLUDED__
#define __JOB_H_INCLUDED__

//=================================


// included external dependencies

#include <string>  
#include <sys/resource.h>

//=================================


// included internal dependencies

#include "simple_gui_sqlite.h" 
#include "simple_gui_logging.h"

//=================================


//  Fork job daemon from master and disconnect

void forkSimpleJobDaemon(){  
	
    pid_t				pid;
    pid_t				sid;   
	struct rlimit* 		fds;
	int 				fdscount;
	int 				r;
	struct stat statbuf;
	pid = fork();
    
    if (pid > 0){   
        return;
    }  
    
    if (pid < 0){     
        log("Failed to fork process");
        return;  
    }     
  
    sid = setsid();  
    
    if (sid < 0){  
        log("Failed to obtain new session id");
        return;  
    }  
  
    if ((chdir("/tmp")) < 0){  //make filesystem function
        exit(EXIT_FAILURE);  
    }  
    
	getrlimit(RLIMIT_NOFILE, fds); // get number possible open handles
	
	for(fdscount = 1; fdscount < fds->rlim_max; fdscount++){
		fstat(fdscount, &statbuf);
		if(S_ISSOCK(statbuf.st_mode)){
			close(fdscount); // Close all socket handles
		}
	}
}  

//=================================


// close include guard

#endif

//=================================
