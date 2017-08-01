/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 08th of April 2011. 15:01am
 *
 * Name:
 * StartLogs.cpp - basic definitions used in all modules.
 *
 * Description:
 * routine to start the loags on from the card.
 *******************************************************************************
 */
#include <stdio.h>
#include <stdlib.h>

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <oclUtils.h>
#include <CL/cl.h>
#include <shrUtils.h>

extern "C" void StartLogs (char*);

void StartLogs (char* filenamelog)
{

  // start the logs

  shrSetLogFileName(filenamelog);
  shrLog("%s Starting...\n\n"); 

}

extern "C" void startlogs_gpu_(char*) __attribute__((weak,alias("StartLogs")));

#endif /* OPENCL */
