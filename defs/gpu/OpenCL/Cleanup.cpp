/*******************************************************************************
 *     Author: Frederic D.R. Bonnet date: 28th of July 2013.13:19am
 *
 * Name:
 * Cleanup.cpp - basic definitions used in all modules.
 *
 * Description:
 * method to cleanup the ressources
 *******************************************************************************
 */

#if defined (OPENCL) /*preprossing for the OPENCL environment */

#include <oclUtils.h>
#include <CL/cl.h>
#include <shrUtils.h>

#include <simul.h>

extern "C" void Cleanup(int);

void Cleanup(int iExitCode)
{

  shrLog("\nStarting Cleanup...\n\n");

  shrLog("Ending...\n Exiting\n");
  exit (iExitCode);
}

#endif /* OPENCL */
