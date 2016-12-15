package edu.kit.mmm.wrapper.simplePrime;

import java.io.*;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import org.apache.log4j.Logger;
import org.openmolgrid.process.ProcessWrapper;
import org.openmolgrid.wrapper.ApplicationWrapper;


public class SimplePrimeCalculation extends ApplicationWrapper 
{

  private static Logger log = Logger.getLogger(SimplePrimeCalculation.class.getName());

  private String workingDirectory = null;
  private String mpiCommand       = null;
  private int numberProcessors    = 1;


  /** Constructor */
  public SimplePrimeCalculation(String[] args) 
  {
    /** Get environment variables ********************************************/

    printArgs(args);

    // WorkingDir
    if (javaParameter.containsKey("--workingDirectory")) 
    {
      workingDirectory = javaParameter.get("--workingDirectory");
    } 
    else
    {
      log.error("Missing workingDirectory");  
      System.exit(1);
    }

    /**    
    // mpiCommand
    if (javaParameter.containsKey("--mpiCommand")) 
    {
      mpiCommand = javaParameter.get("--mpiCommand");
    } 
    else
    {
      log.error("Missing mpiCommand");  
      System.exit(1);
    }

    // numberProcessors
    if (javaParameter.containsKey("--numberProcessors")) 
    {
      numberProcessors = Integer.valueOf(javaParameter.get("--numberProcessors"));
    } 
    else 
    {
      log.error("Missing numberProcessors");  
      System.exit(1);
    }
    */    

  }


  /** Run application ********************************************************/

  @Override
  protected void preprocess() 
  {
    /** TODO Auto-generated method stub */
  }

  @Override
  protected void execute() throws Throwable 
  {
    log.info("Start run method");
    invokesimplePrime();
  }

  private void invokesimplePrime() throws IOException 
  {

    ProcessWrapper process = new ProcessWrapper(new File(workingDirectory));

    String run_application = "df";
    process.execute(run_application);
    log.info("Execute simplePrime: " + run_application);

    log.info(process.getStdout());
    log.info(process.getStderr());

    /**
    -------------------------------------------------------------------------------------
    Commands to run application on workstation
    String run_application = "/home/unicore/Inputs_and_executables/dlpoly/DLPOLY.Z";

    Command to launch DLPOLY.Z on CINECA PLX
    String run_application = "mpirun DLPOLY.Z";

    Command to launch DLPOLY.Z on CSC Louhi
    String run_application = mpiCommand  + " -n " + numberProcessors + " " + "DLPOLY.Z";
    -------------------------------------------------------------------------------------
    */

  }

  @Override
  protected void postprocess() 
  {
    /** TODO Auto-generated method stub */
  }

}

