package edu.kit.mmm.wrapper.simplePreProc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.log4j.Logger;
import org.openmolgrid.process.ProcessWrapper;
import org.openmolgrid.wrapper.ApplicationWrapper;


public class simplePreProcCalculation extends ApplicationWrapper {
  private static Logger log = Logger.getLogger(simplePreProcCalculation.class.getName());

  private String workingDirectory = null;
  private String mpiCommand       = null;
  private int numberProcessors    = 1;

  /** start the parameters for the AimplePreProc */
  private SimplePreProcParameters parameters;
  
  /** command file name */
  private final static String INPUT_SIMPLE_STRIDATA_INPUT_PM = "simple_StripData_input.pm";

  /** Constructor */
  public simplePreProcCalculation(String[] args) 
  {
    /** Get environment variables ********************************************/
      parameters = new SimplePreProcParameters();
      printArgs(args);
      parameters.mapParameter(javaParameter);

      // WorkingDir
      if (javaParameter.containsKey("--workingDirectory")) {
          workingDirectory = javaParameter.get("--workingDirectory");
      } else {
          log.error("Missing workingDirectory");  
          System.exit(1);
      }

    /**    
    // mpiCommand
    if (javaParameter.containsKey("--mpiCommand"))  {
      mpiCommand = javaParameter.get("--mpiCommand");
    } else {
      log.error("Missing mpiCommand");  
      System.exit(1);
    }

    // numberProcessors
    if (javaParameter.containsKey("--numberProcessors")) {
      numberProcessors = Integer.valueOf(javaParameter.get("--numberProcessors"));
    } else {
      log.error("Missing numberProcessors");  
      System.exit(1); }
    */    

  }

  /** creating the required inut files to run the preprocess */
  private void createSimpleStripDataInput_pmFile() throws IOException {
      log.info("Creating the simple_StripData_input.pm file");
      
      //file description      
      String path = parameters.getWorkingDirectory() != null ? parameters.getWorkingDirectory() + "/"
              + INPUT_SIMPLE_STRIDATA_INPUT_PM : INPUT_SIMPLE_STRIDATA_INPUT_PM;
      File commandFile = new File(path);

      FileWriter fw = null;
      BufferedWriter bw = null;

      try {
          fw = new FileWriter(commandFile);
          bw = new BufferedWriter(fw);

          bw.append("#!/usr/bin/perl\n");
          bw.append("################################################################################\n");
          bw.append("# Input variables for the compilation driving perl script                      #\n");
          bw.append("#                                                                              #\n");
          bw.append("#       Package simple_user_input                                              #\n");
          bw.append("#                                                                              #\n");
          bw.append("# perl script that lets the user input the global variables that defines the   #\n");
          bw.append("# global variables for the entire SIMPLE Application. This is the only place   #\n");
          bw.append("# where the user needs to intervene in the compilation step of the library.    #\n");
          bw.append("################################################################################\n");
          bw.append("#\n");
          bw.append("package simple_StripData_input;  # saved as simple_StripData_input.pm\n");
          bw.append("use lib './';\n");
          bw.append("use Exporter;\n");
          bw.append("use warnings;\n");
          bw.append("use Config;\n");
          bw.newLine();
          bw.append("@ISA = ('Exporter');\n");
          bw.append("@EXPORT = qw($NFRAMES $NCONFIG $SIMPLE_DATA_PATH $TARGET_FILE $CREATE_DIR $FILE_PRE_HANDLER $FILE_ORGANISE_DIR $DELETE_DIR $UNBLUR $CTFFIND $UNBLUR_DIR $CTFFIND_DIR $FILE_PST_HANDLER $HAVE_BOXFILE $SPH_ABE $AMP_CON $SZE_PWR_SPC $VERBOSE $DATA_ROT);\n");

      } finally {
          if (bw != null) {
              bw.close();
          }
          if (fw != null) {
              fw.close();
          }
      }

  }
  
  /** Run application ********************************************************/

  @Override
  protected void preprocess() 
  {
      try {
          createSimpleStripDataInput_pmFile();
      } catch (IOException e) {
          log.info("One of the input file variable has not been created properly:");
          log.info(INPUT_SIMPLE_STRIDATA_INPUT_PM/*getInputSimpleStridataInputPm()*/+"\n");
          e.printStackTrace();
      }
  }

  /* (non-Javadoc)
   * @see org.openmolgrid.wrapper.ApplicationWrapper#execute()
   */
  @Override
  protected void execute() throws Throwable 
  {
      log.info("Starting the run() method from the execute() method");
      run();
      log.info("The run() method from the execute() method has been executed");
  }
  /**
   * 
   */
  private void run() {
      log.info("Start of the SimplePreProc calculation");
      try {
          invokesimplePreProc();
          log.info("End of the SimplePreProc calculation");
      } catch (IOException e) {
          log.error("Error while running SimplePreProc job", e);
          e.printStackTrace();
      }
  }
  /**
   * @throws IOException
   */
  private void invokesimplePreProc() throws IOException 
  {
      log.info("#######################################");
      log.info("run SimplePreProc");
      log.info("#######################################");


      ProcessWrapper process = new ProcessWrapper(new File(workingDirectory));

      String run_application = "df";
      process.execute(run_application);
      log.info("Execute simplePreProc: " + run_application);

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

  //Setters and getters
//  /**
//   * @return the workingDirectory
//   */
//  public String getWorkingDirectory() {
//      return workingDirectory;
//  }
//  /**
//   * @param workingDirectory the workingDirectory to set
//   */
//  public void setWorkingDirectory(String workingDirectory) {
//      this.workingDirectory = workingDirectory;
//  }
//  /**
//   * @return the mpiCommand
//   */
//  public String getMpiCommand() {
//      return mpiCommand;
//  }
//  /**
//   * @param mpiCommand the mpiCommand to set
//   */
//  public void setMpiCommand(String mpiCommand) {
//      this.mpiCommand = mpiCommand;
//  }
//  /**
//   * @return the numberProcessors
//   */
//  public int getNumberProcessors() {
//      return numberProcessors;
//  }
//  /**
//   * @param numberProcessors the numberProcessors to set
//   */
//  public void setNumberProcessors(int numberProcessors) {
//      this.numberProcessors = numberProcessors;
//  }
//  /**
//   * @return the parameters
//   */
//  public simplePreProcParameters getParameters() {
//      return parameters;
//  }
//  /**
//   * @param parameters the parameters to set
//   */
//  public void setParameters(simplePreProcParameters parameters) {
//      this.parameters = parameters;
//  }
  /**
   * @return the inputSimpleStridataInputPm
   */
  public static String getInputSimpleStridataInputPm() {
      return INPUT_SIMPLE_STRIDATA_INPUT_PM;
  }

}

