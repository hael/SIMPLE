package edu.kit.gridbeans.simplepreproc;

import java.net.URL;
import java.util.List;

import javax.xml.namespace.QName;

import com.intel.gpe.clients.api.Job;
import com.intel.gpe.clients.api.jsdl.gpe.GPEJob;
import com.intel.gpe.gridbeans.AbstractGridBean;
import com.intel.gpe.gridbeans.GPEConstants;
import com.intel.gpe.gridbeans.GridBeanException;
import com.intel.gpe.gridbeans.IGridBean;
import com.intel.gpe.gridbeans.parameters.IGridBeanParameterValue;
import com.intel.gpe.gridbeans.parameters.processing.ParameterUtils;

public class SimplePreProcGridBean extends AbstractGridBean implements IGridBean 
{
    /**
     * Add UID
     */
    private static final long serialVersionUID = -730192220105271132L;

    private final static String QNAME_PREFIX = "SimplePreProc";
    /** Application name */
    private static String APPLICATION_NAME = "SimplePreProc";
    /** Application version */
    private static String APPLICATION_VERSION = "2.0";
    /** Project NAMESPACE */
    public final static String NAMESPACE = "http://edu.kit/gridbeans/SimplePreProc";

    /** Declare input files */
    public final static QName INPUT_FILE = new QName(NAMESPACE, "INPUT_FILE", QNAME_PREFIX);
    public final static QName FILESET_IN   = new QName(NAMESPACE, "FILESET_IN",   QNAME_PREFIX);

    /** Declare output files */
    public final static QName OUTPUT_FILE = new QName(NAMESPACE, "OUTPUT_FILE", QNAME_PREFIX);
    public final static QName FILESET_OUT   = new QName(NAMESPACE, "FILESET_OUT",   QNAME_PREFIX);

    /** QName */
    public final static QName TITLE             = new QName(NAMESPACE, "TITLE"             , QNAME_PREFIX);
    /** QNames for the Welcome panel */
    public final static QName ISTART            = new QName(NAMESPACE, "ISTART"            , QNAME_PREFIX); 
    /** QName for the PreProc*/
    
    public final static QName SIMPLE_DATA_PATH  = new QName(NAMESPACE, "SIMPLE_DATA_PATH"  , QNAME_PREFIX);
    public final static QName NFRAMES           = new QName(NAMESPACE, "NFRAMES"           , QNAME_PREFIX);
    public final static QName NCONFIG           = new QName(NAMESPACE, "NCONFIG"           , QNAME_PREFIX);
    public final static QName SPH_ABE           = new QName(NAMESPACE, "SPH_ABE"           , QNAME_PREFIX);
    public final static QName AMP_CON           = new QName(NAMESPACE, "AMP_CON"           , QNAME_PREFIX);
    public final static QName SZE_PWR_SPC       = new QName(NAMESPACE, "SZE_PWR_SPC"       , QNAME_PREFIX);
    public final static QName DATA_ROT          = new QName(NAMESPACE, "DATA_ROT"          , QNAME_PREFIX);
    public final static QName TARGET_FILE       = new QName(NAMESPACE, "TARGET_FILE"       , QNAME_PREFIX);
    public final static QName CREATE_DIR        = new QName(NAMESPACE, "CREATE_DIR"        , QNAME_PREFIX);
    public final static QName FILE_PRE_HANDLER  = new QName(NAMESPACE, "FILE_PRE_HANDLER"  , QNAME_PREFIX);
    public final static QName FILE_ORGANISE_DIR = new QName(NAMESPACE, "FILE_ORGANISE_DIR" , QNAME_PREFIX);
    public final static QName UNBLUR            = new QName(NAMESPACE, "UNBLUR"            , QNAME_PREFIX);
    public final static QName CTFFIND           = new QName(NAMESPACE, "CTFFIND"           , QNAME_PREFIX);
    public final static QName UNBLUR_DIR        = new QName(NAMESPACE, "UNBLUR_DIR"        , QNAME_PREFIX);
    public final static QName CTFFIND_DIR       = new QName(NAMESPACE, "CTFFIND_DIR"       , QNAME_PREFIX);
    public final static QName FILE_PST_HANDLER  = new QName(NAMESPACE, "FILE_PST_HANDLER"  , QNAME_PREFIX);
    public final static QName HAVE_BOXFILE      = new QName(NAMESPACE, "HAVE_BOXFILE"      , QNAME_PREFIX);
    public final static QName VERBOSE           = new QName(NAMESPACE, "VERBOSE"           , QNAME_PREFIX);
    public final static QName DELETE_DIR        = new QName(NAMESPACE, "DELETE_DIR"        , QNAME_PREFIX);

    /** Constructor */
    public SimplePreProcGridBean() 
    {
        set(JOBNAME, APPLICATION_NAME);
        createInputEnvironmentVariables();
        createInputFiles();
        createOutputFiles();
    }

    /**
     * create input environment variables (appears in the Variables tab of the GridBean)
     * --workingDirectory=$WORKING_DIR \
     * --simple_data_path=$SIMPLE_DATA_PATH \ #='path to data';
     * --nframes=$NFRAMES \           #= 7;
     * --nconfig=$NCONFIG \           #= 452; #not used for checking only, real nconfig extracted
     * --sph_abe=$SPH_ABE \           #= 2.7; #sphererical aberration
     * --amp_con=$AMP_CON \           #= 0.07;#Amplitude constrast
     * --sze_pwr_spc=$SZE_PWR_SPC \   #= 512; #Size of prower spectrum 
     * --data_rot=$DATA_ROT \         #= -90; #Rotation angle for the Size of prower spectrum 
     * --target_file=$TARGET_FILE \   #= 'target.txt';#Name of output file list
     * --create_dir=$CREATE_DIR \     #= 0;  #set to 1 if they not already created otherwise set to 0 
     * --file_pre_handler=$FILE_PRE_HANDLER \  # = "preprocess";
     * --file_organise_dir=$FILE_ORGANISE_DIR \ #= "or";
     * --unblur=$UNBLUR \             #= 1;
     * --ctffind=$CTFFIND \           #= 1;
     * --unblur_dir=$UNBLUR_DIR \     #= "/opt/Devel_tools/unblur_1.0.2/bin";
     * --ctffind_dir=$CTFFIND_DIR \   #= "/opt/Devel_tools/ctffind-4.0.16/build";
     * --file_pst_handler=$FILE_PST_HANDLER \  #= "extract_data";
     * --have_boxfile=$HAVE_BOXFILE \ #= "0";
     * --verbose=$VERBOSE \           #= "1";
     * --delete_dir=$DELETE_DIR       #= 0;
     *
     */
    private void createInputEnvironmentVariables() {
        QName[] envVariables = {
                TITLE,
                SIMPLE_DATA_PATH,
                NFRAMES,
                NCONFIG,
                SPH_ABE,
                AMP_CON,
                SZE_PWR_SPC,
                DATA_ROT,
                TARGET_FILE,
                CREATE_DIR,
                FILE_PRE_HANDLER,
                FILE_ORGANISE_DIR,
                UNBLUR,
                CTFFIND,
                UNBLUR_DIR,
                CTFFIND_DIR,
                FILE_PST_HANDLER,
                HAVE_BOXFILE,
                VERBOSE,
                DELETE_DIR };
        String[] initialValues = {
                "Title",
                "Path to data",
                "7",
                "452", //not used for checking only, real nconfig extracted
                "2.7", //#sphererical aberration
                "0.07",//#Amplitude constrast
                "512", //#Size of prower spectrum 
                "-90", //#Rotation angle for the Size of prower spectrum 
                "target.txt",//#Name of output file list
                "0", //#set to 1 if they not already created otherwise set to 0 
                "preprocess",
                "or",
                "0",
                "0",
                "/opt/Devel_tools/unblur_1.0.2/bin",
                "/opt/Devel_tools/ctffind-4.0.16/build",
                "extract_data",
                "0",
                "1",
        "0"};

        getInputParameters().addAll(ParameterUtils.createEnvParameters(envVariables));

        List<IGridBeanParameterValue> values = ParameterUtils.createEnvParameterValues(envVariables, initialValues);

        // initialize input environment variables
        for (int i = 0; i < initialValues.length; i++) {
            set(envVariables[i], values.get(i));
        }

        //setting up the input files for the files in tab in
        //      getInputParameters().addAll(ParameterUtils.createFileParameters(new QName[] { INPUT_XYZ_FILE }, true));
        //      set(INPUT_XYZ_FILE, new InputFileParameterValue("posinp.xyz"));

    }

  /** Note: in createInputFiles and createOutputFiles, the InputFileParameterValues and
      OutputFileParameterValues are the specific names of the actual files as they appear
      in the working directory.
  */

  private void createInputFiles() 
  {
//    set(INPUT_FILE, new InputFileParameterValue("simple_StripData_input.pm"));
//    
//    QName[] qnameArray = new QName[]{INPUT_FILE};
//    getInputParameters().addAll(ParameterUtils.createFileParameters(qnameArray, true));
//
//    set(FILESET_IN, new FileSetParameterValue());
//	getInputParameters().add(new GridBeanParameter(FILESET_IN, GridBeanParameterType.FILE_SET, true));
  }

  private void createOutputFiles() 
  {
//    set(OUTPUT_FILE, new OutputFileParameterValue("Output_File"));
//    QName[] qnameArray = new QName[]{OUTPUT_FILE};
//    getOutputParameters().addAll(ParameterUtils.createFileParameters(qnameArray, false));
//
//    set(FILESET_OUT, new FileSetParameterValue());
//    getOutputParameters().add(new GridBeanParameter(FILESET_OUT, GridBeanParameterType.FILE_SET, false));
  }


  /** Standard code for setting up job definition */
  public void setupJobDefinition(Job job) throws GridBeanException 
  {
    super.setupJobDefinition(job);
    if (job instanceof GPEJob) 
    {
      GPEJob gpeJob = (GPEJob) job;
      gpeJob.setApplicationName(APPLICATION_NAME);
      gpeJob.setApplicationVersion(APPLICATION_VERSION);
      gpeJob.setWorkingDirectory(GPEConstants.JobManagement.TEMPORARY_DIR_NAME);
    }
    else
    {
      throw new GridBeanException("Unsupported job class: " + job.getClass().getName());
    }
  }


  public String getName() {return APPLICATION_NAME;}
  public URL getIconURL() {return getIconURL(SimplePreProcGridBean.class, "images/StarFish_helix_cropped.png");}


}
