package edu.kit.gridbeans.simpleeovolassemble;

import java.net.URL;
import java.util.List;
import javax.xml.namespace.QName;
import com.intel.gpe.clients.api.Job;
import com.intel.gpe.clients.api.jsdl.gpe.GPEJob;
import com.intel.gpe.gridbeans.AbstractGridBean;
import com.intel.gpe.gridbeans.GPEConstants;
import com.intel.gpe.gridbeans.GridBeanException;
import com.intel.gpe.gridbeans.IGridBean;
import com.intel.gpe.gridbeans.parameters.GridBeanParameter;
import com.intel.gpe.gridbeans.parameters.GridBeanParameterType;
import com.intel.gpe.gridbeans.parameters.IGridBeanParameterValue;
import com.intel.gpe.gridbeans.parameters.InputFileParameterValue;
import com.intel.gpe.gridbeans.parameters.OutputFileParameterValue;
import com.intel.gpe.gridbeans.parameters.FileSetParameterValue;
import com.intel.gpe.gridbeans.parameters.processing.ParameterUtils;

public class SimpleEoVolAssembleGridBean extends AbstractGridBean implements IGridBean 
{

    /**
     * Add UID
     */
    private static final long serialVersionUID = -1701544599641820475L;

    private final static String QNAME_PREFIX = "SimpleEoVolAssemble";
  
    private static String APPLICATION_NAME = "SimpleEoVolAssemble";
    
    private static String APPLICATION_VERSION = "2.0";
  
    public final static String NAMESPACE = "http://edu.kit/gridbeans/SimpleEoVolAssemble";

    /** Declare input files */
    public final static QName INPUT_FILE = new QName(NAMESPACE, "INPUT_FILE", QNAME_PREFIX);
    public final static QName FILESET_IN   = new QName(NAMESPACE, "FILESET_IN",   QNAME_PREFIX);

    /** Declare output files */
    public final static QName OUTPUT_FILE = new QName(NAMESPACE, "OUTPUT_FILE", QNAME_PREFIX);
    public final static QName FILESET_OUT   = new QName(NAMESPACE, "FILESET_OUT",   QNAME_PREFIX);

    /**
     * QName
     */
    public final static QName TITLE   = new QName(NAMESPACE, "TITLE"  , QNAME_PREFIX);
    public final static QName STK     = new QName(NAMESPACE, "STK"    , QNAME_PREFIX);
    public final static QName NPART   = new QName(NAMESPACE, "NPART"  , QNAME_PREFIX);
    public final static QName NSTATES = new QName(NAMESPACE, "NSTATES", QNAME_PREFIX);
    public final static QName MSK     = new QName(NAMESPACE, "MSK"    , QNAME_PREFIX);
    public final static QName MW      = new QName(NAMESPACE, "MW"     , QNAME_PREFIX);
    public final static QName NTHR    = new QName(NAMESPACE, "NTHR"   , QNAME_PREFIX);
    public final static QName LPSTOP  = new QName(NAMESPACE, "LPSTOP" , QNAME_PREFIX);
    public final static QName INNER   = new QName(NAMESPACE, "INNER"  , QNAME_PREFIX);
    public final static QName WIDTH   = new QName(NAMESPACE, "WIDTH"  , QNAME_PREFIX);

    /** Constructor */
    public SimpleEoVolAssembleGridBean() 
    {
        set(JOBNAME, APPLICATION_NAME);
        createInputEnvironmentVariables();
        createInputFiles();
        createOutputFiles();
    }
  
    /**
     * create input environment variables (appears in the Variables tab of the GridBean)
     * --workingDirectory=$WORKING_DIR \
     *--stk=$STK \                      #stk=<ptcls.ext>
     *--npart=$NPART \                  #npart=<number of partitions to assemble>
     *--nstates=$NSTATES \              #nstates=<nr of states> 
     *--msk=$MSK \                      #msk=<mask radius(in pixels)> smpd=<sampling distance(in A)>
     *--mw=$MW \                        #[mw=<molecular weight(in kD)>]
     *--lpstop=$LPSTOP \                #[lpstop=<stay at this low-pass limit(in A)>]
     *--inner=$INNER \                  #[inner=<inner mask radius(in pixels)>]
     *--width=$WIDTH \                  #[width=<pixels falloff inner mask{10}>]'
     *--nthr=$NTHR                      #[nthr=<nr of openMP threads{1}>]
     *
     */
    private void createInputEnvironmentVariables() {
        QName[] envVariables = { TITLE,
                                 STK,
                                 NPART,
                                 NSTATES,
                                 MSK,
                                 MW,
                                 NTHR,
                                 LPSTOP,
                                 INNER,
                                 WIDTH };
        String[] initialValues = { "Title",
                                   "sumstack.mrc",
                                   "1",
                                   "1",
                                   "0.0",
                                   "0.0",
                                   "1",
                                   "7.0",
                                   "0.0",
                                   "10.0" };
        getInputParameters().addAll(ParameterUtils.createEnvParameters(envVariables));
        
        List<IGridBeanParameterValue> values = ParameterUtils.createEnvParameterValues(envVariables, initialValues);

        // initialize input environment variables
        for (int i = 0; i < initialValues.length; i++) {
            set(envVariables[i], values.get(i));
        }

    }
  
  /** Note: in createInputFiles and createOutputFiles, the InputFileParameterValues and
      OutputFileParameterValues are the specific names of the actual files as they appear
      in the working directory.
  */

  private void createInputFiles() 
  {
//    set(INPUT_FILE, new InputFileParameterValue("Input_File"));
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
  public URL getIconURL() {return getIconURL(SimpleEoVolAssembleGridBean.class, "images/StarFish_single_cropped.png");}


}
