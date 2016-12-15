package edu.kit.gridbeans.simple;

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

/**
 * GridBean for the Simple library
 * 
 * @author Frederic Bonnet
 * @date 12th of November 2015
 *
 */

public class SimpleGridBean extends AbstractGridBean implements IGridBean  {

	/**
	 * Add UID
	 */
	private static final long serialVersionUID = 5732367727306222357L;

  private final static String QNAME_PREFIX = "Simple";
  /** Application name */
  private static String APPLICATION_NAME = "Simple";
  /** Application version */ 
  private static String APPLICATION_VERSION = "2.0";
  /** Project NAMESPACE */
  public final static String NAMESPACE = "http://edu.kit/gridbeans/Simple";

  /** Declare input files */
  public final static QName INPUT_FILE = new QName(NAMESPACE, "INPUT_FILE", QNAME_PREFIX);
  public final static QName FILESET_IN   = new QName(NAMESPACE, "FILESET_IN", QNAME_PREFIX);

  /** Declare output files */
  public final static QName OUTPUT_FILE = new QName(NAMESPACE, "OUTPUT_FILE", QNAME_PREFIX);
  public final static QName FILESET_OUT   = new QName(NAMESPACE, "FILESET_OUT", QNAME_PREFIX);

  /**
   * QName
   */
  public final static QName TITLE = new QName(NAMESPACE, "TITLE", QNAME_PREFIX);
  //QNames for the Welcome panel
  public final static QName ISTART = new QName(NAMESPACE, "ISTART" , QNAME_PREFIX); 
  //QNames for the SimpleJPanel

  /** Constructor */
  public SimpleGridBean() 
  {
    set(JOBNAME, APPLICATION_NAME);
    createInputEnvironmentVariables();
    createInputFiles();
    createOutputFiles();
  }

  
	/**
	 * create input environment variables (appears in the Variables tab of the GridBean)
	 */
	private void createInputEnvironmentVariables() {
		QName[] envVariables = { TITLE };
		String[] initialValues = { "Title" };

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
    set(INPUT_FILE, new InputFileParameterValue("Input_File"));
    
    QName[] qnameArray = new QName[]{INPUT_FILE};
    getInputParameters().addAll(ParameterUtils.createFileParameters(qnameArray, true));

    set(FILESET_IN, new FileSetParameterValue());
	getInputParameters().add(new GridBeanParameter(FILESET_IN, GridBeanParameterType.FILE_SET, true));
  }

  private void createOutputFiles() 
  {
    set(OUTPUT_FILE, new OutputFileParameterValue("Output_File"));
    QName[] qnameArray = new QName[]{OUTPUT_FILE};
    getOutputParameters().addAll(ParameterUtils.createFileParameters(qnameArray, false));

    set(FILESET_OUT, new FileSetParameterValue());
    getOutputParameters().add(new GridBeanParameter(FILESET_OUT, GridBeanParameterType.FILE_SET, false));
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
  public URL getIconURL() {return getIconURL(SimpleGridBean.class, "images/StarFish_single_cropped.png");}


}
