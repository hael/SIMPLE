package org.chemomentum.gridbeans.common;

import javax.xml.namespace.QName;

import com.intel.gpe.clients.api.Application;
import com.intel.gpe.clients.api.Job;
import com.intel.gpe.clients.api.jsdl.gpe.GPEJob;
import com.intel.gpe.clients.common.ApplicationImpl;
import com.intel.gpe.gridbeans.AbstractGridBean;
import com.intel.gpe.gridbeans.GPEConstants;
import com.intel.gpe.gridbeans.GridBeanException;
import com.intel.gpe.gridbeans.apps.FileParameterMetaDataImpl;
import com.intel.gpe.gridbeans.parameters.GridBeanParameter;
import com.intel.gpe.gridbeans.parameters.GridBeanParameterType;
import com.intel.gpe.gridbeans.parameters.InputFileParameterValue;
import com.intel.gpe.gridbeans.parameters.OutputFileParameterValue;

/**
 * Base class for gridbeans related to QSAR modeling. It takes care of setting
 * up a JSDL job object and setting up file stagein/stageout.
 *  
 * @author Sulev Sild
 */
public class GridBeanHelper extends AbstractGridBean {

	private static final long serialVersionUID = -550510884846085658L;

	private static final String prefix = "c9mapp";
	private static final String NS = "http://chemomentum.org/gridbeans";
	public static final QName APPLICATION = new QName(NS, "APPLICATION", prefix);
	public static final QName INPUT   = new QName(NS, "INPUT", prefix);
	public static final QName OUTPUT  = new QName(NS, "OUTPUT", prefix);
	public static final QName ERROR   = new QName(NS, "ERRLOG", prefix);
	
	private String gbName;
	
	/**
	 * Constructor. 
	 * 
	 * @param gbName - GridBean name
	 * @param application - application name
	 * @param version - application version
	 */
	public GridBeanHelper(String gbName, String application, String version) {
		this.gbName = gbName;
		
		set(APPLICATION, new ApplicationImpl(application, version));
		set(JOBNAME, "noname");
		InputFileParameterValue input = new InputFileParameterValue("input.slf");
		input.setMetaData(new FileParameterMetaDataImpl("text/slf"));
		set(INPUT, input);
		OutputFileParameterValue output = new OutputFileParameterValue("output.slf");
		output.setMetaData(new FileParameterMetaDataImpl("text/slf"));
		set(OUTPUT, output);
		OutputFileParameterValue errorlog = new OutputFileParameterValue("errorlog");
		errorlog.setMetaData(new FileParameterMetaDataImpl("text/plain"));
		set(ERROR, errorlog);
		
		getInputParameters().add(new GridBeanParameter(INPUT, GridBeanParameterType.FILE));
		getOutputParameters().add(new GridBeanParameter(OUTPUT, GridBeanParameterType.FILE, false));
		getOutputParameters().add(new GridBeanParameter(ERROR, GridBeanParameterType.FILE, false));
	}

	/* (non-Javadoc)
	 * @see com.intel.gpe.gridbeans.IGridBean#getName()
	 */
	public String getName() {
		return gbName;
	}

	/**
	 * generates JSDL document from the data in the model.
	 * 
	 * @see com.intel.gpe.gridbeans.IGridBean#setupJobDefinition(com.intel.gpe.clients.api.Job)
	 */
	@Override
	public void setupJobDefinition(Job job) throws GridBeanException {
		super.setupJobDefinition(job);
		
		if (job instanceof GPEJob) {
			GPEJob j = (GPEJob) job;

			j.setId((String) get(JOBNAME));

			Application app = (Application)get(APPLICATION);
			j.setApplicationName(app.getName());
			j.setApplicationVersion(app.getApplicationVersion());
			j.setWorkingDirectory(GPEConstants.JobManagement.TEMPORARY_DIR_NAME);
			
			for (QName key: keySet()) {
                if (key.equals(JOBNAME) || key.equals(APPLICATION)) {
                    continue;
                }
                
                Object val = get(key);
                if (val == null) {
                	continue;
                } else if (val instanceof String) {
                	j.addField(key.getLocalPart(), (String)val);
                }
            }
		} else {
			throw new GridBeanException("Unsupported job class: "
					+ job.getClass().getName());
		}		
	}
	
	/** 
	 * Transfers values from the JSDL back into the model
	 * 
	 * @see com.intel.gpe.gridbeans.AbstractGridBean#parseJobDefinition(com.intel.gpe.clients.api.Job)
	 */
	@Override
	public void parseJobDefinition(Job job) throws GridBeanException {
		super.parseJobDefinition(job);
		if (job instanceof GPEJob) {
			GPEJob j = (GPEJob) job;
			for (QName key: keySet()) {
				if (get(key) instanceof String) {
					String val = j.getField(key.getLocalPart());
					if (val != null) {
						set(key, val);
					}
				}
			}
			
			Application app = new ApplicationImpl(
				j.getApplicationName(), j.getApplicationVersion());
			set(APPLICATION, app);
			set(JOBNAME, j.getId());
		}
	}
}
