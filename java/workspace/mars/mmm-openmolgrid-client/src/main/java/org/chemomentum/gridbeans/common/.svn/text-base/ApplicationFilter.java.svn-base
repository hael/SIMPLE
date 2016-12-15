package org.chemomentum.gridbeans.common;

import java.util.ArrayList;
import java.util.List;

import com.intel.gpe.clients.api.Application;
import com.intel.gpe.clients.api.TargetSystemClient;
import com.intel.gpe.clients.api.exceptions.GPEException;

/**
 * This class provides a filter for selecting supported applications from the 
 * target system.
 * 
 * @author Sulev Sild
 */
public class ApplicationFilter {

	private String[] patterns;

	/**
	 * Constructor.
	 * 
	 * @param patterns - array with application name patterns
	 */
	public ApplicationFilter(String ... patterns) {
		this.patterns = patterns;
	}
	
	/**
	 * Filters out supported applications from the given target system.
	 * 
	 * @param tsc
	 * @return
	 */
	public List<Application> filter(TargetSystemClient tsc) {
		List<Application> result = new ArrayList<Application>();
		
		if (tsc != null) {
			List<Application> applications;
			try {
				applications = tsc.getApplications();
			} catch (GPEException e) {
				return result;
			}
			for (String pat: patterns) {
				for (Application a : applications) {
					if (pat.equals(a.getName())) {
						result.add(a);
					}
				}	
			}
		}
		
		return result;
	}
	
}
