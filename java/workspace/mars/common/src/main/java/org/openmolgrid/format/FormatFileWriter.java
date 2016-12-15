package org.openmolgrid.format;

import java.io.IOException;

/**
 * Interface for all FileWriter that write chemical formats
 * 
 * @author Stefan Bozic
 */
public interface FormatFileWriter {

	/**
	 * Writes the file
	 *  
	 * @throws IOException An exception that might occurs when writing the file
	 */
	public void writeFile() throws IOException;
	
}
