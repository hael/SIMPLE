package org.openmolgrid.process;

/***
 * Interface for classes that want to be informed if a specific String has been read by a StreamReader. 
 * 
 * @author Stefan Bozic
 */
public interface StreamReaderHandler {

	/***
	 * Is called by a StreamReader if a specific {@link String} has been read. 
	 */
	public void stringRead();
}
