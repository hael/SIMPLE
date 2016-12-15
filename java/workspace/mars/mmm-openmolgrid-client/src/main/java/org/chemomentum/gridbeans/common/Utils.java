package org.chemomentum.gridbeans.common;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.openmolgrid.format.slf.StructureWriter;
import org.openmolgrid.model.CStructure;


/**
 * Common utilities for GridBeans.
 *
 * @author Sulev Sild
 */
public class Utils {

	/**
	 * Writes the given list of structures to a temporary file.
	 *
	 * @param sl - ArrayList with CStructure objects
	 * @return filename for the temporary file created
	 * @throws IOException on I/O error
	 */
	public static String writeTempStructureList(ArrayList<CStructure> sl)
		throws IOException {

		String fname = File.createTempFile("structures", ".slf").getAbsolutePath();

		StructureWriter sw = new StructureWriter(fname);
		for (int i=0; i<sl.size(); i++) {
			CStructure str=sl.get(i);
			str.setId(Integer.toString(i+1));
			sw.addStructure(str);
		}
		sw.close();

		return fname;
	}

}
