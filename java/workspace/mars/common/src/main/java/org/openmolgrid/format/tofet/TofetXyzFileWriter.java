package org.openmolgrid.format.tofet;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.apache.log4j.Logger;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.CSystem;
import org.openmolgrid.util.CenterOfMass;
import org.openmolgrid.util.CenterOfMassCalculator;


/**
 * Class that writes a ToFeT XYZ file based on a {@link CSystem}
 * 
 * @author Stefan Bozic
 */
public class TofetXyzFileWriter {

	/** Logger */
	private static Logger log = Logger.getLogger(TofetXyzFileWriter.class.getName());
	
	/**
	 * Creates a mopac input file based on a given {@link CSystem}
	 * 
	 * @param system the {@link CSystem} to create an xyz file for the sites 
	 * @param fileName the path for the file
	 */
	public static void writeFile(CSystem system, String fileName) {
		
		log.info("Write ToFeT xyz file to " + fileName);
		
		PrintWriter outWriter = null;
		File f = new File(fileName);
		try {
			outWriter = new PrintWriter(new FileWriter(f));

			for (CStructure structure : system.getStructures()) {
				CenterOfMass com = CenterOfMassCalculator.calculatesCenterOfMass(structure);
				outWriter.println(com.getX() + " " + com.getY() + " " + com.getZ() + " " + "-");
			}

		} catch (IOException e) {
			throw new RuntimeException(e.getMessage());
		} finally {
			outWriter.close();
		}
	}

}
