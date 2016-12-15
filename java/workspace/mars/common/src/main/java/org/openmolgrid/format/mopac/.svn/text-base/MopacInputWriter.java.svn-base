package org.openmolgrid.format.mopac;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.openmolgrid.format.FormatFileWriter;
import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.util.CFormatter;

/**
 * Writes a MOPAC input file for a given {@link CStructure}.
 * 
 * @author Stefan Bozic
 */
public class MopacInputWriter implements FormatFileWriter {

	/** The structure to write as MOPAC input file */
	private CStructure structure;
	
	/** The name of the file to create */
	private String fileName;
	
	/** The MOPAC keywords */
	private String keywords;
	
	/** Geometry optimization flag */
	private boolean geometryOptimization;
	
	/**
	 * Constructor
	 * 
	 * @param structure The structure to write into the input file 
	 * @param fileName The name of the file to create
	 * @param keywords The MOPAC keywords
	 * @param geometryOptimization indicates if MOPAC should perform a geometry optimization
	 */
	public MopacInputWriter(CStructure structure, String fileName, String keywords, boolean geometryOptimization) {
		this.structure = structure;
		this.fileName = fileName;
		this.keywords = keywords;
		this.geometryOptimization = geometryOptimization;
	}
	
	/**
	 * Creates a MOPAC input file based on this structure
	 * 
	 * @throws IOException An exception that might occurs during file writing 
	 */
	public void writeFile() throws IOException {
		String optimizationFlag = geometryOptimization ? "1" : "0";
		
		PrintWriter outWriter = null;
		File f = new File(fileName);
		
		outWriter = new PrintWriter(new FileWriter(f));
		
		outWriter.println(keywords);
		outWriter.println();
		outWriter.println();
		for (int i = 0; i < structure.getAtomCount(); i++) {
			CAtom a = (CAtom) structure.getAtoms().get(i);
			StringBuilder line = new StringBuilder();
			line.append(CFormatter.formatStringLeft(a.getElement(), 2));
			line.append(CFormatter.formatDouble(a.getX(), 11, "0.00000"));
			line.append(" " + optimizationFlag);
			line.append(CFormatter.formatDouble(a.getY(), 11, "0.00000"));
			line.append(" " + optimizationFlag);
			line.append(CFormatter.formatDouble(a.getZ(), 11, "0.00000"));
			line.append(" " + optimizationFlag);
			outWriter.println(line.toString());
		}
		outWriter.close();
	}
}
