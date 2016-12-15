/**
 * 
 */
package org.openmolgrid.format.bigdft;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import org.openmolgrid.format.yaml.GeometryYamlParse;
import org.openmolgrid.model.COrbitals;

/**
 * @author Frederic Bonnet
 * 
 * Class to write the input files from either parsed or other input files.
 *
 */
public class CBigDFTInputWriter implements BigDFTInputWritter {

	/** Logger */
	private static Logger log = Logger.getLogger(CBigDFTInputWriter.class.getName());
	//the new written filename 
	private final static String NEW_INPUT_KPT_FILE_NAME = "input.kpt";
	private final static String NEW_INPUT_XYZ_FILE_NAME = "posinp_opti.xyz";
	//The different possible boundary conditions for the system
	private final static String FREE_BOUNDARY_CONDITION = "Free";
	private final static String SURFACE_BOUNDARY_CONDITION = "Surface";
	private final static String PERIODIC_BOUNDARY_CONDITION = "Periodic";

	private static final int X_INDEX = 0;
	private static final int Y_INDEX = 1;
	private static final int Z_INDEX = 2;
	
	/** Start parameters for BIGDFT GridBean */
	//private BigDFTParameters parameters;
	//local variables
	private PrintWriter pw;
	//Global variables
	private COrbitals corbitals;
	/**
	 * Simple constructor
	 */
	public CBigDFTInputWriter() {
	}
	/**
	 * constructor with PrintWriter 
	 * @param pw
	 */
	public CBigDFTInputWriter(PrintWriter pw) {
		this.pw = pw;
	}
	/**
	 * constructor to take in class orbitals
	 * @param corbitals
	 */
	public CBigDFTInputWriter(COrbitals corbitals) {
		this.corbitals = corbitals;
	}
	/**
	 * Method used in interface
	 * it includes:
	 * 	writeCMLString(String path, String cmlString) throws IOException
	 * 	kptInputWritter()
	 * 	dftInputWritter()
	 * 	xyzInputWritter()
	 */
	public void writeCMLString(String path, String cmlString) throws IOException {
		FileWriter fw = null;
		BufferedWriter bw = null;
		File testIOFile = new File(path);
		try {
			fw = new FileWriter(testIOFile);
			bw = new BufferedWriter(fw);
			
			bw.append(cmlString);
			
		} finally {
			if (bw != null) {
				bw.close();
			}
			if (fw != null) {
				fw.close();
			}
		}
	}
	/**
	 * writer method to write out the input.kpt file.
	 * @param pth
	 * @throws IOException
	 */
	public void kptInputWritter(String pth) throws IOException {
//		String path = parameters.getWorkingDirectory() != null ? parameters.getWorkingDirectory()+ "/"
//				+ NEW_INPUT_KPT_FILE_NAME : NEW_INPUT_KPT_FILE_NAME;
		//TODO: need to fix the path of the file output writter
		String path =
				"src" + File.separator + 
				"test" + File.separator +
				"resources" + File.separator +
				"yaml"
				!= null ?
						"src" + File.separator +
						"test" + File.separator +
						"resources" + File.separator +
						"yaml"+ File.separator +
						NEW_INPUT_KPT_FILE_NAME : NEW_INPUT_KPT_FILE_NAME;
		File commandFile = new File(path);
		
		FileWriter fw = null;
		BufferedWriter bw = null;
		
		try {
			fw = new FileWriter(commandFile);
			bw = new BufferedWriter(fw);
			
			bw.append("This is the new kpt file"); //TODO: line needs to be removed once settle
			
		} finally {
			if (bw != null) {
				bw.close();
			}
			if (fw != null) {
				fw.close();
			}
		}
	
	}
	public void dftInputWritter() {
	}
	/**
	 * writer method to write out the input_opti.xyz file.
	 * @param pth
	 * @throws IOException
	 */
	public void xyzInputWritter(String path) throws IOException {

		File commandFile = new File(path);
		
		FileWriter fw = null;
		BufferedWriter bw = null;
		
		try {
			fw = new FileWriter(commandFile);
			bw = new BufferedWriter(fw);
			
			bw.append("This is the new xyz opti file"); //TODO: line needs to be removed once settle
			
		} finally {
			if (bw != null) {
				bw.close();
			}
			if (fw != null) {
				fw.close();
			}
		}
		
	}
	/**
	 * method to write the optimized posinp.xyz file
	 * 
	 * @param path
	 * @param XcoordArrayList
	 * @param YcoordArrayList
	 * @param ZcoordArrayList
	 * @param AtomNamecoordArrayList
	 * @param rigidShiftList
	 * @param poisson_kernelGeoYamlParse
	 * @param simulation_domainGeoYamlParse
	 * @throws IOException
	 */
	public void xyzInputWritter(String path,
			ArrayList<Double> XcoordArrayList, ArrayList<Double> YcoordArrayList, ArrayList<Double> ZcoordArrayList,
			ArrayList<String> AtomNamecoordArrayList, List<?> rigidShiftList,
			GeometryYamlParse poisson_kernelGeoYamlParse,
			GeometryYamlParse simulation_domainGeoYamlParse) throws IOException {

		poisson_kernelGeoYamlParse.ParsePoissonKernel();
		simulation_domainGeoYamlParse.ParseSimulationdomain();
		
		File commandFile = new File(path);
		
		FileWriter fw = null;
		BufferedWriter bw = null;
		
		try {
			fw = new FileWriter(commandFile);
			bw = new BufferedWriter(fw);

			bw.append(XcoordArrayList.size()+" "+"atomic");
			bw.newLine();
			String iBC = poisson_kernelGeoYamlParse.getBoundaryConditions();
			if (iBC.equalsIgnoreCase(FREE_BOUNDARY_CONDITION)) {
				bw.append(iBC.toLowerCase());
				bw.newLine();
			} else if (iBC.equalsIgnoreCase(SURFACE_BOUNDARY_CONDITION)) {
				bw.append(iBC.toLowerCase()+" "+
						simulation_domainGeoYamlParse.getSimulationDomain_AU_List().get(X_INDEX)+" "+
						"0 "+
						simulation_domainGeoYamlParse.getSimulationDomain_AU_List().get(Z_INDEX));
				bw.newLine();
			} else if (iBC.equalsIgnoreCase(PERIODIC_BOUNDARY_CONDITION)) {
				bw.append(iBC.toLowerCase()+" "+
						simulation_domainGeoYamlParse.getSimulationDomain_AU_List().get(X_INDEX)+" "+
						simulation_domainGeoYamlParse.getSimulationDomain_AU_List().get(Y_INDEX)+" "+
						simulation_domainGeoYamlParse.getSimulationDomain_AU_List().get(Z_INDEX));
				bw.newLine();				
			}
			for (int i = 0; i < XcoordArrayList.size(); i++) {
				bw.append(AtomNamecoordArrayList.get(i)+" "+XcoordArrayList.get(i)+" "+YcoordArrayList.get(i)+" "+ZcoordArrayList.get(i));
				bw.newLine();
			}
		} finally {
			if (bw != null) {
				bw.close();
			}
			if (fw != null) {
				fw.close();
			}
		}
	}
	/**
	 * Helper methods to print info
	 * @param orbitals
	 * @return sw.toString():String
	 */
	public static String orbitalsAsString(COrbitals orbitals) {
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		pw.println("#The new input.kpt file");
		pw.println("<?xml version=\"1.0\"?>");
		new CBigDFTInputWriter(pw).writeOrbitals(orbitals);
		return sw.toString();
	}
	/**
	 * method to write out the orbitals in a xml format
	 * @param orbitals
	 */
	private void writeOrbitals(COrbitals orbitals) {
		int i;
		pw.println("<orbitals id=\""+"some id"+"\">");
		pw.println("   <eigenvaluesArray>");
		for ( i = 0 ; i < orbitals.getE().size() ; i++ ) {
			String line ="     <eigenvalues id=\"e\""+ " value=\""+orbitals.getE().get(i)+"\">"+"</eigenvalues>";
			pw.println(line);
			line += "</e>";
		}
		pw.println("   </eigenvaluesArray>");
		pw.println("   <occupational number Array>");
		for ( i = 0 ; i < orbitals.getE().size() ; i++ ) {
			String line ="     <occupational number id=\"f\""+ " value=\""+orbitals.getF().get(i)+"\">"+"</occupational number>";
			pw.println(line);
			line += "</e>";
		}
		pw.println("   </occupational number Array>");
		pw.println("</orbitals>");
	}
	/**
	 * closing the PrintWriter 
	 */
	public void close() {
		if (pw != null) {
			pw.close();
			pw = null;
		}
	}
	/**
	 * @return the pw
	 */
	public PrintWriter getPw() {
		return pw;
	}
	/**
	 * @param pw the pw to set
	 */
	public void setPw(PrintWriter pw) {
		this.pw = pw;
	}
	/**
	 * @return the corbitals
	 */
	public COrbitals getCorbitals() {
		return corbitals;
	}
	/**
	 * @param corbitals the corbitals to set
	 */
	public void setCorbitals(COrbitals corbitals) {
		this.corbitals = corbitals;
	}

	
}
