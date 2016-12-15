package org.openmolgrid.format.cml;

import java.io.PrintWriter;
import java.io.StringWriter;

import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CBond;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.CSystem;
import org.openmolgrid.model.CBond.Stereo;

/**
 * Class that writes given {@link CSystem} and {@link CStructure} classes as CML
 * 
 * @author Stefan Bozic
 */
public class CMLWriter {

	/** XML Tag String */
	public final static String XML_TAG = "<?xml version=\"1.0\"?>";
	
	/** Cml Tag String with all namespaces and convention */
	public final static String CML_TAG_WITH_NAMESPACES = "<cml xmlns=\"http://www.xml-cml.org/schema\" " +
			"xmlns:cc=\"http://www.xml-cml.org/dictionary/compchem/\" " +
			"xmlns:nonsi=\"http://www.xml-cml.org/unit/nonSi/\" " +
			"xmlns:si=\"http://www.xml-cml.org/unit/si/\" " +
			"xmlns:convention=\"http://www.xml-cml.org/convention/\" " +
			"convention=\"convention:compchem\" " + 					
			">";

	/**
	 * A {@link PrintWriter} to print the CML markup to
	 */
	private PrintWriter pw;

	/**
	 * Constructor
	 * 
	 * @param pw The {@link PrintWriter} instance for this object
	 */
	public CMLWriter(PrintWriter pw) {
		this.pw = pw;
	}

	/**
	 * Returns the CML molecule markup representation for the given {@link CStructure} as {@link String}
	 * 
	 * @param str the {@link CStructure} to convert
	 * 
	 * @return CML markup
	 */
	public static String getCmlMolecule(CStructure str) {
		StringWriter sw = new StringWriter();
		new CMLWriter(new PrintWriter(sw)).write(str);

		return sw.toString();
	}

	/**
	 * Returns a full CML document for the given {@link CStructure} as {@link String}
	 * 
	 * @param str full CML document for the {@link CStructure}
	 * @return a full CML document
	 */
	public static String getCmlDocument(CStructure str) {
		StringWriter sw = new StringWriter();

		PrintWriter pw = new PrintWriter(sw);
		pw.println(XML_TAG);
		pw.println(CML_TAG_WITH_NAMESPACES);
		new CMLWriter(pw).write(str);
		pw.println("</cml>");

		return sw.toString();
	}
	/**
	 * Returns a full CML document for the given {@link CStructure}
	 * @param system
	 * @return sw.toString()
	 */
	public static String getCmlDocument(CSystem system) {
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);

		pw.println(XML_TAG);
		pw.println(CML_TAG_WITH_NAMESPACES);
		new CMLWriter(pw).write(system);
		pw.println("</cml>");

		return sw.toString();
	}

	/**
	 * Writes a CML representation for the given {@link CSystem}
	 * 
	 * @param system the system to convert to CML
	 */
	public void write(CSystem system) {
		for (CStructure structure : system.getStructures()) {
			write(structure);
		}
	}
	/**
	 * Writes a CML representation for the given {@link CSystem}
	 * @param str
	 */
	public void write(CStructure str) {
		pw.println("<molecule id=\"" + str.getId() + "\">");
		writeAtomArray(str);
		writeBondArray(str);
		pw.println("</molecule>");
	}

	/**
	 * Writes the BondArray for the given {@link CStructure} to the {@link PrintWriter}
	 * 
	 * @param str the {@link CStructure} with the bonds.
	 */
	private void writeAtomArray(CStructure str) {
		if (str.getAtomCount() > 0) {
			pw.println("<atomArray>");
			for (int i = 1; i <= str.getAtomCount(); i++) {
				CAtom a = str.getAtom(i);
				String line = "<atom id=\"a" + a.getId() + "\" elementType=\"" + a.getElement() + "\" x3=\"" + a.getX()
						+ "\" y3=\"" + a.getY() + "\" z3=\"" + a.getZ();
				if (a.getFormalCharge() != 0)
					line += "\" formalCharge=\"" + a.getFormalCharge();
				line += "\">";
				pw.println(line);

				if (a.getPartialCharge() != 0) {
					pw.println("<scalar dictRef=\"cc:charge\" dataType=\"xsd:double\" units=\"nonsi:elementaryCharge\">"
							+ a.getPartialCharge() + "</scalar>");
				}

				pw.println("</atom>");
			}
			pw.println("</atomArray>");
		} else {
			close();
			throw new IllegalArgumentException("A CStructure must have a least 1 atom.");
		}
	}

	/**
	 * Writes the BondArray for the given {@link CStructure} to the {@link PrintWriter}
	 * 
	 * @param str the {@link CStructure} with the bonds.
	 */
	private void writeBondArray(CStructure str) {
		if (str.getBondCount() > 0) {
			pw.println("<bondArray>");
			for (int i = 1; i <= str.getBondCount(); i++) {
				CBond b = str.getBond(i);
				String ln = "<bond id=\"b" + i + "\" atomRefs2=\"a" + b.getAtom1Id() + " a" + b.getAtom2Id()
						+ "\" order=\"" + getCmlBondOrder(b.getType());
				if (b.getStereo() != Stereo.NONE) {
					ln += "\">";
					pw.println(ln);
					String stereo = getCmlBondStereo(b.getStereo());
					pw.println(String.format("<bondStereo dictRef=\"cml:%s\">%s</bondStereo>", stereo, stereo));
					pw.println("</bond>");
				} else {
					ln += "\"/>";
					pw.println(ln);
				}
			}
			pw.println("</bondArray>");
		}
	}

	/**
	 * Closes the {@link PrintWriter}
	 */
	public void close() {
		if (pw != null) {
			pw.close();
			pw = null;
		}
	}

	/**
	 * Returns the CML bond order for a given value
	 * 
	 * @param order the value to convert
	 * 
	 * @return the cml bond order
	 */
	private String getCmlBondOrder(int order) {
		switch (order) {
		case 1:
			return "S";
		case 2:
			return "D";
		case 3:
			return "T";
		case 4:
			return "A";
		}
		return "";
	}

	/**
	 * Return the CML bond stereo
	 * 
	 * @param stereo the value to convert
	 * @return the CML bond stereo
	 */
	@SuppressWarnings("incomplete-switch")
	private String getCmlBondStereo(Stereo stereo) {
		switch (stereo) {
		case UP:
			return "W";
		case DOWN:
			return "H";
		}
		return "";
	}
}
