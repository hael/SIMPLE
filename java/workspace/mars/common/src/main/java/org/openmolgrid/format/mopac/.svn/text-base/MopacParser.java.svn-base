package org.openmolgrid.format.mopac;

import java.io.IOException;
import java.io.LineNumberReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.apache.log4j.Logger;
import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CStructure;

/**
 * Parser for Mopac Output files. Creates a corresponding {@link CStructure} of the Mopac Output
 * 
 * @author Stefan Bozic
 */
public class MopacParser {

	/** Logger */
	private static Logger log = Logger.getLogger(MopacParser.class.getName());

	// Represents the final geometry calculated by Mopac
	private CStructure finalGeometry;

	// returns the number of filledLevels calculated by Mopac
	private int filledLevels = -1;

	// stores the list of eigenvalues calculated by Mopac
	private List<Double> eigenvalues;

	// stores the mopac version
	private int mopacVersion;

	/**
	 * Constructor
	 * 
	 * @param reader a Reader instance for the mopac output file
	 */
	public MopacParser(Reader reader) {
		try {
			finalGeometry = new CStructure();
			eigenvalues = new ArrayList<Double>();

			LineNumberReader lnr = new LineNumberReader(reader);
			parse(lnr);
		} catch (IOException e) {
			throw new RuntimeException(e.getMessage());
		}
	}

	/**
	 * The resulting {@link CStructure}
	 * 
	 * @return the final geometry
	 */
	public CStructure getFinalGeometry() {
		return finalGeometry;
	}

	/**
	 * Parses a mopac output file. The result is stored in the {@link CStructure}
	 * 
	 * @param lnr the reader for the file
	 * @return <code>true</code> if the file has been parses successfully.
	 * 
	 * @throws IOException an exception that might occurs
	 */
	private boolean parse(LineNumberReader lnr) throws IOException {
		double x;
		double y;
		double z;
		int index;

		String line = null;
		CAtom atom = null;
		List<String> vs = new ArrayList<String>();
		boolean coordReaded = false;
		boolean chargeReaded = false;

		// Deterine the mopac version

		while ((line = lnr.readLine()) != null) {

			// get the mopac version
			if (line.indexOf("MOPAC2009") != -1) {
				mopacVersion = 2009;
				log.debug("Mopac Version: " + 2009);
			} else if (line.indexOf("MOPAC2007") != -1) {
				mopacVersion = 2007;
				log.debug("Mopac Version: " + 2007);
			} else if (line.indexOf("MOPAC 7.0") != -1) {
				mopacVersion = 7;
				log.debug("Mopac Version: " + 7.0);
			}

			// extract the no of filled levels
			else if (line.indexOf("NO. OF FILLED LEVELS") != -1) {
				log.debug("");
				log.debug("Parse FILLED LEVELS");
				log.debug("----------------------------");
				String[] split = line.split("=");
				filledLevels = Integer.parseInt(split[1].trim());
				log.debug("Filled Levels=" + filledLevels);
			}

			// read in the eigenvalues
			else if (line.indexOf("EIGENVALUES") != -1) {
				log.debug("");
				log.debug("Parse Eigenvalues");
				log.debug("----------------------------");

				// Mopac 7.0 and 2007 contains a blank row after EIGENVALUES
				if ((mopacVersion == 2007 || mopacVersion == 7) && (line = lnr.readLine()) == null) {
					return false;
				}

				while (true) {
					try {
						if ((line = lnr.readLine()) == null || line.isEmpty()) {
							break;
						}

						log.debug("Eigenvalues: " + line);

						tokenize(vs, line, " ");

						for (String value : vs) {
							Double eigenValue = Double.parseDouble(value);
							eigenvalues.add(eigenValue);
						}
					} catch (NumberFormatException nfe) {
						break;
					}
				}
				log.debug("Eigenvalues: " + eigenvalues);
			}

			// The CARTESIAN COORDINATES can occurs two times in the mopac output file
			// The first appearance reflects the original input coordinates
			// The second appearance reflects the optimized coordinates by mopac

			// Im not sure if always both coordinates are rendered in the output file. So parse them both
			else if (line.indexOf("CARTESIAN COORDINATES") != -1) {
				log.debug("");
				log.debug("Parse Cartesian coordinates");
				log.debug("----------------------------");
				// blank
				if ((line = lnr.readLine()) == null) {
					return false;
				}

				// column headings
				if ((line = lnr.readLine()) == null) {
					return false;
				}

				// blank
				if ((line = lnr.readLine()) == null) {
					return false;
				}

				if ((line = lnr.readLine()) == null) {
					return false;
				}

				tokenize(vs, line);

				// Check the right size of the coordinates column!
				if (vs.size() != 5) {
					break;
				}

				int coordinatesColumnCount = 5;
				while (vs.size() == coordinatesColumnCount) {
					try {
						// Parse the current one
						x = Double.parseDouble((String) vs.get(2));
						y = Double.parseDouble((String) vs.get(3));
						z = Double.parseDouble((String) vs.get(4));
						index = Integer.valueOf(vs.get(0));

						// Created new atom instance if neither coordinates or charge has been read
						if (!coordReaded && !chargeReaded) {
							log.debug("Add new atom for " + vs.get(0));
							log.debug("with coordinates xyz " + x + " " + y + " " + z);
							atom = new CAtom(vs.get(1), x, y, z, index);
							finalGeometry.addAtom(atom);
						} else {
							atom = finalGeometry.getAtom(index);
							log.debug("Set coordinates for " + index + " : " + x + " " + y + " " + z);
							atom.setX(x);
							atom.setY(y);
							atom.setZ(z);
						}

						if ((line = lnr.readLine()) == null) {
							break;
						}

						tokenize(vs, line);
					} catch (NumberFormatException e) {
						// A NumberFormatException indicates that the first column is not a atom number
						// and so not a line of the CARTESIAN COORDINATES table
						break;
					}
				}
				coordReaded = true;
			} else if (line.indexOf("NET ATOMIC CHARGES") != -1) {
				log.debug("");
				log.debug("Parse net atomic charges");
				log.debug("------------------------");

				// blank
				if ((line = lnr.readLine()) == null) {
					return false;
				}

				// column headings
				if ((line = lnr.readLine()) == null) {
					return false;
				}

				if ((line = lnr.readLine()) == null) {
					return false;
				}

				tokenize(vs, line);

				// The tokensize is the column count of the charges table. This count can vary!
				while (true) {
					try {
						index = Integer.parseInt((String) vs.get(0));

						// Created new atom instance if neither coordinates or charge has been read
						if (!coordReaded && !chargeReaded) {
							log.debug("Add new atom for " + vs.get(0));
							atom = new CAtom(vs.get(1), 0, 0, 0, index);
							finalGeometry.addAtom(atom);
						} else {
							atom = finalGeometry.getAtom(index);
						}

						double partialCharge = Double.parseDouble((String) vs.get(2));
						atom.setPartialCharge(partialCharge);
						log.debug("set partial charge for " + atom.getId() + ": " + partialCharge);
						if ((line = lnr.readLine()) == null) {
							break;
						}

						tokenize(vs, line);
					} catch (NumberFormatException e) {
						// A NumberFormatException indicates that the first column is not a atom number
						// and so not a line of the charges table
						break;
					}
				}

				chargeReaded = true;
			}
		}

		return true;
	}

	/**
	 * @param vcr {@link java.util.Vector} of <tt>String</tt>
	 * @param buf Description of the Parameter
	 * @return Description of the Return Value
	 */
	public boolean tokenize(List<String> vcr, String buf) {
		return tokenize(vcr, buf, " \t\n");
	}

	/**
	 * @param vcr {@link java.util.Vector} of <tt>String</tt>
	 * @param buf Description of the Parameter
	 * @param delimstr Description of the Parameter
	 * @return Description of the Return Value
	 */
	public static boolean tokenize(List<String> vcr, String buf, String delimstr) {
		vcr.clear();
		buf = buf + "\n";

		StringTokenizer st = new StringTokenizer(buf, delimstr);

		while (st.hasMoreTokens()) {
			vcr.add(st.nextToken());
		}

		return true;
	}

	/**
	 * @param vcr {@link java.util.Vector} of <tt>String</tt>
	 * @param s Description of the Parameter
	 * @param delimstr Description of the Parameter
	 * @param limit Description of the Parameter
	 * @return Description of the Return Value
	 */
	public static boolean tokenize(List<String> vcr, String s, String delimstr, int limit) {
		vcr.clear();
		s = s + "\n";

		int endpos = 0;
		int matched = 0;

		StringTokenizer st = new StringTokenizer(s, delimstr);

		while (st.hasMoreTokens()) {
			String tmp = st.nextToken();
			vcr.add(tmp);

			matched++;

			if (matched == limit) {
				endpos = s.lastIndexOf(tmp);
				vcr.add(s.substring(endpos + tmp.length()));

				break;
			}
		}

		return true;
	}

	/**
	 * @return the mopacVersion
	 */
	public int getMopacVersion() {
		return mopacVersion;
	}

	/**
	 * @return the filledLevels
	 */
	public int getFilledLevels() {
		return filledLevels;
	}

	/**
	 * @return the eigenvalues
	 */
	public List<Double> getEigenvalues() {
		return eigenvalues;
	}

	/**
	 * Returns the HOMO value from the mopac calculation. The HOMO value is the value from the eigenvalues with the
	 * index equal to filledLevels
	 * 
	 * @return the HOMO value from the mopac calculation.
	 */
	public double getHomo() {
		if (filledLevels < 0 || getEigenvalues() == null || getEigenvalues().isEmpty()
				|| filledLevels > getEigenvalues().size()) {
			throw new RuntimeException("HOMO not available!");
		} else {
			return getEigenvalues().get(filledLevels);
		}
	}

	/**
	 * Returns the HOMO- value from the mopac calculation. HOMO- is the value from the eigenvalues with the index equal
	 * to filledLevels-1
	 * 
	 * @return the HOMO- value from the mopac calculation
	 */
	public double getHomoMinus() {
		if (filledLevels < 0 || getEigenvalues() == null || getEigenvalues().isEmpty()
				|| filledLevels > getEigenvalues().size()) {
			throw new RuntimeException("HOMO- not available!");
		} else {
			return getEigenvalues().get(filledLevels - 1);
		}
	}

	/**
	 * Returns the LUMO value from the mopac calculation. The LUMO value is the value from the eigenvalues with the
	 * index equal to filledLevels+1
	 * 
	 * @return the LUMO value from the mopac calculation
	 */
	public double getLumo() {
		if (filledLevels < 0 || getEigenvalues() == null || getEigenvalues().isEmpty()
				|| filledLevels > getEigenvalues().size()) {
			throw new RuntimeException("lumo value not available!");
		} else {
			return getEigenvalues().get(filledLevels + 1);
		}
	}

	/**
	 * Returns the LUMO+ value from the mopac calculation. LUMO+ is the value from the eigenvalues with the index equal
	 * to filledLevels+2
	 * 
	 * @return the LUMO+ value from the mopac calculation.
	 */
	public double getLumoPlus() {
		if (filledLevels < 0 || getEigenvalues() == null || getEigenvalues().isEmpty()
				|| filledLevels > getEigenvalues().size()) {
			throw new RuntimeException("homo value not available!");
		} else {
			return getEigenvalues().get(filledLevels + 2);
		}
	}

	/**
	 * Returns the difference between HOMO and HOMO-1
	 * 
	 * @return the difference between HOMO and HOMO-1
	 */
	public double getHomoSplit() {
		return getHomo() - getHomoMinus();
	}

	/**
	 * Returns the difference between LUMO+1 and LUMO.
	 * 
	 * @return the difference between LUMO+1 and LUMO
	 */
	public double getLumoSplit() {
		return getLumoPlus() - getLumo();
	}
}
