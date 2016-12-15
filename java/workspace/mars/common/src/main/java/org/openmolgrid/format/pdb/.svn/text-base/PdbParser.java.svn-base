package org.openmolgrid.format.pdb;

import java.io.IOException;
import java.io.LineNumberReader;
import java.io.Reader;

import org.apache.log4j.Logger;
import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CStructure;

/**
 * Parser for pdb files. Creates a {@link CStructure} from a pdb file.
 * 
 * @author Stefan Bozic
 */
public class PdbParser {

	/** Logger */
	private static Logger log = Logger.getLogger(PdbParser.class.getName());

	/** The final geometry after parsing */
	private CStructure finalGeometry;

	/**
	 * Constructor
	 * 
	 * @param reader a Reader instance for a pdb file
	 */
	public PdbParser(Reader reader) {
		try {
			finalGeometry = new CStructure();
			LineNumberReader lnr = new LineNumberReader(reader);
			parse(lnr);
		} catch (IOException e) {
			throw new RuntimeException(e.getMessage());
		}
	}

	/**
	 * Constructor
	 * 
	 * @param reader a Reader instance for the mopac output file
	 * @param moleculeSize the number of atoms for a molecule
	 */
	public PdbParser(Reader reader, int moleculeSize) {
		try {
			finalGeometry = new CStructure();
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
	 * Parses a pdb file. The result is stored in the {@link CStructure}
	 * 
	 * @param lnr the reader for the file
	 * @return <code>true</code> if the file has been parses successfully.
	 * 
	 * @throws IOException an exception that might occurs
	 */
	private boolean parse(LineNumberReader lnr) throws IOException {
		String line;
		CAtom atom;

		while ((line = lnr.readLine()) != null) {
			if (line.startsWith("ATOM") || line.startsWith("HETATM")) {
				PdbAtomRecord record = new PdbAtomRecord(line);
				int index = Integer.parseInt(record.getSerial());
				double x = Double.parseDouble(record.getX());
				double y = Double.parseDouble(record.getY());
				double z = Double.parseDouble(record.getZ());

				atom = new CAtom(record.getElement(), x, y, z, index);
				finalGeometry.addAtom(atom);
			}
		}

		return true;
	}
}