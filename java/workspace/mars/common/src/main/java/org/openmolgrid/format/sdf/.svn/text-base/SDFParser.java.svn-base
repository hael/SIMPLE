/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.format.sdf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;

import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CBond;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.CBond.Stereo;
import org.openmolgrid.qsar.DataDictionary;

/**
 * Parses MDL's SD and MOL files.
 * 
 * @author Andre Lomaka
 * @author Sulev Sild
 */
public class SDFParser
{
	private BufferedReader inReader;
	private String fileName;

	public SDFParser(String inFileName) {
		fileName = new File(inFileName).getName();
		try {
			inReader = new BufferedReader(new FileReader(inFileName));
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}
	
	public SDFParser(Reader reader) {
		if (reader instanceof BufferedReader) {
			inReader = (BufferedReader) reader;
		} else {
			inReader = new BufferedReader(reader);
		}
	}
	
	public CStructure getNextStructure() {
		SDFileEntry mdl = next();
		if (mdl == null) {
			return null;
		} else {
			return createCStructure(mdl);
		}
	}
	
	public void close() {
		if (inReader != null) {
			try {
				inReader.close();
			} catch (IOException e) {
				throw new RuntimeException("Error closing SDFParser");
			}
			inReader = null;
		}
	}

	public static class SDFileEntry {
		final public String[] molfile;
		final public DataDictionary data;

		public SDFileEntry(String[] molfile, DataDictionary data) {
			this.molfile = molfile;
			this.data = data;
		}
	};
	
	public SDFileEntry next() {
		String[] molfile = readMolfilePart();
		if (molfile == null) {
			return null;
		} else {
			return new SDFileEntry(molfile, readDataItems());
		}
	}
	
	private String[] readMolfilePart() {
		ArrayList<String> l = new ArrayList<String>();
		while (true) {
			String line = readLine();
			if (line == null) {
				return null;
			}
			l.add(line);
			if ("M  END".equals(line)) {
				break;
			}
		}
		return l.toArray(new String[l.size()]);
	}
	
	private DataDictionary readDataItems() {
		String name = null;
		StringBuilder value = new StringBuilder();
		DataDictionary data = new DataDictionary();
		if (fileName != null) {
			data.put("Filename (w/o extension)", fileName.split("\\.")[0]);
		}

		while (true) {
			String line = readLine();
			if (line == null || "$$$$".equals(line))
				break;

			// start of a data field
			if (line.startsWith("> ")) {
				name = line.split("[<>]")[2];
				continue;
			}
			// content of the data field
			if (name != null) {
				if (line.length() == 0) {
					data.put(name, value.toString());
					value.delete(0, value.length());
					name = null;
				}
				if (value.length() > 0) {
					value.append("\n");
				}
				value.append(line);
			}
		}

		return data;
	}
	
	private String readLine() {
		try {
			return inReader.readLine();
		} catch (IOException e) {
			throw new RuntimeException();
		}
	}
	
	private CStructure createCStructure(SDFileEntry mdl) {
		CStructure str = new CStructure();
	
		String[] lines = mdl.molfile;
		
		mdl.data.put("TITLE", lines[0].trim());
		
		int cur = 3;
		int na = Integer.parseInt(lines[cur].substring(0, 3).trim());
		int nb = Integer.parseInt(lines[cur].substring(3, 6).trim());
		for (int i = 0; i < na; i++) {
			++cur;
			double x = Double.parseDouble(lines[cur].substring(0, 10).trim());
			double y = Double.parseDouble(lines[cur].substring(10, 20).trim());
			double z = Double.parseDouble(lines[cur].substring(20, 30).trim());
			String element = lines[cur].substring(31, 34).trim();
			int fc = Integer.parseInt(lines[cur].substring(36, 39).trim());
			CAtom a = new CAtom(element, x, y, z, i + 1);
			if (fc != 0) {
				a.setFormalCharge(SDFParser.translateFromOldStyleCharge(fc));
			}
			str.addAtom(a);
		}
		
		for (int i = 0; i < nb; i++) {
			++cur;
			int bnda1 = Integer.parseInt(lines[cur].substring(0, 3).trim());
			int bnda2 = Integer.parseInt(lines[cur].substring(3, 6).trim());
			int mlp = Integer.parseInt(lines[cur].substring(6, 9).trim());
			int st = Integer.parseInt(lines[cur].substring(9, 12).trim());
			Stereo stereo = translateMdlStereo(st);
			CBond bond = new CBond(str.getAtom(bnda1), str.getAtom(bnda2), mlp, stereo);
			str.addBond(bond);
		}
		
		for (;cur<lines.length;++cur) {
			if ("M  END".equals(lines[cur])) {
				break;
			} else if (lines[cur].startsWith("M  CHG")) {
				// this property supersedes *all* charge values in the atom block
				for (int i = 0; i < na; i++) {
					str.getAtom(i+1).setFormalCharge(0);
				}
				int ii = Integer.parseInt(lines[cur].substring(6, 9).trim());
				for (int i = 0; i < ii; i++) {
					int an = Integer.parseInt(lines[cur].substring(10 + i * 8, 13 + i * 8).trim());
					int fc = Integer.parseInt(lines[cur].substring(14 + i * 8, 17 + i * 8).trim());
					str.getAtom(an).setFormalCharge(fc);
				}
			}
		}

		if (mdl.data.size() > 0) {
			str.addAttribute(mdl.data);
		}

		return str;
	}

	public static Stereo translateMdlStereo(int mdlStereo) {
		switch (mdlStereo) {
		case 1:
			return Stereo.UP;
		case 6:
			return Stereo.DOWN;
		}
		return Stereo.NONE;
	}

	private static int translateFromOldStyleCharge(int c) {
		switch (c) {
		case 1:
			return 3;
		case 2:
			return 2;
		case 3:
			return 1;
		case 5:
			return -1;
		case 6:
			return -2;
		case 7:
			return -3;
		}
		return 0;
	}

}
