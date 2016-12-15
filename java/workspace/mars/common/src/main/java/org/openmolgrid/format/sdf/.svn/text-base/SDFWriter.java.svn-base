/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.format.sdf;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;


import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CBond;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.CBond.Stereo;
import org.openmolgrid.util.CFormatter;

public class SDFWriter {
	protected PrintWriter outWriter = null;

	public SDFWriter(String outFileName) {
		this(new File(outFileName));
	}

	public SDFWriter(File outFile) {
		try {
			outWriter = new PrintWriter(new FileWriter(outFile));
		} catch (IOException e) {
			throw new RuntimeException(e.toString());
		}
	}

	public SDFWriter(PrintWriter pw) {
		outWriter = pw;
	}

	public void addStructure(CStructure s) {
		int i;
		int[][] cp = new int[12][2];
		int ncharges = 0;
		outWriter.println();
		outWriter.println("   -OMG-");
		outWriter.println();
		outWriter.println(CFormatter.formatStringRight(String.valueOf(s.getAtomCount()), 3)
				+ CFormatter.formatStringRight(String.valueOf(s.getBondCount()), 3)
				+ "  0  0  0  0  0  0  0  0999 V2000");
		for (i = 1; i <= s.getAtomCount(); i++) {
			CAtom a = s.getAtom(i);
			if (a.getFormalCharge() != 0) {
				if (ncharges == cp.length)
					return;
				cp[ncharges][0] = i;
				cp[ncharges][1] = a.getFormalCharge();
				ncharges++;
			}
			outWriter.println(CFormatter.formatDouble(a.getX(), 10, "0.0000")
					+ CFormatter.formatDouble(a.getY(), 10, "0.0000")
					+ CFormatter.formatDouble(a.getZ(), 10, "0.0000")
					+ " "
					+ CFormatter.formatStringLeft(a.getElement(), 3)
					+ " 0"
					+ CFormatter.formatStringRight(String.valueOf(SDFWriter.translateToOldStyleCharge(a
							.getFormalCharge())), 3) + "  0  0  0");
		}
		for (i = 1; i <= s.getBondCount(); i++) {
			CBond b = s.getBond(i);
			outWriter.println(CFormatter.formatStringRight(String.valueOf(b.getAtom1Id()), 3)
					+ CFormatter.formatStringRight(String.valueOf(b.getAtom2Id()), 3)
					+ CFormatter.formatStringRight(String.valueOf(b.getType()), 3)
					+ CFormatter.formatStringRight(translateMdlStereo(b.getStereo()), 3) + "  0  0");
		}
		if (ncharges > 0) {
			outWriter.print("M  CHG" + CFormatter.formatStringRight(String.valueOf(ncharges), 3));
			for (i = 0; i < ncharges; i++) {
				outWriter.print(CFormatter.formatStringRight(String.valueOf(cp[i][0]), 4)
						+ CFormatter.formatStringRight(String.valueOf(cp[i][1]), 4));
			}
			outWriter.println();
		}
		outWriter.println("M  END");

		outWriter.println("$$$$");
	}

	public void close() {
		if (outWriter != null) {
			outWriter.close();
			outWriter = null;
		}
	}

	public static String translateMdlStereo(Stereo stereo) {
		switch (stereo) {
		case UP:
			return "1";
		case DOWN:
			return "6";
		}
		return "0";
	}

	private static int translateToOldStyleCharge(int c) {
		switch (c) {
		case 3:
			return 1;
		case 2:
			return 2;
		case 1:
			return 3;
		case -1:
			return 5;
		case -2:
			return 6;
		case -3:
			return 7;
		}
		return 0;
	}

}
