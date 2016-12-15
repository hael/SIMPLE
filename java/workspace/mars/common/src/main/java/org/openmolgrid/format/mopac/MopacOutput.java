/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.format.mopac;

import java.io.PrintWriter;

import org.openmolgrid.qsar.ApplicationOutput;
import org.openmolgrid.qsar.CAttribute;

public class MopacOutput extends ApplicationOutput implements CAttribute {

	public MopacOutput(ApplicationOutput m) {
		super(m);
	}

	public MopacOutput(String fp, String content) {
		super(fp, content);
	}

	public CAttribute copy() {
		return new MopacOutput(this);
	}

	public void writeCDR(PrintWriter pw) {
		pw.print("    <mopacOutput>");
		super.writeCDR(pw);
		pw.println("    </mopacOutput>");
	}

}
