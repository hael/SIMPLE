/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.format.gamess;

import java.io.PrintWriter;
import java.io.Serializable;

import org.openmolgrid.qsar.ApplicationOutput;
import org.openmolgrid.qsar.CAttribute;

public class GamessOutput extends ApplicationOutput implements Serializable, CAttribute {

	public GamessOutput(String filePath, String content) {
		super(filePath, content);
	}

	public GamessOutput(ApplicationOutput m) {
		super(m);
	}

	public CAttribute copy() {
		return new GamessOutput(this);
	}

	public void writeCDR(PrintWriter pw) {
		pw.print("    <gamessOutput>");
		super.writeCDR(pw);
		pw.println("    </gamessOutput>");
	}

}
