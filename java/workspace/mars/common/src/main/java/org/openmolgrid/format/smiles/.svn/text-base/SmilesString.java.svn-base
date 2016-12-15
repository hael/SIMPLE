/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.format.smiles;

import java.io.PrintWriter;
import java.io.Serializable;

import org.openmolgrid.qsar.CAttribute;

public class SmilesString implements CAttribute, Serializable {

	private String smiles;

	public SmilesString(String smiles) {
		this.smiles = smiles;
	}

	public SmilesString(SmilesString source) {
		smiles = source.smiles;
	}

	public CAttribute copy() {
		return new SmilesString(this);
	}

	public void writeCDR(PrintWriter pw) {
		pw.println("    <smiles>" + smiles + "</smiles>");
	}

	@Override
	public String toString() {
		return smiles;
	}

}
