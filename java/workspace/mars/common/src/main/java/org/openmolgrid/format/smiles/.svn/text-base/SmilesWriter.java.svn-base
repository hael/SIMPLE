/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.format.smiles;

import java.io.IOException;
import java.io.StringWriter;

import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.Convertor;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.io.SMILESWriter;

public class SmilesWriter {

	public String getSmiles(CStructure s) {
		StringWriter smiles = new StringWriter();
		SMILESWriter writer = new SMILESWriter(smiles);
		Molecule mol = Convertor.convertMolecule(s);
		try {
			writer.write(mol);
			writer.close();
			return smiles.toString();
		} catch (CDKException e) {
			return "";
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

}
