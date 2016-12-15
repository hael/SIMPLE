/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.format.smiles;

import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.ChemLibException;
import org.openmolgrid.model.Convertor;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.smiles.SmilesParser;

/**
 * Reads SMILES strings.
 * 
 * @author Sulev Sild
 */
public class SmilesReader {

	private StructureDiagramGenerator sdg = new StructureDiagramGenerator();

	public SmilesReader() {
	}

	public CStructure parseSmiles(String smiles) throws ChemLibException {
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

		try {
			IMolecule mol = sp.parseSmiles(smiles);
			sdg.setMolecule(mol);
			sdg.generateCoordinates();
			mol = sdg.getMolecule();
			CStructure str = Convertor.convert(new AtomContainer(mol), true);
			str.setId(smiles);
			str.setName(smiles);
			return str;
		} catch (InvalidSmilesException e) {
			throw new ChemLibException("Invalid SMILES string: " + e.getMessage());
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

}
