package org.openmolgrid.util;

import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CStructure;

/**
 * This class contains two instances of {@link CStructure} that forms a pair. Two structures forms a pair if there
 * distance is below a given barrier.
 * 
 * @author Stefan Bozic
 */
public class StructurePair {

	/**
	 * The first structure of the pair
	 */
	private CStructure structure1;

	/**
	 * The second structure of the pair
	 */
	private CStructure structure2;

	/**
	 * A Structure that contains all atoms from structure1 and structure2.
	 */
	private CStructure combinedStructure = null;

	/**
	 * Creates a new {@link StructurePair}.
	 * 
	 * @param structure1 a CStructure
	 * @param structure2 a CStructure
	 */
	public StructurePair(CStructure structure1, CStructure structure2) {
		this.structure1 = structure1;
		this.structure2 = structure2;
	}

	/**
	 * Returns the first structure of the pair
	 * 
	 * @return the first structure of the pair
	 */
	public CStructure getStructure1() {
		return structure1;
	}

	/**
	 * Returns the second structure of the pair
	 * 
	 * @return the second structure of the pair
	 */
	public CStructure getStructure2() {
		return structure2;
	}

	@Override
	public String toString() {
		return this.structure1.getId() + " - " + this.structure2.getId();
	}

	/**
	 * Returns a {@link CStructure} that contains all atoms from structure1 and structure2. Use lazy initialization for
	 * creating the object.
	 * 
	 * @return a {@link CStructure} that contains all atoms from structure1 and structure2.
	 */
	public CStructure getSingleStructureRepresentation() {
		CStructure combinedStructure = new CStructure();
		combinedStructure.setId(structure1.getId() + "_" + structure2.getId());
		
		int atomId = 0;
		
		for (CAtom atom : structure1.getAtoms()) {		
			atomId++;
			combinedStructure.addAtom(new CAtom(atom.getElement(), atom.getX(),atom.getY(),atom.getZ(), atomId));
		}

		for (CAtom atom : structure2.getAtoms()) {
			atomId++;
			combinedStructure.addAtom(new CAtom(atom.getElement(), atom.getX(),atom.getY(),atom.getZ(), atomId));
		}
		
		return combinedStructure;
	}

}
