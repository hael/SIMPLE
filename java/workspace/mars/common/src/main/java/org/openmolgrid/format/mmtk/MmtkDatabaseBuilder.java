package org.openmolgrid.format.mmtk;

import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CBond;
import org.openmolgrid.model.CStructure;

/**
 * This class builds MMMTK databases. MmtkDatabaseBuilder not intended for building instances of it, so its declared as
 * abstract.
 * 
 * @author Stefan Bozic
 */
public abstract class MmtkDatabaseBuilder {

	/**
	 * Creates a new instance of {@link MmtkDatabaseBuilder}
	 * 
	 * @param moleculeName the name of the molecule
	 * @param structure the structure to build the database of
	 * @return a string representation of the mmtk database
	 */
	public static String createMoleculeDefinition(String moleculeName, CStructure structure) {
		StringBuffer buffer = new StringBuffer();

		buffer.append("name = '" + moleculeName + "'");
		buffer.append("\n");
		buffer.append(createAtomsForMoleculeDefinition(structure));
		buffer.append(createBondsForMoleculeDefinition(structure));
		buffer.append(createConfigurationForMoleculeDefinition(structure));
		buffer.append(createMolChargeForMoleculeDefinition(structure));

		return buffer.toString();
	}

	/**
	 * Creates a {@link String} that represents the atoms in a molecule definition of a MMTK database
	 * 
	 * @param structure the structure to parse
	 * 
	 * @return a {@link String} that represents the atoms in a molecule definition
	 */
	private static String createAtomsForMoleculeDefinition(CStructure structure) {
		StringBuffer buffer = new StringBuffer();
		try {
			for (CAtom curAtom : structure.getAtoms()) {
				buffer.append("a" + curAtom.getId() + " = Atom('" + curAtom.getElement() + "')");
				buffer.append("\n");
			}

		} catch (Exception e) {
			e.printStackTrace();
		}

		buffer.append("\n");
		buffer.append("\n");
		return buffer.toString();
	}

	/**
	 * Creates a {@link String} that represents the bonds in a molecule definition of a MMTK database
	 * 
	 * @param structure the structure to parse
	 * 
	 * @return a {@link String} that represents the bonds in a molecule definition
	 */
	private static String createBondsForMoleculeDefinition(CStructure structure) {
		StringBuffer buffer = new StringBuffer();
		try {
			buffer.append("bonds = [");

			int bondsAppended = 0;
			for (CBond curBond : structure.getBonds()) {
				buffer.append("Bond(");
				buffer.append("a" + curBond.getAtom(1).getId());
				buffer.append(", ");
				buffer.append("a" + curBond.getAtom(2).getId());
				buffer.append(")");
				bondsAppended++;

				if (bondsAppended < structure.getBondCount()) {
					buffer.append(",");
				}
			}
			buffer.append("]");

		} catch (Exception e) {
			e.printStackTrace();
		}
		buffer.append("\n");
		buffer.append("\n");
		return buffer.toString();
	}

	/**
	 * Creates a {@link String} that represents the configurations in a molecule definition of a MMTK database
	 * 
	 * @param structure the structure to parse
	 * 
	 * @return a {@link String} that represents the configurations in a molecule definition
	 */
	private static String createConfigurationForMoleculeDefinition(CStructure structure) {
		StringBuffer buffer = new StringBuffer();

		try {

			int atomsAppended = 0;

			buffer.append("configurations = {");
			buffer.append("\n");
			buffer.append("'default':Cartesian({");
			for (CAtom atom : structure.getAtoms()) {

				// Recalculate the coordinates from Angstrom to nm
				double xCoord = atom.getX() * 0.1;
				double yCoord = atom.getY() * 0.1;
				double zCoord = atom.getZ() * 0.1;

				buffer.append("a" + atom.getId() + ": (" + xCoord + ", " + yCoord + ", " + zCoord + ")");

				atomsAppended++;
				if (atomsAppended < structure.getAtomCount()) {
					buffer.append(",");
					buffer.append("\n");
				}

			}
			buffer.append("})");
			buffer.append("\n");
			buffer.append("}");

		} catch (Exception e) {
			e.printStackTrace();
		}
		buffer.append("\n");
		buffer.append("\n");
		return buffer.toString();
	}

	/**
	 * Creates a {@link String} that represents the mol_charges in a molecule definition of a MMTK database
	 * 
	 * @param structure the structure to parse
	 * 
	 * @return a {@link String} that represents the mol_charges in a molecule definition
	 */
	private static String createMolChargeForMoleculeDefinition(CStructure structure) {
		StringBuffer buffer = new StringBuffer();
		try {
			int atomsAppended = 0;

			buffer.append("mol_charge = {");
			for (CAtom atom : structure.getAtoms()) {
				buffer.append("a" + atom.getId() + ": ");

				if (atom.getPartialCharge() != 0) {
					buffer.append(atom.getPartialCharge());
				}
				// if no charge is defined set the value to 0
				else {
					buffer.append("0.000000");
				}

				atomsAppended++;
				if (atomsAppended < structure.getAtomCount()) {
					buffer.append(",");
					buffer.append("\n");
				}
			}
			buffer.append("}");

		} catch (Exception e) {
			e.printStackTrace();
		}
		buffer.append("\n");
		buffer.append("\n");
		return buffer.toString();
	}
}
