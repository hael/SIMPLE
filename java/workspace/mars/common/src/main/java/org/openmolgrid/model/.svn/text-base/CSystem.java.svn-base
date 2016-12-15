package org.openmolgrid.model;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import org.openmolgrid.format.cml.CMLWriter;

/**
 * This class represents a chemical system. A system is build of several chemical structures.
 * 
 * @author Stefan Bozic
 */
public class CSystem {	
	
	/** A {@link List} of structures that builds a chemical system. */
	List<CStructure> structures = null;

	/**
	 * Creates an empty system.
	 */
	public CSystem() {
		structures = new ArrayList<CStructure>();
	}

	/**
	 * Creates a system based on the given {@link CStructure}
	 * 
	 * @param structures the structures that buils the system
	 */
	public CSystem(List<CStructure> structures) {
		this.structures = structures;
	}

	/**
	 * @return the structures
	 */
	public List<CStructure> getStructures() {
		return structures;
	}

	/**
	 * @param structures the structures to set
	 */
	public void setStructures(List<CStructure> structures) {
		this.structures = structures;
	}

	/**
	 * Adds a {@link CStructure} to the system.
	 * 
	 * @param structure the structure to add
	 */
	public void addStructure(CStructure structure) {
		this.structures.add(structure);
	}

	/**
	 * Removes a {@link CStructure} from the system
	 * 
	 * @param structure the structure to remove
	 * 
	 * @return <code>true</code> if the structure has been removed successfully
	 */
	public boolean removeStructure(CStructure structure) {
		boolean retVal = this.structures.remove(structure);

		return retVal;
	}

	/**
	 * Returns a structure of the system that corresponds to the given id.
	 *  
	 * @param structureId The id of a structure to search for
	 * @return a structure that fits to the id
	 */
	public CStructure getStructurebyId(String structureId) {
		CStructure retVal = null;
		
		for (CStructure structure : structures) {
			if (structureId.equals(structure.getId())) {
				retVal = structure;
				break;
			}
		}
		
		return retVal;
	}
	
	
	/**
	 * Returns a {@link String} that contains a full CML document representation of the system.
	 * 
	 * @return a CML document representation of the system
	 */
	public String writeCML2String() {		
		return CMLWriter.getCmlDocument(this);
	}

}
