package org.openmolgrid.util;

import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CStructure;

/**
 * Class for calculating the center of mass of a given {@link CStructure}
 * 
 * @author Stefan Bozic
 */
public class CenterOfMassCalculator {

	/**
	 * Calculates the center of mass for a given {@link CStructure}.
	 * 
	 * @param structure the structure for which to calculate the center of mass
	 * 
	 * @return the center of mass for this structure
	 */
	public static CenterOfMass calculatesCenterOfMass(CStructure structure) {
		CenterOfMass retVal = null;

		double sumFactorMassAndX = 0;
		double sumFactorMassAndY = 0;
		double sumFactorMassAndZ = 0;
		double sumMass = 0;

		for (CAtom atom : structure.getAtoms()) {
			AtomDatabase atomdata = AtomDatabase.valueOf(atom.getElement());

			sumFactorMassAndX += (atom.getX() * atomdata.getAtomicWeight());
			sumFactorMassAndY += (atom.getY() * atomdata.getAtomicWeight());
			sumFactorMassAndZ += (atom.getZ() * atomdata.getAtomicWeight());
			sumMass += atomdata.getAtomicWeight();
		}

		double centerOfMassX = sumFactorMassAndX / sumMass;
		double centerOfMassY = sumFactorMassAndY / sumMass;
		double centerOfMassZ = sumFactorMassAndZ / sumMass;

		retVal = new CenterOfMass(structure, centerOfMassX, centerOfMassY, centerOfMassZ);

		return retVal;
	}
	
}
