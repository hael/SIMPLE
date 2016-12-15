package org.openmolgrid.util;

import java.util.List;

import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CConstants;
import org.openmolgrid.model.CStructure;

/**
 * This class calculates the distance between two {@link CStructure} based on the center of mass.
 * 
 * @author Stefan Bozic
 */
public class DistanceCalculator {

	/**
	 * Calculates the distance between a given element of two structures. The element must have an occurrence of 1 in
	 * each structure.
	 * 
	 * @param struct1 structure 1
	 * @param struct2 structure 2
	 * @param element the element for the distance calculation
	 * 
	 * @return the distance between a given element of two structures
	 */
	public static double calculateDistanceBetweenUniqueElement(CStructure struct1, CStructure struct2, CConstants element) {
		double distance = 0;

		List<CAtom> atoms1 = struct1.getAtomsOfElement(element);
		List<CAtom> atoms2 = struct2.getAtomsOfElement(element);

		if (atoms1.size() != 1 || atoms2.size() != 1) {
			throw new RuntimeException(
					"The structures must contain only one atom of the given element! The first structure has "
							+ atoms1.size() + " and the second structure has " + atoms2.size() + " atoms of " + element);
		}

		distance = calculateDistance(atoms1.get(0), atoms2.get(0));

		return distance;
	}

	/**
	 * Calculates the distance between a given element of two structures. The element must have an occurrence of 1 in
	 * each structure.
	 * 
	 * @param struct1 structure 1
	 * @param struct2 structure 2 O *
	 * @return the distance between a given element of two structures
	 */
	public static double calculateDistanceBetweenCenterOfMass(CStructure struct1, CStructure struct2) {
		double distance = 0;

		final CenterOfMass com1 = CenterOfMassCalculator.calculatesCenterOfMass(struct1);
		final CenterOfMass com2 = CenterOfMassCalculator.calculatesCenterOfMass(struct2);

		distance = DistanceCalculator.calculateDistance(com1, com2);

		return distance;
	}

	/**
	 * Calculates the distance between the two atoms.
	 * 
	 * @param com1 the first center of mass for the calculation
	 * @param com2 the second center of mass for the calculation
	 * @return the distance between the center of mass of two structure
	 */
	public static double calculateDistance(CenterOfMass com1, CenterOfMass com2) {
		double distance = 0;

		distance = Math
				.sqrt((Math.pow((com1.getX() - com2.getX()), 2) + Math.pow((com1.getY() - com2.getY()), 2) + Math.pow(
						(com1.getZ() - com2.getZ()), 2)));

		return distance;
	}
	
	/**
	 * Calculates the distance between the two atoms.
	 * 
	 * @param atom1 the first atom for the calculation
	 * @param atom2 the first atom for the calculation
	 * @return the distance between the two atoms.
	 */
	public static double calculateDistance(CAtom atom1, CAtom atom2) {
		double distance = 0;

		distance = Math
				.sqrt((Math.pow((atom1.getX() - atom2.getX()), 2) + Math.pow((atom1.getY() - atom2.getY()), 2) + Math
						.pow((atom1.getZ() - atom2.getZ()), 2)));

		return distance;
	}
}
