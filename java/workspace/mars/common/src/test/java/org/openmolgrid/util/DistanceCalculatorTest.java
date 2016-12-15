package org.openmolgrid.util;

import org.junit.Assert;
import org.junit.Test;
import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CConstants;
import org.openmolgrid.model.CStructure;

/**
 * Test class for the {@link DistanceCalculator}
 * 
 * @author Stefan Bozic
 */
public class DistanceCalculatorTest {

	/**
	 * Tests the calculation of the distance between two Atoms. All atoms in this test has a distance of one to each
	 * other.
	 */
	@Test
	public void calculateDistanceBetweenAtomsWithDistanceOneCenterOfMass() {
		try {

			// Test1
			CAtom atom1 = new CAtom(CConstants.B.name(), 0.0, 0.0, 0.0, 0);
			CStructure struct1 = new CStructure();
			struct1.addAtom(atom1);

			CAtom atom2 = new CAtom(CConstants.B.name(), 0.0, 0.0, 1.0, 0);
			CStructure struct2 = new CStructure();
			struct2.addAtom(atom2);

			double distance = DistanceCalculator.calculateDistanceBetweenCenterOfMass(struct1, struct2);

			Assert.assertEquals(1.0, distance, 0.0);
		} catch (Exception e) {
			Assert.fail();
		}
	}

	/**
	 * Tests the calculation of the distance between two Atoms. The distance should b o, because both atoms have the
	 * same coordinates.
	 */
	@Test
	public void calculateDistanceBetweenAtomsWithSameCoordinates() {
		try {
			CAtom atom1 = new CAtom(CConstants.B.name(), 0.0, 0.0, 0.0, 0);
			CAtom atom2 = new CAtom(CConstants.B.name(), 0.0, 0.0, 0.0, 0);

			double distance = DistanceCalculator.calculateDistance(atom1, atom2);

			Assert.assertEquals(0.0, distance, 0.0);
		} catch (Exception e) {
			Assert.fail();
		}
	}

	/**
	 * Tests the calculation of the distance between two Atoms. All atoms in this test has a distance of one to each
	 * other.
	 */
	@Test
	public void calculateDistanceBetweenAtomsWithDistanceOne() {
		try {

			// Test1
			CAtom atom1 = new CAtom(CConstants.B.name(), 0.0, 0.0, 0.0, 0);
			CAtom atom2 = new CAtom(CConstants.B.name(), 0.0, 0.0, 1.0, 0);

			double distance = DistanceCalculator.calculateDistance(atom1, atom2);

			Assert.assertEquals(1.0, distance, 0.0);

			// Test2
			atom1 = new CAtom(CConstants.B.name(), 0.0, 0.0, 0.0, 0);
			atom2 = new CAtom(CConstants.B.name(), 0.0, 1.0, 0.0, 0);

			distance = DistanceCalculator.calculateDistance(atom1, atom2);

			Assert.assertEquals(1.0, distance, 0.0);

			// Test3
			atom1 = new CAtom(CConstants.B.name(), 0.0, 0.0, 0.0, 0);
			atom2 = new CAtom(CConstants.B.name(), 1.0, 0.0, 0.0, 0);

			distance = DistanceCalculator.calculateDistance(atom1, atom2);

			Assert.assertEquals(1.0, distance, 0.0);

			// Test4
			atom1 = new CAtom(CConstants.B.name(), 0.0, 0.0, 1.0, 0);
			atom2 = new CAtom(CConstants.B.name(), 0.0, 0.0, 0.0, 0);

			distance = DistanceCalculator.calculateDistance(atom1, atom2);

			Assert.assertEquals(1.0, distance, 0.0);

			// Test5
			atom1 = new CAtom(CConstants.B.name(), 0.0, 1.0, 0.0, 0);
			atom2 = new CAtom(CConstants.B.name(), 0.0, 0.0, 0.0, 0);

			distance = DistanceCalculator.calculateDistance(atom1, atom2);

			Assert.assertEquals(1.0, distance, 0.0);

			// Test6
			atom1 = new CAtom(CConstants.B.name(), 1.0, 0.0, 0.0, 0);
			atom2 = new CAtom(CConstants.B.name(), 0.0, 0.0, 0.0, 0);

			distance = DistanceCalculator.calculateDistance(atom1, atom2);

			Assert.assertEquals(1.0, distance, 0.0);
		} catch (Exception e) {
			Assert.fail();
		}
	}

	/**
	 * Tests the calculation of the distance between two Atoms.
	 */
	@Test
	public void calculateDistanceBetweenAtoms() {
		try {

			// Test1
			CAtom atom1 = new CAtom(CConstants.B.name(), 1.0, 1.0, 0.0, 0);
			CAtom atom2 = new CAtom(CConstants.B.name(), 0.0, 0.0, 0.0, 0);

			double distance = DistanceCalculator.calculateDistance(atom1, atom2);

			Assert.assertEquals(1.4142135623730950488016887242097, distance, 0.000001);

			// Test2
			atom1 = new CAtom(CConstants.B.name(), 1.0, 1.0, 0.0, 0);
			atom2 = new CAtom(CConstants.B.name(), 0.0, 0.8, 0.5, 0);

			distance = DistanceCalculator.calculateDistance(atom1, atom2);

			Assert.assertEquals(1.1357816691600547221784675967985, distance, 0.000001);
		} catch (Exception e) {
			Assert.fail();
		}
	}

	/**
	 * Tests the calculation of the distance between two Structures. The structures
	 */
	@Test
	public void calculateDistanceBetweenStructures() {
		// Test 1 should fail because structure 2 has no B element
		try {
			CAtom atom1 = new CAtom(CConstants.B.name(), 1.0, 1.0, 0.0, 0);

			CStructure structure1 = new CStructure();
			structure1.addAtom(atom1);

			CAtom atom2 = new CAtom(CConstants.Ac.name(), 0.0, 0.0, 0.0, 0);

			CStructure structure2 = new CStructure();
			structure2.addAtom(atom2);

			DistanceCalculator.calculateDistanceBetweenUniqueElement(structure1, structure2, CConstants.B);

			Assert.fail();
		} catch (Exception e) {
			// This point must be reached because atom1 and atom2 are different
			Assert.assertTrue(true);
		}

		// Test2
		try {
			CAtom atom1 = new CAtom(CConstants.B.name(), 1.0, 1.0, 0.0, 0);
			CAtom atom2 = new CAtom(CConstants.Ac.name(), 0.0, 0.0, 0.0, 0);

			CStructure structure1 = new CStructure();
			structure1.addAtom(atom1);
			structure1.addAtom(atom2);

			CAtom atom3 = new CAtom(CConstants.Ac.name(), 1.0, 1.0, 0.0, 0);
			CAtom atom4 = new CAtom(CConstants.B.name(), 0.0, 0.0, 0.0, 0);

			CStructure structure2 = new CStructure();
			structure2.addAtom(atom3);
			structure2.addAtom(atom4);

			// calculate the distance between the B atoms
			double distance = DistanceCalculator.calculateDistanceBetweenUniqueElement(structure1, structure2,
					CConstants.B);
			Assert.assertEquals(1.4142135623730950488016887242097, distance, 0.000001);

			// calculate the distance between the Ac atoms
			distance = DistanceCalculator.calculateDistanceBetweenUniqueElement(structure1, structure2,
					CConstants.Ac);
			Assert.assertEquals(1.4142135623730950488016887242097, distance, 0.000001);

		} catch (Exception e) {
			e.printStackTrace();
		}

		// Test2
		try {
			CAtom atom1 = new CAtom(CConstants.Al.name(), 10.0, 8.0, 9.0, 0);
			CAtom atom2 = new CAtom(CConstants.Ac.name(), 0.0, 0.0, 0.0, 0);

			CStructure structure1 = new CStructure();
			structure1.addAtom(atom1);
			structure1.addAtom(atom2);

			CAtom atom3 = new CAtom(CConstants.Ac.name(), 1.0, 1.0, 0.0, 0);
			CAtom atom4 = new CAtom(CConstants.Al.name(), 4.0, 3.0, 2.0, 0);

			CStructure structure2 = new CStructure();
			structure2.addAtom(atom3);
			structure2.addAtom(atom4);

			// calculate the distance between the B atoms
			double distance = DistanceCalculator.calculateDistanceBetweenUniqueElement(structure1, structure2,
					CConstants.Al);
			Assert.assertEquals(10.488088481701515, distance, 0.0000000000000000001);
		} catch (Exception e) {
			Assert.fail();
		}
	}
}
