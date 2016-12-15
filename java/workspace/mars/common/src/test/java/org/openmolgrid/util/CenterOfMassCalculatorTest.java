package org.openmolgrid.util;

import org.junit.Assert;
import org.junit.Test;
import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CConstants;
import org.openmolgrid.model.CStructure;

public class CenterOfMassCalculatorTest {
	/**
	 * Calculates the Center of Mass of water.
	 */
	@Test
	public void calculateCenterOfMassOfWater() {
		try {
			double L = 0.91;

			CAtom atom1 = new CAtom(CConstants.H.name(), L, 0.0, 0.0, 0);
			CAtom atom2 = new CAtom(CConstants.H.name(), L * Math.cos(Math.toRadians(105)), L
					* Math.sin(Math.toRadians(105)), 0.0, 0);
			CAtom atom3 = new CAtom(CConstants.O.name(), 0.0, 0.0, 0.0, 0);

			CStructure structure = new CStructure();

			structure.addAtom(atom1);
			structure.addAtom(atom2);
			structure.addAtom(atom3);

			CenterOfMass com = CenterOfMassCalculator.calculatesCenterOfMass(structure);

			Assert.assertNotNull(com);
			Assert.assertEquals(0.0377, com.getX(), 0.0001);
			Assert.assertEquals(0.0492, com.getY(), 0.0001);
			Assert.assertEquals(0.0, com.getZ(), 0.001);

		} catch (Exception e) {
			Assert.fail();
		}
	}

	/**
	 * Calculates the Center of Mass of water.
	 */
	@Test
	public void calculateCenterOfSulfurOxid() {
		try {
			double L = 0.143;

			CAtom atom1 = new CAtom(CConstants.S.name(), 0.0, L / 2, 0.0, 0);
			CAtom atom2 = new CAtom(CConstants.O.name(), 0.0, -L * Math.sin(Math.toRadians(60)), 0.0, 0);
			CAtom atom3 = new CAtom(CConstants.O.name(), 0.0, L * Math.sin(Math.toRadians(60)), 0.0, 0);

			CStructure structure = new CStructure();

			structure.addAtom(atom1);
			structure.addAtom(atom2);
			structure.addAtom(atom3);

			CenterOfMass com = CenterOfMassCalculator.calculatesCenterOfMass(structure);

			Assert.assertNotNull(com);
			Assert.assertEquals(0.0, com.getX(), 0.0);
			Assert.assertEquals(0.03575, com.getY(), 0.0001);
			Assert.assertEquals(0.0, com.getZ(), 0.0);

		} catch (Exception e) {
			Assert.fail();
		}
	}
}
