package org.openmolgrid.common;

import junit.framework.Assert;

import org.junit.Test;
import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CStructure;

/**
 * Test class for the {@link CStructure}
 * 
 * @author Stefan Bozic
 */
public class CStructureTest {

	/**
	 * The methods creates a CML {@link String} based on a {@link CStructure}
	 */
	@Test
	public void createCml() {
		CStructure structure = new CStructure();
		structure.setId("m1");
		CAtom atom1 = new CAtom("H", 1.0, 1.0, 1.0, 1);
		atom1.setPartialCharge(0.9);
		structure.addAtom(atom1);
		
		CAtom atom2 = new CAtom("H", -1.0, -1.0, -1.0, 2);
		atom2.setPartialCharge(0.567);
		structure.addAtom(atom2);		
		
		String cml = structure.writeCML2String();
		System.out.println(cml);
		Assert.assertNotNull(cml);
		boolean hasCharge = cml.contains("<scalar dictRef=\"cc:charge\" dataType=\"xsd:double\" units=\"nonsi:elementaryCharge\">0.567</scalar>");
		Assert.assertTrue(hasCharge);
	}
	
}
