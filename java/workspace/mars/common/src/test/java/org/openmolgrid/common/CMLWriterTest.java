/**
 * 
 */
package org.openmolgrid.common;

import junit.framework.Assert;

import org.junit.Test;
import org.openmolgrid.format.cml.CMLWriter;
import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.CSystem;

/**
 * Test class for the {@link CMLWriter}
 * 
 * @author Frederic Bonnet
 */
public class CMLWriterTest {
	/**
	 * write a simple CML file
	 */
	@Test
	public void writeSimpleCmlFromStructure() {
		
		try {
			CSystem system = new CSystem();
			CStructure structure = new CStructure();
			structure.setId("m1");
			CAtom atom1 = new CAtom("H", 1.0, 1.0, 1.0, 1);
			atom1.setPartialCharge(0.9);
			structure.addAtom(atom1);
			
			CAtom atom2 = new CAtom("H", -1.0, -1.0, -1.0, 2);
			atom2.setPartialCharge(0.567);
			structure.addAtom(atom2);		
			system.addStructure(structure);
			
			String cml = system.writeCML2String();
			System.out.println(cml);
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		}
	}
}
