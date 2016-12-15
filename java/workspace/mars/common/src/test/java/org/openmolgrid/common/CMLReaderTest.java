package org.openmolgrid.common;

import java.io.File;

import org.junit.Assert;
import org.junit.Test;
import org.openmolgrid.format.cml.CMLReader;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.CSystem;

/**
 * Test class for the {@link CMLReader}.
 *  
 * @author Stefan Bozic
 */
public class CMLReaderTest {

	/**
	 * reads a simple CML file and provides a {@link CStructure}
	 */
	@Test
	public void readSimpleCmlToStructure() {
		try {
			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "cml" + File.separator + "Estron.cml");
			
			CStructure structure = new CStructure();
			CMLReader cmlreader = new CMLReader(testFile.getAbsolutePath(), null, structure);
			cmlreader.provideStructures();

			Assert.assertEquals(-0.68264902, structure.getAtom(1).getPartialCharge(), 0.00000001);
			
			Assert.assertNotNull(structure);			
			Assert.assertTrue(structure.getAtomCount() == 168);
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		}
	}

	/**
	 * reads a simple CML file and provides a {@link CStructure}
	 */
	@Test
	public void readSingleMoleculeCmlToSystem() {
		try {
			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "cml" + File.separator + "Estron.cml");
			
			CSystem system = new CSystem();
			CMLReader cmlreader = new CMLReader(testFile.getAbsolutePath(), null, system);
			cmlreader.provideStructures();

			//System.out.println(system.writeCML2String());
			
			Assert.assertNotNull(system);
			Assert.assertTrue(system.getStructures().size() > 0);
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		}
	}
	
	/**
	 * reads a simple CML file and provides a {@link CStructure}
	 */
	@Test
	public void readMultiMoleculeCmlToSystem() {
		try {
			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "cml" + File.separator + "MultiMolecule.xml");
			
			CSystem system = new CSystem();
			CMLReader cmlreader = new CMLReader(testFile.getAbsolutePath(), null, system);
			cmlreader.provideStructures();

			//System.out.println(system.writeCML2String());
			
			Assert.assertNotNull(system);
			Assert.assertEquals(2, system.getStructures().size());
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		}
	}
	
	/**
	 * reads a simple CML file and provides a {@link CStructure}
	 */
	@Test
	public void readSi4hCmlToSystem() {
		try {
			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "cml" + File.separator + "Si4H.cml");
			
			CSystem system = new CSystem();
			CMLReader cmlreader = new CMLReader(testFile.getAbsolutePath(), null, system);
			cmlreader.provideStructures();
			
			Assert.assertNotNull(system);
			Assert.assertEquals(1, system.getStructures().size());
			Assert.assertEquals(5, system.getStructures().get(0).getAtomCount());
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		}
	}
	
}
