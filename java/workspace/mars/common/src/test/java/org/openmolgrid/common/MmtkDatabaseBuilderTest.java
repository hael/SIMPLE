package org.openmolgrid.common;

import java.io.File;

import org.junit.Assert;
import org.junit.Test;
import org.openmolgrid.format.cml.CMLReader;
import org.openmolgrid.format.mmtk.MmtkDatabaseBuilder;
import org.openmolgrid.model.CStructure;

/**
 * Test class for the {@link MmtkDatabaseBuilder}
 * 
 * @author Stefan Bozic
 */
public class MmtkDatabaseBuilderTest {

	/**
	 * The methods creates a MmtkDatabase on the bas of a CML file.
	 */
	@Test
	public void createMoleculeDefinition() {
		try {

			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "cml" + File.separator + "Estron.cml");
//			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
//					+ "cml" + File.separator + "gs-4z_2.xml");

			
			CStructure structure = new CStructure();
			CMLReader cmlreader = new CMLReader(testFile.getAbsolutePath(), null, structure);
			cmlreader.provideStructures();
			String moleculeDefinition = MmtkDatabaseBuilder.createMoleculeDefinition("estron", structure);

			Assert.assertTrue(moleculeDefinition.contains("name = 'estron'"));
			System.out.println(moleculeDefinition);
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		}
	}
}