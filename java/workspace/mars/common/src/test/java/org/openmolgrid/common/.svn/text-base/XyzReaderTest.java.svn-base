package org.openmolgrid.common;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;
import org.openmolgrid.format.xyz.XyzReader;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.util.AtomDatabase;

public class XyzReaderTest {
	
	@Test
	public void testReadXyzFile() {
		try {
			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "xyz" + File.separator + "Estron.xyz");
			Assert.assertNotNull(testFile);
			CStructure structure = XyzReader.readXyzFile(testFile);
			Assert.assertNotNull(structure);
			Assert.assertTrue(structure.getAtomCount() == 168);
			Assert.assertEquals(7.66923, structure.getAtom(1).getZ(), 0.00000001);
			
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		}
	}
	
	
	@Test
	public void testReadXyz() {
		try {
			String testString = "2 \n \nC 0.000 1.000 3.000\nH 12.0 33.33 5544.3\n";
			CStructure structure = XyzReader.readXyz(testString);
			Assert.assertNotNull(structure);
			Assert.assertTrue(structure.getAtomCount() == 2);
			Assert.assertEquals(5544.3, structure.getAtom(2).getZ(), 0.00000001);
			
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		}
	}
	
	@Test
	public void testReadXyzToAtoms() {
		try {
			String testString = "2 \n \nC 0.000 1.000 3.000\nH 12.0 33.33 5544.3\n";
			ArrayList<AtomDatabase> aD = XyzReader.readXyzToAtoms(testString);
			Assert.assertNotNull(aD);
			Assert.assertTrue(aD.size() == 2);
			Assert.assertEquals(AtomDatabase.H.name(), aD.get(1).name());
			
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		}
	}
	
	@Test
	public void failtestReadXyzFile() {
		try {
			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "xyz" + File.separator + "benzene3.xyz");
			Assert.assertNotNull(testFile);
			CStructure structure = XyzReader.readXyzFile(testFile);
			Assert.fail("You should not get this far. An Exception had to be thrown.");
		} catch (Exception e) {
			Assert.assertTrue(e.getLocalizedMessage().contains("Could not read the file"));
		}
	}
	
	
}
