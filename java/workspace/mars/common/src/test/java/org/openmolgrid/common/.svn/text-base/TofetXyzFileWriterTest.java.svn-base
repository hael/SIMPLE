package org.openmolgrid.common;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import junit.framework.Assert;

import org.junit.Test;
import org.openmolgrid.format.tofet.TofetXyzFileWriter;
import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.CSystem;

public class TofetXyzFileWriterTest {

	/**
	 * This test writes a Tofet XYZ file to the tmp directory
	 * 
	 * @throws IOException
	 */
	@Test
	public void testFileWriter() throws IOException {

		CStructure structure = new CStructure();
		CStructure structure2 = new CStructure();

		CAtom atom1 = new CAtom("H", 0.534, 1.1, 1.0, 1);
		atom1.setPartialCharge(0.9);
		structure.addAtom(atom1);

		CAtom atom2 = new CAtom("H", -1.0, -1.0, -1.0, 2);
		atom2.setPartialCharge(0.567);
		structure.addAtom(atom2);

		CAtom atom3 = new CAtom("H", -0.5, 1.1, 1.0, 1);
		atom1.setPartialCharge(0.9);
		structure2.addAtom(atom3);

		CAtom atom4 = new CAtom("H", 1.0, -1.0, 1.0, 2);
		atom2.setPartialCharge(0.567);
		structure2.addAtom(atom4);

		List<CStructure> structureList = new ArrayList<CStructure>();
		structureList.add(structure);
		structureList.add(structure2);

		CSystem system = new CSystem(structureList);

		File tempFile = File.createTempFile("tofet", "xyz");

		TofetXyzFileWriter.writeFile(system, tempFile.getAbsolutePath());

		Assert.assertEquals(true, tempFile.exists());
		Assert.assertEquals(true, tempFile.getTotalSpace() > 0);
		FileReader fr = null;
		BufferedReader reader = null;
		try {
			fr = new FileReader(tempFile);
			reader = new BufferedReader(fr);

			boolean hasCharge = false;
			while (true) {
				String input = reader.readLine();
				if (input == null) {
					break;
				}

				if (input.contains("0.05")) {
					hasCharge = true;
				}
				Assert.assertTrue(hasCharge);
			}
		} finally {
			reader.close();
			fr.close();
		}
	}
}
