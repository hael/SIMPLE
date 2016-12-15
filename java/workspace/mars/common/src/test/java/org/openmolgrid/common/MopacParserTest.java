package org.openmolgrid.common;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;

import junit.framework.Assert;

import org.junit.Test;
import org.openmolgrid.format.mopac.MopacParser;
import org.openmolgrid.model.CStructure;

/**
 * Test class for the {@link MopacParser}.
 * 
 * @author Stefan Bozic
 */
public class MopacParserTest {

	/**
	 * This test parses a mopac result file.
	 */
	@Test
	public void testMopacParser() {
		Reader reader = null;
		try {
			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "mopac" + File.separator + "o2h.mopout");

			reader = new FileReader(testFile);
			MopacParser parser = new MopacParser(reader);
			CStructure structure = parser.getFinalGeometry();

			
			Assert.assertNotNull(structure);
			
			Assert.assertEquals(7, parser.getMopacVersion());
			
			Assert.assertEquals(9, structure.getAtomCount());

			Assert.assertEquals(10, parser.getFilledLevels());

			Assert.assertEquals(18, parser.getEigenvalues().size());

		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * This test parses a mopac result file.
	 */
	@Test
	public void test2MopacParser() {
		Reader reader = null;
		try {

			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "mopac" + File.separator + "Formaldehyde.out");

			reader = new FileReader(testFile);
			MopacParser parser = new MopacParser(reader);
			CStructure structure = parser.getFinalGeometry();

			Assert.assertNotNull(structure);
			
			Assert.assertEquals(2007, parser.getMopacVersion());
			
			Assert.assertEquals(4, structure.getAtomCount());

			Assert.assertEquals(6, parser.getFilledLevels());

			Assert.assertEquals(10, parser.getEigenvalues().size());

		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * This test parses a mopac result file.
	 */
	@Test
	public void test3MopacParser() {
		Reader reader = null;
		try {
			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "mopac" + File.separator + "pdt_15_R.out");

			reader = new FileReader(testFile);
			MopacParser parser = new MopacParser(reader);
			CStructure structure = parser.getFinalGeometry();

			Assert.assertNotNull(structure);
			
			Assert.assertEquals(2007, parser.getMopacVersion());
			
			Assert.assertEquals(55, structure.getAtomCount());

			Assert.assertEquals(0, parser.getFilledLevels());

			Assert.assertEquals(433, parser.getEigenvalues().size());

		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * This test parses a mopac result file.
	 */
	@Test
	public void test4MopacParser() {
		Reader reader = null;
		try {
			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "mopac" + File.separator + "ammoniak.mop");

			reader = new FileReader(testFile);
			MopacParser parser = new MopacParser(reader);
			CStructure structure = parser.getFinalGeometry();

			Assert.assertNotNull(structure);
			
			Assert.assertEquals(2009, parser.getMopacVersion());
			
			Assert.assertEquals(4, structure.getAtomCount());

			Assert.assertEquals(4, parser.getFilledLevels());

			Assert.assertEquals(7, parser.getEigenvalues().size());

		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
}
