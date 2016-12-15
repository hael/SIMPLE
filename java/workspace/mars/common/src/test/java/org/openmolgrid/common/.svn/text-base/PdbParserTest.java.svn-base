package org.openmolgrid.common;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;

import junit.framework.Assert;

import org.junit.Test;
import org.openmolgrid.format.mopac.MopacParser;
import org.openmolgrid.format.pdb.PdbAtomRecord;
import org.openmolgrid.format.pdb.PdbParser;
import org.openmolgrid.model.CStructure;

/**
 * Test class for the {@link MopacParser}.
 * 
 * @author Stefan Bozic
 */
public class PdbParserTest {

	/**
	 * This test parses a mopac result file and creates an instance of {@link CStructure}
	 */
	@Test
	public void testPdbParser() {
		Reader reader = null;

		try {
			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "deposit" + File.separator + "last.pdb");

			reader = new FileReader(testFile);
			PdbParser parser = new PdbParser(reader);

			CStructure structure = parser.getFinalGeometry();
			Assert.assertNotNull(structure);
			Assert.assertEquals(208, structure.getAtomCount());			
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
	 * This tests parses a pdb atom record.
	 */
	@Test
	public void testPdbRecord1() {
		// PdbAtomRecord record = new
		// PdbAtomRecord("ATOM     23  C23 MOL     1       2.028   4.571  -2.499  1.00  0.00              ");

		PdbAtomRecord record = new PdbAtomRecord(
				"HETATM 1357 MG    MG   168       4.669  34.118  19.123  1.00  3.16          MG2+");

		System.out.println("getRecname " + record.getRecname());
		System.out.println("getSerial " + record.getSerial());
		System.out.println("getAtom " + record.getAtom());
		System.out.println("getAltLocation " + record.getAltLocation());
		System.out.println("getResName " + record.getResName());
		System.out.println("getChainId " + record.getChainId());
		System.out.println("getSeqNo " + record.getSeqNo());
		System.out.println("getY " + record.getY());
		System.out.println("getZ " + record.getZ());
		System.out.println("getOccupancy " + record.getOccupancy());
		System.out.println("getTempFactor " + record.getTempFactor());
		System.out.println("getRecId " + record.getRecId());
		System.out.println("getSegId " + record.getSegId());
		System.out.println("getElement " + record.getElement());
		System.out.println("getCharge " + record.getCharge());

		Assert.assertEquals("HETATM", record.getRecname());
		Assert.assertEquals("1357", record.getSerial());
		Assert.assertEquals("MG", record.getAtom());
		Assert.assertEquals("", record.getAltLocation());
		Assert.assertEquals("MG", record.getResName());
		Assert.assertEquals("", record.getChainId());
		Assert.assertEquals("168", record.getSeqNo());
		Assert.assertEquals("4.669", record.getX());
		Assert.assertEquals("34.118", record.getY());
		Assert.assertEquals("19.123", record.getZ());
		Assert.assertEquals("1.00", record.getOccupancy());
		Assert.assertEquals("3.16", record.getTempFactor());
		Assert.assertEquals("MG", record.getElement());
		Assert.assertEquals("2+", record.getCharge());
	}

	/**
	 * This tests parses a pdb atom record.
	 */
	@Test
	public void testPdbRecord2() {
		PdbAtomRecord record = new PdbAtomRecord(
				"ATOM     23  C23 MOL     1       2.028   4.571  -2.499  1.00  0.00              ");

		System.out.println("getRecname " + record.getRecname());
		System.out.println("getSerial " + record.getSerial());
		System.out.println("getAtom " + record.getAtom());
		System.out.println("getAltLocation " + record.getAltLocation());
		System.out.println("getResName " + record.getResName());
		System.out.println("getChainId " + record.getChainId());
		System.out.println("getSeqNo " + record.getSeqNo());
		System.out.println("getY " + record.getY());
		System.out.println("getZ " + record.getZ());
		System.out.println("getOccupancy " + record.getOccupancy());
		System.out.println("getTempFactor " + record.getTempFactor());
		System.out.println("getRecId " + record.getRecId());
		System.out.println("getSegId " + record.getSegId());
		System.out.println("getElement " + record.getElement());
		System.out.println("getCharge " + record.getCharge());

		Assert.assertEquals("getRecname", "ATOM", record.getRecname());
		Assert.assertEquals("getSerial", "23", record.getSerial());
		Assert.assertEquals("getAtom", "C23", record.getAtom());
		Assert.assertEquals("getAltLocation", "", record.getAltLocation());
		Assert.assertEquals("getResName", "MOL", record.getResName());
		Assert.assertEquals("getChainId", "", record.getChainId());
		Assert.assertEquals("getSeqNo", "1", record.getSeqNo());
		Assert.assertEquals("getX", "2.028", record.getX());
		Assert.assertEquals("getY", "4.571", record.getY());
		Assert.assertEquals("getZ", "-2.499", record.getZ());
		Assert.assertEquals("getOccupancy", "1.00", record.getOccupancy());
		Assert.assertEquals("getTempFactor", "0.00", record.getTempFactor());
		Assert.assertNull("getRecId", record.getRecId());
		Assert.assertEquals("getSegId", "", record.getSegId());
		Assert.assertEquals("getElement", "", record.getElement());
		Assert.assertEquals("getCharge", "", record.getCharge());
	}

	/**
	 * This tests parses a pdb atom record.
	 */
	@Test
	public void testPdbRecord3() {
		PdbAtomRecord record = new PdbAtomRecord(
				"ATOM     14  C     C    14      16.252   3.241   0.158  0.00  0.00           C");

		System.out.println("getRecname " + record.getRecname());
		System.out.println("getSerial " + record.getSerial());
		System.out.println("getAtom " + record.getAtom());
		System.out.println("getAltLocation " + record.getAltLocation());
		System.out.println("getResName " + record.getResName());
		System.out.println("getChainId " + record.getChainId());
		System.out.println("getSeqNo " + record.getSeqNo());
		System.out.println("getY " + record.getY());
		System.out.println("getZ " + record.getZ());
		System.out.println("getOccupancy " + record.getOccupancy());
		System.out.println("getTempFactor " + record.getTempFactor());
		System.out.println("getRecId " + record.getRecId());
		System.out.println("getSegId " + record.getSegId());
		System.out.println("getElement " + record.getElement());
		System.out.println("getCharge " + record.getCharge());

		Assert.assertEquals("ATOM", record.getRecname());
		Assert.assertEquals("14", record.getSerial());
		Assert.assertEquals("C", record.getAtom());
		Assert.assertEquals("", record.getAltLocation());
		Assert.assertEquals("C", record.getResName());
		Assert.assertEquals("", record.getChainId());
		Assert.assertEquals("14", record.getSeqNo());
		Assert.assertEquals("16.252", record.getX());
		Assert.assertEquals("3.241", record.getY());
		Assert.assertEquals("0.158", record.getZ());
		Assert.assertEquals("0.00", record.getOccupancy());
		Assert.assertEquals("0.00", record.getTempFactor());
		Assert.assertEquals("", record.getSegId());
		Assert.assertEquals("C", record.getElement());
		Assert.assertEquals(null, record.getCharge());
		Assert.assertEquals(null, record.getRecId());
	}
}