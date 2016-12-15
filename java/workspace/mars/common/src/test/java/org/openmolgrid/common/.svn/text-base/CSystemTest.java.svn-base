package org.openmolgrid.common;

import java.util.ArrayList;
import java.util.List;

import junit.framework.Assert;

import org.junit.Test;
import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.CSystem;

/**
 * 
 * @author Stefan Bozic
 */
public class CSystemTest {

	/**
	 * The methods creates a {@link CSystem}
	 */
	@Test
	public void createEmptySystem() {
		CSystem system = new CSystem();
		Assert.assertNotNull(system);
	}
	
	/**
	 * The methods creates a {@link CSystem}
	 */
	@Test
	public void createSystemWithStructureList() {
		CStructure structure = new CStructure();
		CStructure structure2 = new CStructure();

		CAtom atom1 = new CAtom("H", 1.0, 1.0, 1.0, 1);
		atom1.setPartialCharge(0.9);
		structure.addAtom(atom1);
		structure2.addAtom(atom1);

		CAtom atom2 = new CAtom("H", -1.0, -1.0, -1.0, 2);
		atom2.setPartialCharge(0.567);
		structure.addAtom(atom2);
		structure2.addAtom(atom2);

		List<CStructure> structureList = new ArrayList<CStructure>();
		structureList.add(structure);
		structureList.add(structure2);
		
		CSystem system = new CSystem(structureList);
		Assert.assertNotNull(system);
		Assert.assertEquals(2, system.getStructures().size());
	}

	/**
	 * The methods creates a CML {@link String} based on a {@link CStructure}
	 */
	@Test
	public void createCml() {
		CStructure structure = new CStructure();
		structure.setId("m1");
		
		CAtom atom1 = new CAtom("H", 1.0, 1.0, 1.0, 1);
		structure.addAtom(atom1);
		
		CAtom atom2 = new CAtom("H", -1.0, -1.0, -1.0, 2);
		structure.addAtom(atom2);		
		
		CSystem system = new CSystem();
		system.addStructure(structure);
		
		String cml = system.writeCML2String();
		System.out.println(cml);
		Assert.assertNotNull(cml);
		Assert.assertTrue(cml.contains("<cml"));
		Assert.assertTrue(cml.contains("<molecule"));
	}	
	
	/**
	 * Test for {@link CSystem#getStructurebyId(String)}
	 */
	@Test
	public void testGetStructureById() {
		CStructure structure = new CStructure();
		structure.setId("m1");
		CStructure structure2 = new CStructure();
		structure2.setId("m2");
		CStructure structure3 = new CStructure();
		structure3.setId("m3");

		
		List<CStructure> structureList = new ArrayList<CStructure>();
		structureList.add(structure);
		structureList.add(structure2);
		structureList.add(structure3);
		
		CSystem system = new CSystem(structureList);		
		
		Assert.assertNotNull(system);
		Assert.assertEquals(3, system.getStructures().size());
		
		
		CStructure structureById = system.getStructurebyId("m2");
		
		Assert.assertEquals(structure2, structureById);
	}
}
