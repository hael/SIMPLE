/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.model;

import org.openscience.cdk.Atom;
import org.openscience.cdk.AtomContainer;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.interfaces.IAtom;

import javax.vecmath.Point3d;

/**
 * Data conversion utilities for the data exchange between OpenMolGRID and CDK.
 *
 * @author Sulev Sild
 */
public class Convertor {

	/**
	 * Converst CStructure to CDK's AtomContainer
	 *
	 * @param str is a  CStructure object to be converted
	 * @return a converted AtomContainer object
	 */
	public static AtomContainer convert(CStructure str)
	{
		if (str == null)
			return null;
		
		AtomContainer ac = new AtomContainer();
		
		// copy atoms 
		for (int i=1; i<=str.getAtomCount(); ++i)
		{
			ac.addAtom(convert(str.getAtom(i)));
		}
		
		// copy bonds
		for (int i=1; i<=str.getBondCount(); ++i)
		{
			CBond b = str.getBond(i);
			ac.addBond(b.getAtom1Id()-1, b.getAtom2Id()-1, b.getType());
		}
		
		return ac;
	}

	/**
	 * Converts CAtom to CDK's Atom
	 *
	 * @param atom is a CAtom object to be converted
	 * @return a converted Atom object
	 */
	public static Atom convert(CAtom atom)
	{
		Atom cdkAtom = new Atom(atom.getElement());
		cdkAtom.setPoint3d(new Point3d(atom.getX(), atom.getY(), atom.getZ()));
		cdkAtom.setFormalCharge(atom.getFormalCharge());
		cdkAtom.setID(String.valueOf(atom.getId()-1));
		return cdkAtom;
	}
	
	public static Molecule convertMolecule(CStructure s) {
		return new Molecule(convert(s));
	}

	/**
	 * Converts CDK's AtomContainer object to CStructure. 
	 *
	 * @param ac is an AtomContainer object to be converted
	 * @param coord2D is true when 2D coordinates from AtomContainer object are used
	 * @return str is a converted CStructure object
	 */
	public static CStructure convert(AtomContainer ac, boolean coord2D)
	{
		CStructure str = new CStructure();

		// copy atoms
		for (int i=0; i<ac.getAtomCount(); ++i) {
			str.addAtom(convert(ac.getAtom(i), i+1, coord2D));
		}

		// copy bonds
		for (int i=0; i<ac.getBondCount(); ++i) {
			CAtom a1 = str.getAtom(ac.getAtomNumber(ac.getBond(i).getAtom(0)) + 1);
			CAtom a2 = str.getAtom(ac.getAtomNumber(ac.getBond(i).getAtom(1)) + 1);
			CBond b = new CBond(a1, a2, (int)ac.getBond(i).getOrder());
			str.addBond(b);
		}

		return str;
	}

	/**
	 * Converst CDK's Atom to CAtom
	 *
	 * @param atom is a CDK's Atom object to be converted
	 * @param id is an atom id that will be assigned to new CAtom object
	 * @param coord2D is true when 2D coordinates from Atom object are used
	 * @return a converted CAtom object
	 */
	public static CAtom convert(IAtom atom, int id, boolean coord2D)
	{
		double x, y, z;

		if (coord2D)
		{
			x = atom.getPoint2d().x;
			y = atom.getPoint2d().y;
			z = 0.0;
		} else
		{
			x = atom.getPoint3d().x;
			y = atom.getPoint3d().y;
			z = atom.getPoint3d().z;
		}
		CAtom result = new CAtom(atom.getSymbol(), x, y, z, id);
		result.setFormalCharge(atom.getFormalCharge());
		return result;
	}

}
