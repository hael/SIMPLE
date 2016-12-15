/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.model;

import java.io.Serializable;

/**
 * Representation of atoms.
 * 
 * @author Andre Lomaka
 * @author Sulev Sild
 * @author Stefan Bozic
 */
public class CAtom implements Serializable {

	/**
	 * generated serialVersionUID
	 */
	private static final long serialVersionUID = 8709229821925686333L;

	/**
	 * the id of the atom
	 */
	private int id;

	/**
	 * the xyz coordinates
	 */
	private double X = 0.0, Y = 0.0, Z = 0.0;

	/**
	 * The element type.
	 * 
	 * @see CConstants
	 */
	private int element;

	/**
	 * <p>
	 * In chemistry, a formal charge (FC) is the charge assigned to an atom in a molecule, assuming that electrons in a
	 * chemical bond are shared equally between atoms, regardless of relative electronegativity.
	 * </p>
	 * 
	 * FC = V - N - B/2<br>
	 * V = number of valence electrons<br>
	 * N = number of non-bonding valence electrons on this atom in the molecule<br>
	 * B = is the total number of electrons shared in covalent bonds with other atoms in the molecule
	 * 
	 * @see <a href="http://en.wikipedia.org/wiki/Formal_charge">Wikipedia entry</a>
	 */
	private int formalCharge = 0;

	/**
	 * Partial charges are created due to the asymmetric distribution of electrons in chemical bonds. The resulting
	 * partial charges are a property only of zones within the distribution, and not the assemblage as a whole. For
	 * example, chemists often choose to look at a small space surrounding the nucleus of an atom: When an electrically
	 * neutral atom bonds chemically to another neutral atom that is more electronegative, its electrons are partially
	 * drawn away. This leaves the region about that atom's nucleus with a partial positive charge, and it creates a
	 * partial negative charge on the atom to which it is bonded.
	 * 
	 * @see <a href="http://en.wikipedia.org/wiki/Partial_charge">Wikipedia entry</a>
	 */
	private double partialCharge = 0.0;

	/**
	 * Copy Constructor
	 * 
	 * @param a the atom to copy
	 */
	public CAtom(CAtom a) {
		id = a.id;
		X = a.X;
		Y = a.Y;
		Z = a.Z;
		element = a.element;
		formalCharge = a.formalCharge;
		partialCharge = a.partialCharge;
	}

	/**
	 * Constructor
	 * 
	 * @param element the chimcal element
	 * @param x the x coordinate
	 * @param y the y coordinate
	 * @param z the z coordinate
	 * @param id the id of the atom
	 */
	public CAtom(String element, double x, double y, double z, int id) {
		X = x;
		Y = y;
		Z = z;
		this.element = CConstants.getElementByString(element);
		this.id = id;
	}

	/**
	 * Returns the id for this atom
	 * 
	 * @return the id for this atom
	 */
	public int getId() {
		return id;
	}

	/**
	 * Sets the formal charge for this atom.
	 * 
	 * @param fc
	 */
	public void setFormalCharge(int fc) {
		formalCharge = fc;
	}

	/**
	 * Returns the formal charge for this atom
	 * 
	 * @return the formal charge for this atom
	 */
	public int getFormalCharge() {
		return formalCharge;
	}

	/**
	 * Checks if the atom is from a specific element type
	 * 
	 * @param el the element type to check for
	 * 
	 * @return <code>true</code> if the atom is from a specific element type
	 */
	public boolean isElement(int el) {
		return el == element;
	}

	/**
	 * Checks if the atom is from a specific element type
	 * 
	 * @param el the element type to check for
	 * 
	 * @return <code>true</code> if the atom is from a specific element type
	 */
	public boolean isElement(CConstants el) {
		return el.ordinal() == element;
	}

	/**
	 * Sets the element type for this atom
	 * 
	 * @param el the value to set
	 */
	public void setElement(int el) {
		element = el;
	}

	/**
	 * Sets the element type for this atom
	 * 
	 * @param el the value to set
	 */
	public void setElement(CConstants el) {
		element = el.ordinal();
	}

	/**
	 * Returns the element type of this atom as string
	 * 
	 * @return the element type of this atom as string
	 */
	public String getElement() {
		return CConstants.getStringByElement(element);
	}

	/**
	 * Returns the atomic number of this atom
	 * 
	 * @return the atomic number of this atom
	 */
	public int getAtomicNumber() {
		return element;
	}

	/**
	 * Returns the x coordinate of this atom
	 * 
	 * @return the x coordinate of this atom
	 */
	public double getX() {
		return X;
	}

	/**
	 * Returns the y coordinate of this atom
	 * 
	 * @return the y coordinate of this atom
	 */
	public double getY() {
		return Y;
	}

	/**
	 * Returns the z coordinate of this atom
	 * 
	 * @return the z coordinate of this atom
	 */
	public double getZ() {
		return Z;
	}

	/**
	 * @param x the x to set
	 */
	public void setX(double x) {
		X = x;
	}

	/**
	 * @param y the y to set
	 */
	public void setY(double y) {
		Y = y;
	}

	/**
	 * @param z the z to set
	 */
	public void setZ(double z) {
		Z = z;
	}

	/**
	 * @return the partialCharge
	 */
	public double getPartialCharge() {
		return partialCharge;
	}

	/**
	 * @param partialCharge the partialCharge to set
	 */
	public void setPartialCharge(double partialCharge) {
		this.partialCharge = partialCharge;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "CAtom [X=" + X + ", Y=" + Y + ", Z=" + Z + ", element=" + CConstants.getStringByElement(element) + ", formalCharge=" + formalCharge
				+ ", id=" + id + ", partialCharge=" + partialCharge + "]";
	}	
}
