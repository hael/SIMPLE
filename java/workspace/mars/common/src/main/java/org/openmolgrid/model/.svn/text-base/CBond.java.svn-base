/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.model;

import java.io.*;

public class CBond implements Serializable {

	public enum Stereo {
		NONE, // no stereo information about the bond
		UP,   // wide end of the bond is above the plane (i.e. wedge bond)
		DOWN; // wide end of the bond is below the plane (i.e. hatch bond)
	}
	
	private int order;
	private Stereo stereo = Stereo.NONE;
	private CAtom atom1, atom2;

	public CBond(CAtom bnda1, CAtom bnda2, int bondOrder) {
		atom1 = bnda1;
		atom2 = bnda2;
		order = bondOrder;
	}

	public CBond(CAtom bnda1, CAtom bnda2, int bondOrder, Stereo bondStereo) {
		atom1 = bnda1;
		atom2 = bnda2;
		order = bondOrder;
		stereo = bondStereo;
	}

	public CAtom getAtom(int i) {
		if (i == 1)
			return atom1;
		else if (i == 2)
			return atom2;
		else
			throw new IllegalArgumentException("Expected 1 or 2 for CBond.getAtom(idx)");
	}

	public void setStereo(Stereo bondStereo) {
		stereo = bondStereo;
	}

	public Stereo getStereo() {
		return stereo;
	}

	public int getType() {
		return order;
	}

	public int getAtom1Id() {
		return atom1.getId();
	}

	public int getAtom2Id() {
		return atom2.getId();
	}

}
