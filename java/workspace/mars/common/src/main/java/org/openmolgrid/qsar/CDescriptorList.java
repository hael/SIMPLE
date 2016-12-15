/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.io.*;
import java.util.*;


public class CDescriptorList implements Serializable, CAttribute {
	private ArrayList<CDescriptor> descriptors = new ArrayList<CDescriptor>();

	public CDescriptorList() {
	}

	public CDescriptorList(CDescriptorList d) {
		for (int i = 0; i < d.descriptors.size(); i++) {
			descriptors.add(new CDescriptor(d.descriptors.get(i)));
		}
	}

	public void addDescriptor(CDescriptor d) {
		descriptors.add(d);
	}

	public void addDescriptor(String id, double value) {
		descriptors.add(new CDescriptor(id, value));
	}

	public void addDescriptor(String id, int pos, double value) {
		descriptors.add(new CDescriptor(id, pos, value));
	}

	public CDescriptor getDescriptor(int i) {
		return descriptors.get(i);
	}

	public CDescriptor getDescriptorById(String id) {
		for (int i = 0; i < descriptors.size(); i++) {
			CDescriptor d = descriptors.get(i);
			if (d.getId().equals(id))
				return d;
		}
		return null;
	}

	public int getNumberOfDescriptors() {
		return descriptors.size();
	}

	public CAttribute copy() {
		return new CDescriptorList(this);
	}
	
	public void writeCDR(PrintWriter pw) {
		pw.println("    <descriptorList>");
		for (int i = 0; i < descriptors.size(); i++) {
			CDescriptor d = descriptors.get(i);
			pw.print("      <descriptor id=\"" + d.getId() + "\"");
			if (d.getPos() != 0)
				pw.print(" pos=\"" + d.getPos() + "\"");
			pw.println(" value=\"" + d.getValue() + "\"/>");
		}
		pw.println("    </descriptorList>");
	}

}
