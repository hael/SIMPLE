/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.util.*;
import java.io.PrintWriter;

/**
 * This class provides mapping between descriptor id and its metadata.
 * 
 * @author Sulev Sild
 */
public class CDescriptorMetaMap implements Iterable<String> {
	private LinkedHashMap<String, CDescriptorMeta> data;

	public CDescriptorMetaMap() {
		data = new LinkedHashMap<String, CDescriptorMeta>();
	}

	/**
	 * Adds mapping between the descriptor Id and its meta data.
	 * 
	 * @param descId - descriptor Id
	 * @param meta - descriptor meta data
	 */
	public void put(String descId, CDescriptorMeta meta) {
		data.put(descId, meta);
	}

	/**
	 * Returns the number of descriptors in the descriptor meta data map.
	 * 
	 * @return number of descriptors
	 */
	public int getNumberOfDescriptors() {
		return data.size();
	}

	/**
	 * Check whether metadata is defined for a given descriptor Id.
	 * 
	 * @param descId - descriptorId
	 * @return true if metadata is defined
	 */
	public boolean hasMetaData(String descId) {
		return data.get(descId) != null;
	}

	public CDescriptorMeta getMetaData(String descId) {
		return data.get(descId);
	}

	/**
	 * Returns descriptor name for the given descriptor id.
	 * 
	 * @param descId - descriptor id
	 * @return descriptor name or null if not defined
	 */
	public String getName(String descId) {
		CDescriptorMeta m = data.get(descId);
		if (m == null)
			return null;
		else
			return m.getAttribute(CDescriptorMeta.NAME);
	}

	/**
	 * Returns descriptor calculation application name for the given descriptor
	 * id.
	 * 
	 * @param descId - descriptor id
	 * @return application name or null if not defined
	 */
	public String getAppName(String descId) {
		CDescriptorMeta m = data.get(descId);
		if (m == null)
			return null;
		else
			return m.getAttribute(CDescriptorMeta.APPNAME);
	}

	/**
	 * Returns descriptor calculation application version for the given
	 * descriptor id.
	 * 
	 * @param descId - descriptor id
	 * @return application version or null if not defined
	 */
	public String getAppVersion(String descId) {
		CDescriptorMeta m = data.get(descId);
		if (m == null)
			return null;
		else
			return m.getAttribute(CDescriptorMeta.APPVERSION);
	}

	/**
	 * Returns descriptor definition for the given descriptor id.
	 * 
	 * @param descId - descriptor id
	 * @return descriptor definition or null if not defined
	 */
	public String getDescDefinition(String descId) {
		CDescriptorMeta m = data.get(descId);
		if (m == null)
			return null;
		else
			return m.getAttribute(CDescriptorMeta.DEFINITION);
	}

	/**
	 * Returns literary reference for the given descriptor id.
	 * 
	 * @param descId - descriptor id
	 * @return descriptor type or null if not defined
	 */
	public String getCitation(String descId) {
		CDescriptorMeta m = data.get(descId);
		if (m == null)
			return null;
		else
			return m.getAttribute(CDescriptorMeta.CITATION);
	}

	/**
	 * Returns descriptor type for the given descriptor id.
	 * 
	 * @param descId - descriptor id
	 * @return descriptor type or null if not defined
	 */
	public String getType(String descId) {
		CDescriptorMeta m = data.get(descId);
		if (m == null)
			return null;
		else
			return m.getAttribute(CDescriptorMeta.TYPE);
	}

	/**
	 * Returns descriptor subtype for the given descriptor id.
	 * 
	 * @param descId
	 *            a descriptor id
	 * @return descriptor subtype or null if not defined
	 */
	public String getSubtype(String descId) {
		CDescriptorMeta m = data.get(descId);
		if (m == null)
			return null;
		else
			return m.getAttribute(CDescriptorMeta.SUBTYPE);
	}

	/**
	 * Returns descriptor category for the given descriptor id.
	 * 
	 * @param descId - descriptor id
	 * @return descriptor type or null if not defined
	 */
	public String getCategory(String descId) {
		CDescriptorMeta m = data.get(descId);
		if (m == null)
			return null;
		else
			return m.getAttribute(CDescriptorMeta.CATEGORY);
	}

	public void mergeMetaData(CDescriptorMetaMap meta) {
		for (String id: meta) {
			if (!data.containsKey(id)) {
				put(id, meta.getMetaData(id));
			}
		}
	}

	/**
	 * Returns Iterator for all the descriptor Ids that are stored in this map.
	 * 
	 * @return Iterator
	 */
	public Iterator<String> iterator() {
		return data.keySet().iterator();
	}

	public void writeCDR(PrintWriter pw) {
		for (String id: this) {
			pw.println("    <descriptor id=\"" + id + "\">");
			data.get(id).writeCDR(pw);
			pw.println("    </descriptor>");
		}
	}

}
