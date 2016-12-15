/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.util.*;
import java.io.PrintWriter;


/**
 * This class provides mapping between property id and its metadata.
 * 
 * @author Sulev Sild
 */
public class CPropertyMetaMap
{
	private HashMap<String, CPropertyMeta> data;
	
	public CPropertyMetaMap()
	{
		data = new HashMap<String, CPropertyMeta>();
	}

	/**
	 * Adds mapping between the property Id and its metadata.
	 *
	 * @param propId is a property Id
	 * @param meta is a metadata
	 */
	public void put(String propId, CPropertyMeta meta)
	{
		data.put(propId, meta);
	}

	/**
	 * Check whether metadata is defined for a given property Id.
	 *
	 * @param propId is a property Id.
	 * @return true if metadata is defined
	 */
	public boolean hasMetaData(String propId)
	{
		return data.get(propId) != null;
	}
	
	public CPropertyMeta getMetaData(String propId)
	{
		return data.get(propId);
	}
	
	/**
	 * Returns the number of properties in the property meta data map.
	 * 
	 * @return number of properties
	 */
	public int getNumberOfProperties()
	{
		return data.size();
	}
	
	/**
	 * Returns property name for the given property id.
	 *
	 * @param propId a property id
	 * @return property name or null if not defined
	 */
	public String getName(String propId)
	{
		CPropertyMeta m = (CPropertyMeta)data.get(propId);
		if (m == null)
			return null;
		else
			return m.getAttribute(CPropertyMeta.NAME);
	}
	
	/**
	 * Returns Iterator for all the descriptor Ids that are stored in this map.
	 *
	 * @return Iterator
	 */
	public Iterator<String> getIterator()
	{
		return data.keySet().iterator();
	}
	
	public void mergeMetaData(CPropertyMetaMap meta)
	{
		Iterator<String> it = meta.getIterator();
		while (it.hasNext()) {
			String id = (String) it.next();
			if (!data.containsKey(id)) {
				put(id, meta.getMetaData(id));
			}
		}
	}


	public void writeCDR(PrintWriter pw)
	{
		Iterator<String> it = data.keySet().iterator();
		if (!it.hasNext())
			return;
		
		while (it.hasNext())
		{
			Object id = it.next();
			Object val = data.get(id);
			pw.println("    <property id=\"" + id + "\">");
			((CPropertyMeta)val).writeCDR(pw);
			pw.println("    </property>");
		}
	}
}
