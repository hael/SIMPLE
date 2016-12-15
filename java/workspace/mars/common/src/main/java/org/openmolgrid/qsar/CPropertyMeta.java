/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.util.ArrayList;
import java.io.PrintWriter;

import org.openmolgrid.util.CUtil;


/**
 * This class in intended for holding metadata about properties.
 * 
 * @author Sulev Sild
 */
public class CPropertyMeta
{
	// supported attributes for properties
	public static final int NAME = 0;
	private static final int numberOfAttributes = 1;
	
	private ArrayList<String> data;

	/**
	 * Default constructor.
	 */
	public CPropertyMeta()
	{
		data = new ArrayList<String>(numberOfAttributes);
		for (int i=0; i<numberOfAttributes; ++i)
			data.add(null);
	}

	/**
	 * Copy constructor.
	 *
	 * @param meta is a CPropertyMeta object
	 */
	public CPropertyMeta(CPropertyMeta meta)
	{
		this();
		for (int i=0; i<numberOfAttributes; ++i)
			setAttribute(i, meta.getAttribute(i));
	}

	/**
	 * Sets attribute to given value.
	 *
	 * @param attr is the attribute ID
	 * @param value is a value of the attribute
	 */
	public void setAttribute(int attr, String value)
	{
		data.set(attr, value);
	}
	
	/**
	 * Checks whether the attribute has a value.
	 *
	 * @param attr is the attribute ID
	 * @return true if value is defined, false otherwise
	 */
	public boolean hasAttribute(int attr)
	{
		return data.get(attr) != null;
	}

	/**
	 * Returns the attribute value.
	 *
	 * @param attr is the attribute ID
	 * @return a String with the attribute value
	 */
	public String getAttribute(int attr)
	{
		return data.get(attr);
	}

	public void writeCDR(PrintWriter pw)
	{
		if (hasAttribute(NAME))
			pw.println("      <name>" + CUtil.escapeXMLString(getAttribute(NAME)) + "</name>");
	}

}
