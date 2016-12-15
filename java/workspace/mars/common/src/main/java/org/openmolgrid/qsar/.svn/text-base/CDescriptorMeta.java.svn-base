/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.util.ArrayList;
import java.io.PrintWriter;

import org.openmolgrid.util.CUtil;


/**
 * This class is intended for holding metadata about descriptors.
 * 
 * @author Sulev Sild
 */
public class CDescriptorMeta
{
	// supported attributes for descriptors
	public static final int NAME = 0;
	public static final int TYPE = 1;
	public static final int SUBTYPE = 2;
	public static final int CATEGORY = 3;
	public static final int APPNAME = 4;
	public static final int APPVERSION = 5;
	public static final int DEFINITION = 6;
	public static final int CITATION = 7;
	private static final int numberOfAttributes = 8;
	
	private ArrayList<String> data;

	/**
	 * Default constructur.
	 */
	public CDescriptorMeta()
	{
		data = new ArrayList<String>(numberOfAttributes);
		for (int i=0; i<numberOfAttributes; ++i)
			data.add(null);
	}

	/**
	 * Copy constructor.
	 *
	 * @param meta is a CDescriptorMeta object
	 */
	public CDescriptorMeta(CDescriptorMeta meta)
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
		return (String)data.get(attr);
	}

	public void writeCDR(PrintWriter pw)
	{
		if (hasAttribute(NAME))
			pw.println("      <name>" + CUtil.escapeXMLString(getAttribute(NAME)) + "</name>");
		if (hasAttribute(TYPE))
			pw.println("      <type>" + CUtil.escapeXMLString(getAttribute(TYPE)) + "</type>");
		if (hasAttribute(SUBTYPE))
			pw.println("      <subType>" + CUtil.escapeXMLString(getAttribute(SUBTYPE)) + "</subType>");
		if (hasAttribute(CATEGORY))
			pw.println("      <category>" + CUtil.escapeXMLString(getAttribute(CATEGORY)) + "</category>");
		if (hasAttribute(APPNAME))
			pw.println("      <appName>" + CUtil.escapeXMLString(getAttribute(APPNAME)) + "</appName>");
		if (hasAttribute(APPVERSION))
			pw.println("      <appVersion>" + CUtil.escapeXMLString(getAttribute(APPVERSION)) + "</appVersion>");
		if (hasAttribute(DEFINITION))
			pw.println("      <definition>" + CUtil.escapeXMLString(getAttribute(DEFINITION)) + "</definition>");
		if (hasAttribute(CITATION))
			pw.println("      <citation>" + CUtil.escapeXMLString(getAttribute(CITATION)) + "</citation>");
	}

}
