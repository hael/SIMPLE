/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.io.*;

import org.openmolgrid.format.slf.StructureWriter;

/**
 * This class is intended for creating fragment list files.
 *
 * @author Sulev Sild
 */
public class FragmentWriter extends StructureWriter
{
	public FragmentWriter(String outFileName)
	{
		super(outFileName);
	}

	public FragmentWriter(File outFile)
	{
		super(outFile);
	}
	
	protected void writeHeader()
	{
		outWriter.println("<?xml version=\"1.0\"?>");
		outWriter.println("<fragmentList>");
	}

	protected void writeFooter()
	{
		outWriter.println("</fragmentList>");
	}

}
