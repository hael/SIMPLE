/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.io.*;


/**
 * This class is used for recording the status information about the 2D to 3D
 * conversion of molecular structures.<p>
 *
 * @author Sulev Sild
 */
public class ErrorWriter
{
	private BufferedWriter writer;
	
	/**
	 * Creates a new ErrorWriter instance.
	 *
	 * @param fileName is a file name for error log file
	 * @exception IOException if I/O error occurs
	 */
	public ErrorWriter(String fileName) throws IOException
	{
		writer = new BufferedWriter(new FileWriter(fileName));
		writer.write("<?xml version=\"1.0\"?>");
		writer.newLine();
		writer.write("<calculationLog>");
		writer.newLine();
	}

	/**
	 * Adds record about the calculation to the log file.
	 *
	 * @param structureID is a structure ID of the molecule
	 * @param res a <code>CalculationResult</code> value
	 */
	public void addResult(String structureID, CalculationResult res)
		throws IOException
	{
		res.writeCDR(writer);
	}


	/**
	 * Adds record about the successful calculation to the log file.
	 *
	 * @param stuctureID is a structure ID of the molecule
	 */
	public void addSuccess(String structureID) throws IOException
	{
		CalculationResult res = new CalculationResult(structureID);
		res.setStatus("OK");
		res.writeCDR(writer);
	}


	/**
	 * Add record about the failed calculation to the log file.
	 *
	 * @param structureID is a structure ID of a molecule
	 * @param message is a error message
	 */
	public void addFailure(String structureID, String message) 
		throws IOException
	{
		CalculationResult res = new CalculationResult(structureID);
		res.setStatus("FAILED");
		res.setError(message);
		res.writeCDR(writer);
	}


	/**
	 * Closes ErrorWriter object.
	 */
	public void close() throws IOException
	{
		writer.write("</calculationLog>");
		writer.close();
	}
	
}
