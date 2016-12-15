/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.io.Writer;
import java.io.IOException;
import java.io.BufferedWriter;
import java.io.File;

import org.openmolgrid.util.CUtil;


/**
 * This class is used to store calculation results for chemical structures.
 *
 * @author Sulev Sild
 */
public class CalculationResult 
{
	private String ID;
	private String status;
	private String error;
	private String output;
	
	/**
	 * Creates a new CalculationResult instance.
	 *
	 * @param ID is the structure ID
	 */
	public CalculationResult(String ID)
	{
		this.ID = ID;
	}

	/**
	 * Sets status of the calculation.
	 *
	 * @param status
	 */
	public void setStatus(String status)
	{
		this.status = status;
	}

	/**
	 * Sets error message
	 *
	 * @param error
	 */
	public void setError(String error)
	{
		this.error = error;
	}

	/**
	 * Sets contents of the output file to the CalculationResult object.
	 *
	 * @param outputData
	 */
	public void setOutput(String outputData)
	{
		output = outputData;
	}


	/**
	 * Reads contents from the output file and sets it to the CalculationResult
	 * object.
	 *
	 * @param fileName
	 */
	public void setOutput(File fileName)
	{
		output = CUtil.file2String(fileName);
	}

	public String getId()
	{
		return ID;
	}

	public String getStatus()
	{
		return status;
	}

	public String getError()
	{
		return error;
	}

	public String getOutput()
	{
		return output;
	}
	

	/**
	 * Returns true, if the calculation was successful.
	 *
	 * @return boolean
	 */
	public boolean isSuccessful()
	{
		if (status.equals("OK"))
			return true;
		else
			return false;
	}

	/**
	 * Writes an XML representation of the CalculationResult object to the
	 * given writer.
	 *
	 * @param writer a <code>Writer</code> value
	 */
	public void writeCDR(Writer outWriter)
		throws IOException
	{
		BufferedWriter writer;
		if (outWriter instanceof BufferedWriter)
			writer = (BufferedWriter)outWriter;
		else
			writer = new BufferedWriter(outWriter);
		
		writer.write("<result structureId=\"");
		writer.write(getId());
		writer.write("\">");

		if (getStatus() != null)
		{
			writer.write("  <status>");
			writer.write(getStatus());
			writer.write("</status>");
			writer.newLine();
		}	

		if (getError() != null) 
		{
			writer.write("  <error>");
			writer.write(getError());
			writer.write("</error>");
			writer.newLine();
		}

		if (getOutput() != null)
		{
			writer.write("  <output><![CDATA[");
			writer.write(output);
			writer.write("]]></output>");
			writer.newLine();
		}		

		writer.write("</result>");
		writer.newLine();
		writer.flush();
	}

}
