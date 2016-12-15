/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.io.InputStream;

import org.openmolgrid.format.xml.XMLGenericParser;
import org.openmolgrid.model.ChemLibException;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

/**
 * This class implements parser for reading calculation results.
 *
 * Sulev Sild
 */
public class CalculationResultReader extends XMLGenericParser{
	private CalculationResultConsumer crc;
	private CalculationResult curResult;
	
	public CalculationResultReader(String inputFile, CalculationResultConsumer crc) {
		super(inputFile, false, null);
		this.crc = crc;
		curResult = null;
	}

	public CalculationResultReader(InputStream inputStream, CalculationResultConsumer crc) {
		super(new InputSource(inputStream), false, null);
		this.crc = crc;
		curResult = null;
	}
	
	/**
	 * Starts results file parser.
	 */
	public void produceData() throws ChemLibException {
		try {
			startParsing();
		} catch (Exception e) {
			throw new ChemLibException(e.getMessage());
		}
	}

	// SAX Document Handler methods
	
	/**
	 * Handle start of a tag.
	 *
	 * @param ns is an XML namespace
	 * @param lName is a tag name
	 * @param qName ia a fully qualified XML tag name
	 * @param atts are tag attributes
	 * @exception SAXException if an error occurs
	 */
	public void startElement(String ns, String lName, String qName, Attributes atts) throws SAXException {
		if ("result".equals(lName)) {
			String id = atts.getValue("structureId");
			curResult = new CalculationResult(id);
		} else if ("status".equals(lName)) {
			collectCharacterData();
		} else if ("error".equals(lName)) {
			collectCharacterData();
		} else if ("output".equals(lName)) {
			collectCharacterData();
		}
	}
	
	/**
	 * Handle end of a tag.
	 *
	 * @param ns is an XML namespace
	 * @param lName is a tag name
	 * @param qName is a fully qualified XML tag name
	 * @exception SAXException if an error occurs
	 */
	public void endElement(String ns, String lName, String qName) throws SAXException {
		if ("result".equals(lName)) {
			crc.processResult(curResult);
		} else if ("status".equals(lName)) {
			curResult.setStatus(getCharacterData());
		} else if ("error".equals(lName)) {
			curResult.setError(getCharacterData());
		} else if ("output".equals(lName)) {
			curResult.setOutput(getCharacterData());
		}
	}

}
