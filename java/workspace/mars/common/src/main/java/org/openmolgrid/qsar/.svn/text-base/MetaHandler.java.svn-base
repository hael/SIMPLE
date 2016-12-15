/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import org.xml.sax.*;
import org.xml.sax.helpers.DefaultHandler;
import java.lang.StringBuffer;


/**
 * SAX handler for metadata. 
 *
 * @author Sulev Sild
 */
public class MetaHandler extends DefaultHandler
{
	// container for descriptor metadata
	private CDescriptorMetaMap descMetaMap = new CDescriptorMetaMap();
 
	// container for property metadata
	private CPropertyMetaMap propMetaMap = new CPropertyMetaMap();

	// current metadata entries for descriptor and property
	private CDescriptorMeta currentDescMeta;
	private CPropertyMeta currentPropMeta;
	private String currentId;

	private String parentElement;

	// Buffer to collect character data. 
	private StringBuffer characters = new StringBuffer();

	// Flag that enables/disables collection of character data
	private boolean charDataCollected = false;

	// Flag that enables/disables meta handler
	private boolean handlerActive = false;


	/**
	 * Returns true if MetaHandler is currently parsing data.
	 *
	 * @return true if handler is active
	 */
	public boolean isActive()
	{
		return handlerActive;
	}	

	public void startElement(String ns, String lName, String qName, Attributes atts)
		throws SAXException
	{
		if ("metadata".equals(lName))
		{
			handlerActive = true;
		} else if ("descriptor".equals(lName))
		{
			currentId = atts.getValue("id");
			currentDescMeta = new CDescriptorMeta();
			parentElement = lName;
		} else if ("property".equals(lName))
		{
			currentId = atts.getValue("id");
			currentPropMeta = new CPropertyMeta();
			parentElement = lName;
		} else if ("name".equals(lName))
		{
			collectCharacterData();
		} else if ("type".equals(lName))
		{
			collectCharacterData();
		} else if ("subtype".equals(lName))
		{
			collectCharacterData();
		} else if ("category".equals(lName))
		{
			collectCharacterData();
		}
	}
	

	public void endElement(String ns, String lName, String qName)
		throws SAXException
	{
		if ("metadata".equals(lName))
		{
			handlerActive = false;
		} else if ("descriptor".equals(lName))
		{
			descMetaMap.put(currentId, currentDescMeta);
		} else if ("property".equals(lName))
		{
			propMetaMap.put(currentId, currentPropMeta);
		} else if ("name".equals(lName))
		{
			if ("descriptor".equals(parentElement))
				currentDescMeta.setAttribute(CDescriptorMeta.NAME, getCharData());
			else if ("property".equals(parentElement))
				currentPropMeta.setAttribute(CPropertyMeta.NAME, getCharData());
		} else if ("type".equals(lName))
		{
			currentDescMeta.setAttribute(CDescriptorMeta.TYPE, getCharData());
		} else if ("subtype".equals(lName))
		{
			currentDescMeta.setAttribute(CDescriptorMeta.SUBTYPE, getCharData());
		} else if ("category".equals(lName))
		{
			currentDescMeta.setAttribute(CDescriptorMeta.CATEGORY, getCharData());
		}		
	}


	public void characters(char buf[], int offset, int len)
		throws SAXException
	{
		if (charDataCollected)
			characters.append(buf, offset, len);
	}
	

	/**
	 * Returns character data that has been collected and stops collection
	 * of character data.
	 *
	 * @return String with collected data
	 */
	private String getCharData()
	{
		charDataCollected = false;
		return characters.toString();
	}

	/**
	 * Instructs SAX parser to collect character data.
	 */
	private void collectCharacterData()
	{
		charDataCollected = true;
		characters.setLength(0);
	}	

	/**
	 * Returns reference to descriptor metadata that has been collected by
	 * the MetaHandler.
	 *
	 * @return CDescriptorMetaMap
	 */
	public CDescriptorMetaMap getDescriptorMetaData()
	{
		return descMetaMap;
	}

	/**
	 * Returns reference to property metadata that has been collected by
	 * the MetaHandler.
	 *
	 * @return CPropertyMetaMap
	 */
	public CPropertyMetaMap getPropertyMetaData()
	{
		return propMetaMap;
	}


}
