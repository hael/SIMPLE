/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.format.xml;

import java.io.File;
import java.io.IOException;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.openmolgrid.model.ChemLibException;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXNotRecognizedException;
import org.xml.sax.SAXNotSupportedException;
import org.xml.sax.helpers.DefaultHandler;

/**
 * Helper class for creating SAX parsers.
 * 
 * @author Andre Lomaka
 * @author Sulev Sild
 */
public class XMLGenericParser extends DefaultHandler
{
   protected String inputFile=null;
   private InputSource isource=null;
   private DefaultHandler handler;
   private SAXParser saxParser;
   private boolean needToReadString=false;
   private StringBuffer textBuffer=new StringBuffer();

   private static final String JAXP_SCHEMA_LANGUAGE =
      "http://java.sun.com/xml/jaxp/properties/schemaLanguage";
   private static final String W3C_XML_SCHEMA =
      "http://www.w3.org/2001/XMLSchema";
   private static final String JAXP_SCHEMA_SOURCE =
      "http://java.sun.com/xml/jaxp/properties/schemaSource";

   public XMLGenericParser(InputSource isource, boolean doValidation, String schemaSource)
   {
      this.isource = isource;
      prepareParser(doValidation, schemaSource);
   }

   public XMLGenericParser(String inputFile, boolean doValidation, String schemaSource)
   {
      this.inputFile = inputFile;
      prepareParser(doValidation, schemaSource);
   }

   /**
    * Constructs a new XMLGenericParser object. This constructor has been
    * added to make children of the XMLGenericParser class useable without
    * input file or InputSource (i.e. SAX events are directly provided by
    * the DocumentHandler interface).
    */
   public XMLGenericParser()
   {
   }

   private void prepareParser(boolean doValidation, String schemaSource)
   {
      handler = this;
      SAXParserFactory factory = SAXParserFactory.newInstance();
      factory.setNamespaceAware(true);
      if (doValidation)
      {
         factory.setValidating(true);
      }
      try
      {
         saxParser = factory.newSAXParser();
      } catch (Throwable t)
      {
         throw new RuntimeException(t.toString());
      }
      if (doValidation)
      {
         try
         {
            saxParser.setProperty(JAXP_SCHEMA_LANGUAGE, W3C_XML_SCHEMA);
         } catch (SAXNotRecognizedException x)
         {
            System.err.println(x.toString());
            throw new RuntimeException(x.toString());
         } catch (SAXNotSupportedException x)
         {
            System.err.println(x.toString());
            throw new RuntimeException(x.toString());
         }
         try
         {
            saxParser.setProperty(JAXP_SCHEMA_SOURCE, new File(schemaSource));
         } catch (SAXNotRecognizedException x)
         {
            System.err.println(x.toString());
            throw new RuntimeException(x.toString());
         } catch (SAXNotSupportedException x)
         {
            System.err.println(x.toString());
            throw new RuntimeException(x.toString());
         }
      }
   }

	protected void startParsing() throws ChemLibException {
		try {
			if (inputFile != null)
				saxParser.parse(new File(inputFile), handler);
			else if (isource != null)
				saxParser.parse(isource, handler);
		} catch (SAXException e) {
			throw new ChemLibException("SAX error during XML parsing: "+e.getMessage(), e);
		} catch (IOException e) {
			throw new ChemLibException("I/O error during XML parsing: "+e.getMessage(), e);
		}
 	}

   public void characters(char buf[], int offset, int len) throws SAXException
   {
      if (needToReadString)
      {
         String s = new String(buf, offset, len);
         textBuffer.append(s);
      }
   }

   /**
    * Instructs parser to collect character data. You can use
    * getCharacterData method to get collected character data.
    * Has no effect when data collection already in progress.
    *
    */
   public void collectCharacterData()
   {
      if (needToReadString) return;
      needToReadString=true;
      textBuffer.setLength(0);
   }

   /**
    * Returns character data that has been collected and ends collection of data.
    *
    * @return String with character data when data collection in progress, null otherwise.
    */
   public String getCharacterData()
   {
      if (!needToReadString) return null;
      needToReadString=false;
      return textBuffer.toString();
   }
}
