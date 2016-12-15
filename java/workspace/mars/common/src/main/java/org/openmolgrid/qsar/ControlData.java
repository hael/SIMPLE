/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.io.*;
import java.util.*;
import javax.xml.parsers.*;
import javax.xml.transform.*;
import javax.xml.transform.stream.*;
import javax.xml.transform.dom.*; 

import org.openmolgrid.model.ChemLibException;
import org.w3c.dom.*;
import org.xml.sax.*;

/**
 * This class is responsible for the management of the control data for the 
 * QSPR/QSAR model building tools. The control data is entered by user with 
 * the help of a UNICORE plugin. This data is saved in the XML format and 
 * later used by the UNICORE wrappers to construct a command line and/or to 
 * prepare input file(s) for the application to be executed at the target 
 * site.<p>
 * 
 * @author Sulev Sild
 */
public class ControlData implements Serializable
{
	// hash map with option and value pairs
	private Hashtable<String,String> options;

	// application name
	private String appName;

	/**
	 * Constructor to create an emty ControlData object.
	 */
	public ControlData()
	{
		appName = null;
		options = new Hashtable<String,String>();
	}

	/**
	 * Constructs an empty ControlData object for a given application. Useful
	 * if you want to set control data for the given UNICORE application and
	 * save it in the XML format.
	 * 
	 * @param appName is the name of the application
	 */
	public ControlData(String appName)
	{
		this();
		setApplication(appName);
	}

	/**
	 * Sets a name for the application.
	 *  
	 * @param name application name
	 */
	public void setApplication(String name)
	{
		appName = name;
	}

	/**
	 * Gets the name of the application.
	 *  
	 * @return application name
	 */
	public String getApplication()
	{
		return appName;
	}

	/**
	 * Sets the given option value (string).
	 * 
	 * @param name is the option name
	 * @param value is the option value
	 */
	public void setOption(String name, String value)
	{
		options.put(name, value);
	}

	/**
	 * Sets the given boolean option value. The boolean value is stored as a string with 
	 * the value of "true" or "false".
	 * 
	 * @param name is the option name
	 * @param value is the option value
	 */
	public void setOption(String name, boolean value)
	{
		options.put(name, String.valueOf(value));
	}

	/**
	 * Returns the option value with a given name.
	 * 
	 * @param name is the option name
	 * @return String with the option value
	 */
	public String getOption(String name)
	{
		return options.get(name);
	}

	/**
	 * Checks whether the named option is defined.
	 * 
	 * @param name is the option name
	 * @return true if option is defined
	 */
	public boolean hasOption(String name)
	{
		return options.containsKey(name);
	}

	/**
	 * Returns the boolean option value with a given name.
	 * 
	 * @param name is the option name
	 * @return boolean value
	 */
	public boolean getBooleanOption(String name)
	{
		return Boolean.parseBoolean(options.get(name));
	}

	/**
	 * Returns the double option value with a given name.
	 * 
	 * @param name is the option name
	 * @return double value
	 */
	public double getDoubleOption(String name) {
		return Double.parseDouble(options.get(name));
	}

	/**
	 * Parses control data from the XML format.
	 * 
	 * @param in a stream with the input data
	 */
	public void read(Reader in) throws ChemLibException
	{
		// reset old data if any
		options.clear();
		appName = null;

		// read and parse the XML data
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		Document doc = null;
		try
		{
			DocumentBuilder db = dbf.newDocumentBuilder();
			doc = db.parse(new InputSource(in));
		} catch (Exception e)
		{
			throw new ChemLibException("Can't parse control data: " +
				e.getMessage());
		}

		// set application name
		Element root = (Element) doc.getDocumentElement();
		appName = root.getNodeName();

		assignValues(root, "");
	}

	private void assignValues(Node node, String path) {
		NodeList nl = node.getChildNodes();
		for (int i=0; i<nl.getLength(); ++i) {
			node = nl.item(i);
			if (node.getNodeType() != Node.ELEMENT_NODE) {
				continue;
			}
			NodeList children = node.getChildNodes();
			if (children.getLength() == 1) {
				String name = path+node.getNodeName();
				String value = node.getFirstChild().getNodeValue();
				options.put(name, value);
			} else if (children.getLength() > 1) {
				assignValues(node, path+node.getNodeName()+".");
			}
		}
	}

	/**
	 * Writes control data to the writer in the XML format.
	 *  
	 * @param out is a destination represented by the Writer object
	 */
	public void write(Writer output) throws IOException
	{
		// initialize builder for DOM Document
		DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
		Document doc;
		try {
			DocumentBuilder builder = factory.newDocumentBuilder();
			doc = builder.newDocument();
		} catch (ParserConfigurationException e) {
			throw new IOException("Unable to serialize to XML format: "+e.getMessage());
		}
		Element root = doc.createElement(appName);
		doc.appendChild(root);

		// store options to DOM tree
		for (String key: options.keySet()) {
			Element name = doc.createElement(key);
			Text value = doc.createTextNode(options.get(key));

			String[] path = key.split("\\.");
			Element cur = root;
			if (path.length > 1) {
				for (String i: path) {
					// is the path element already under current node? If yes,
					// then step forward
					boolean found = false;
					NodeList nl = cur.getChildNodes();
					for (int j=0; j<nl.getLength(); ++j) {
						if (i.equals(nl.item(j).getNodeName())) {
							cur = (Element)nl.item(j);
							found = true;
							break;
						}
					}

					// else append the path element and step forward
					if (!found) {
						name = doc.createElement(i);
						cur.appendChild(name);
						cur = name;
					}
				}

				name.appendChild(value);
			} else {
				name.appendChild(value);
				cur.appendChild(name);
			}
		}

		// serialize DOM tree 
		DOMSource domSource = new DOMSource(doc);
		StreamResult streamResult = new StreamResult(output);
		TransformerFactory tf = TransformerFactory.newInstance();
		tf.setAttribute("indent-number", Integer.valueOf(2));
		try {
			Transformer serializer = tf.newTransformer();
			serializer.setOutputProperty(OutputKeys.INDENT,"yes");
			serializer.setOutputProperty(OutputKeys.ENCODING,"UTF-8");
			serializer.transform(domSource, streamResult); 
		} catch (TransformerException e) {
			throw new IOException("Unable to serialize to XML format: "+e.getMessage());
		}
	}

	/**
	 * Returns control data as XML formatted string.
	 *
	 * @return String with control data
	 */
	public String toString()
	{
		try
		{
			StringWriter sw = new StringWriter();
			write(sw);
			return sw.toString();
		} catch (IOException e)
		{
			throw new RuntimeException(e);
		}
	}

}
