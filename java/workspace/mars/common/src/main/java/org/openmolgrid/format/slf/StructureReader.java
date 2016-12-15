/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.format.slf;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringReader;
import java.lang.reflect.Constructor;
import java.util.concurrent.CancellationException;

import org.openmolgrid.format.gamess.GamessOutput;
import org.openmolgrid.format.mopac.MopacOutput;
import org.openmolgrid.format.smiles.SmilesString;
import org.openmolgrid.format.xml.XMLGenericParser;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.ChemLibException;
import org.openmolgrid.qsar.ApplicationOutput;
import org.openmolgrid.qsar.CDescriptorList;
import org.openmolgrid.qsar.CDescriptorMetaMap;
import org.openmolgrid.qsar.CPropertyList;
import org.openmolgrid.qsar.CPropertyMetaMap;
import org.openmolgrid.qsar.MetaHandler;
import org.openmolgrid.qsar.StructureConsumer;
import org.openmolgrid.qsar.StructureProvider;
import org.openmolgrid.util.CUtil;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;

/**
 * This class implements a SAX handler and parser for structure list files.
 *
 * @author Andre Lomaka
 * @author Sulev Sild
 */
public class StructureReader extends XMLGenericParser implements StructureProvider {
	private CStructure str;
	private CDescriptorList dl;
	private CPropertyList pl;
	private StructureConsumer sc;
	private String id, fmt;
	private String inchi = "";
	private String name = "";
	private int structureCounter;
	private File tempOutputDir;

	// SAX handler that collects and holds metadata
	private MetaHandler metaHandler;

	// true if mopac/gamess outcome is stored on temp file, false if in the RAM
	private boolean swapOutputData = false;

	/**
	 * Constructs a new StructureReader object to parse a given StructureList
	 * file.
	 *
	 * @param inputFile is a String object with input file name.
	 * @param sc is an abject implementing StructureConsumer interface.
	 */
	public StructureReader(String inputFile, StructureConsumer sc) {
		super(inputFile, false, null);
		this.sc = sc;
		metaHandler = new MetaHandler();
	}

	/**
	 * Constructs a new StructureReader object to parse a StructureList from
	 * input stream.
	 *
	 * @param inputStream - input stream for the structure list
	 * @param sc - an abject implementing StructureConsumer interface.
	 */
	public StructureReader(InputStream inputStream, StructureConsumer sc) {
		super(new InputSource(inputStream), false, null);
		this.sc = sc;
		metaHandler = new MetaHandler();
	}

	public void provideStructures() throws ChemLibException {
		structureCounter = 0;
		tempOutputDir = null;
		try {
			startParsing();
		} catch (CancellationException e) {
			// ignored
		}
	}

	/**
	 * Returns reference to the descriptor metadata.
	 *
	 * @return CDescriptorMetaMap
	 */
	public CDescriptorMetaMap getDescriptorMetaData() {
		return metaHandler.getDescriptorMetaData();
	}

	/**
	 * Returns reference to the property metadata.
	 *
	 * @return CPropertyMetaMap
	 */
	public CPropertyMetaMap getPropertyMetaData() {
		return metaHandler.getPropertyMetaData();
	}

	//==============================================
	// SAX DocumentHandler methods
	//==============================================

	public void startDocument() throws SAXException {
		sc.startHandling();
	}

	public void endDocument() throws SAXException {
		sc.endHandling();
	}

	public void startElement(String namespaceURI, String sName, String qName, Attributes attrs) throws SAXException {
		if (metaHandler.isActive() || "metadata".equals(sName)) {
			metaHandler.startElement(namespaceURI, sName, qName, attrs);
		} else if ("structure".equals(sName)) {
			str = new CStructure();
			id = attrs.getValue("id");
			++structureCounter;
		} else if ("name".equals(sName)) {
			collectCharacterData();
		} else if ("inchi".equals(sName)) {
			collectCharacterData();
		} else if ("coordinates".equals(sName)) {
			collectCharacterData();
			fmt = attrs.getValue("format");
		} else if ("mopacOutput".equals(sName)) {
			collectCharacterData();
		} else if ("gamessOutput".equals(sName)) {
			collectCharacterData();
		} else if ("smiles".equals(sName)) {
			collectCharacterData();
		} else if ("descriptorList".equals(sName)) {
			dl = new CDescriptorList();
		} else if ("descriptor".equals(sName)) {
			String did = attrs.getValue("id");
			String dvalue = attrs.getValue("value");
			String dpos = attrs.getValue("pos");
			if (did != null && dvalue != null) {
				if (dpos == null)
					dl.addDescriptor(did, Double.parseDouble(dvalue));
				else
					dl.addDescriptor(did, Integer.parseInt(dpos),
							Double.parseDouble(dvalue));
			}
		} else if ("propertyList".equals(sName)) {
			pl = new CPropertyList();
		} else if ("property".equals(sName)) {
			String pid = attrs.getValue("id");
			String pvalue = attrs.getValue("value");
			String pmi = attrs.getValue("modelId");
			if (pid != null && pvalue != null) {
				if (pmi != null)
					pl.addProperty(pid, Double.parseDouble(pvalue), pmi);
				else
					pl.addProperty(pid, Double.parseDouble(pvalue));
			}
		}
	}

	public void endElement(String namespaceURI, String sName, String qName)
			throws SAXException {
		String data;

		if (metaHandler.isActive()) {
			metaHandler.endElement(namespaceURI, sName, qName);
		} else if ("structure".equals(sName)) {
			if (id != null)
				str.setId(id);
			str.setName(name);
			str.setInChi(inchi);
			sc.consumeStructure(str);
		} else if ("name".equals(sName)) {
			name = getCharacterData().trim();
		} else if ("inchi".equals(sName)) {
			inchi = getCharacterData().trim();
		} else if ("coordinates".equals(sName)) {
			data = getCharacterData();
			if (data != null) {
				if (fmt == null || "chemical/x-mdl-molfile".equals(fmt)) {
					str.loadMOLFile(new BufferedReader(new StringReader(data)));
				} else if ("chemical/x-cml".equals(fmt)) {
					try {
						str.loadCMLFile(data);
					} catch (ChemLibException e) {
						throw new SAXException(e.toString(), e);
					}
				}
			}
		} else if ("mopacOutput".equals(sName)) {
			data = getCharacterData();
			if (data != null) {
				MopacOutput mopout = storeOutput(data, MopacOutput.class);
				str.addAttribute(mopout);
			}
		} else if ("gamessOutput".equals(sName)) {
			data = getCharacterData();
			if (data != null) {
				GamessOutput gamout = storeOutput(data, GamessOutput.class);
				str.addAttribute(gamout);
			}
		} else if ("smiles".equals(sName)) {
			data = getCharacterData();
			if (data != null) {
				str.addAttribute(new SmilesString(data));
			}
		} else if ("descriptorList".equals(sName)) {
			str.addAttribute(dl);
		} else if ("propertyList".equals(sName)) {
			str.addAttribute(pl);
		}
	}

   public void characters(char buf[], int offset, int len) throws SAXException
   {
      if (metaHandler.isActive())
         metaHandler.characters(buf, offset, len);
      else
         super.characters(buf, offset, len);
   }
   
   /**
    * Enables saving memory by storing ApplicationOutput data on temporary 
    * files. Useful when keeping large structure lists from MOPAC or other 
    * calculations.
    * 
    * @param enable - true if enabled
    */
   public void useSwapping(boolean enable)
   {
      swapOutputData = enable;
   }
   
   private <O extends ApplicationOutput> O storeOutput(String data, Class<O> cl) {
		String fname = null;

		if (swapOutputData) {
			try {
				if (tempOutputDir == null) {
					tempOutputDir = File.createTempFile("_slftemp_", "_dir", new File ("."));
					tempOutputDir.delete();
					tempOutputDir.mkdir();
					tempOutputDir.deleteOnExit();
				}
				File f = File.createTempFile("str" + structureCounter + "_", ".dat", tempOutputDir);
				f.deleteOnExit();
				CUtil.string2File(data, f);
				fname = f.getAbsolutePath();
				data = null;
			} catch (IOException e) {
				throw new RuntimeException("Unable to create temporary file", e);
			}
		}

		try {
			Constructor<O> ctor = cl.getConstructor(String.class, String.class);
			return ctor.newInstance(fname, data);
		} catch (RuntimeException e) {
			throw e;
		} catch (Exception e) {
			return null;
		}
	}

}
