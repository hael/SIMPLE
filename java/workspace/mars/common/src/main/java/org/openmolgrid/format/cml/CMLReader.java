/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.format.cml;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.StringTokenizer;

import org.openmolgrid.format.sdf.SDFParser;
import org.openmolgrid.format.xml.XMLGenericParser;
import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CBond;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.CSystem;
import org.openmolgrid.model.ChemLibException;
import org.openmolgrid.model.CBond.Stereo;
import org.openmolgrid.qsar.StructureConsumer;
import org.openmolgrid.qsar.StructureProvider;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;
import org.xml.sax.helpers.DefaultHandler;

/**
 * Parser for chemical markup language (CML).
 * 
 * @author Andre Lomaka
 * @author Stefan Bozic
 */
public class CMLReader extends XMLGenericParser implements StructureProvider {

	private CSystem system;
	private StructureConsumer sc;
	private CStructure currentMolecule;
	private CBond currentBond;
	private CAtom currentAtom;

	private ArrayList<String> atomIds = new ArrayList<String>();
	private int moleculeTagLevel = 0;
	private int bondTagLevel = 0;
	private boolean needToAssignStereo = false;
	private boolean needToAssignCharge = false;

	/**
	 * Creates a new instance of CMLReader.
	 * 
	 * @param inputFile the cml file to read
	 * @param schemaSource the schema to use
	 * @param sc the consumer of the structure
	 */
	public CMLReader(String inputFile, String schemaSource, StructureConsumer sc) {
		super(inputFile, schemaSource != null, schemaSource);
		this.sc = sc;
	}

	/**
	 * Creates a new instance of CMLReader.
	 * 
	 * @param isource an input source
	 * @param schemaSource the schema to use
	 * @param sc the consumer of the structure
	 */
	public CMLReader(InputSource isource, String schemaSource, StructureConsumer sc) {
		super(isource, schemaSource != null, schemaSource);
		this.sc = sc;
	}

	/**
	 * Creates a new instance of CMLReader.
	 * 
	 * @param inputFile the cml file to read
	 * @param schemaSource the schema to use
	 * @param system the system that represents the CML
	 */
	public CMLReader(String inputFile, String schemaSource, CSystem system) {
		super(inputFile, schemaSource != null, schemaSource);
		this.system = system;
	}

	/**
	 * @see StructureProvider#provideStructures()
	 */
	public void provideStructures() throws ChemLibException {
		startParsing();
	}

	// ==============================================
	// SAX DocumentHandler methods
	// ==============================================

	/**
	 * @see DefaultHandler#startDocument()
	 */
	public void startDocument() throws SAXException {
		if (null != sc) {
			sc.startHandling();
		}
	}

	/**
	 * @see DefaultHandler#endDocument()
	 */
	public void endDocument() throws SAXException {
		if (null != sc) {
			sc.endHandling();
		}
	}

	/**
	 * @see DefaultHandler#startElement(String, String, String, Attributes)
	 */
	public void startElement(String namespaceURI, String sName, String qName, Attributes attrs) throws SAXException {
		double x, y, z;
		String id, el, refs, fc, dr, convval;
		if ("molecule".equals(sName)) {
			moleculeTagLevel++;
			if (moleculeTagLevel == 1) {
				currentMolecule = new CStructure();
				atomIds.clear();
				currentMolecule.setId(attrs.getValue("id"));
			}
		} else if ("atom".equals(sName)) {
			currentAtom = null;
			id = attrs.getValue("id");
			el = attrs.getValue("elementType");
			fc = attrs.getValue("formalCharge");

			if (id != null && el != null && !atomIds.contains(id)) {
				atomIds.add(id);
				if (attrs.getValue("x3") != null && attrs.getValue("y3") != null && attrs.getValue("z3") != null) {
					x = Double.valueOf(attrs.getValue("x3")).doubleValue();
					y = Double.valueOf(attrs.getValue("y3")).doubleValue();
					z = Double.valueOf(attrs.getValue("z3")).doubleValue();
				} else if (attrs.getValue("x2") != null && attrs.getValue("y2") != null) {
					x = Double.valueOf(attrs.getValue("x2")).doubleValue();
					y = Double.valueOf(attrs.getValue("y2")).doubleValue();
					z = 0.0;
				} else {
					x = 0.0;
					y = 0.0;
					z = 0.0;
				}
				CAtom a = new CAtom(el, x, y, z, atomIds.size());
				currentAtom = a;
				if (fc != null) {
					if (fc.charAt(0) == '+')
						fc = fc.substring(1);
					a.setFormalCharge(Integer.parseInt(fc));
				}
				currentMolecule.addAtom(a);
			}
		} else if ("bond".equals(sName)) {
			bondTagLevel++;
			refs = attrs.getValue("atomRefs2");
			currentBond = null;
			if (refs != null) {
				StringTokenizer st = new StringTokenizer(refs);
				int bnda1 = atomIds.indexOf(st.nextToken());
				int bnda2 = atomIds.indexOf(st.nextToken());
				if (bnda1 == -1 || bnda2 == -1) {
					System.err.println("Unknown atom id during bond reading");
				} else {
					int mlp = 1;
					if (attrs.getValue("order") != null) {
						mlp = translateCMLBondOrder(attrs.getValue("order"));
					}
					currentBond = new CBond(currentMolecule.getAtom(bnda1 + 1), currentMolecule.getAtom(bnda2 + 1), mlp);
				}
			}
		} else if ("scalar".equals(sName)) {
			dr = attrs.getValue("dictRef");
			if (bondTagLevel == 1 && currentBond != null && "mdl:stereo".equals(dr)) {
				collectCharacterData();
				needToAssignStereo = true;
			}

			if (currentAtom != null && "cc:charge".equals(dr)) {
				collectCharacterData();
				needToAssignCharge = true;
			}
		} else if ("bondStereo".equals(sName)) {
			collectCharacterData();
			if (bondTagLevel == 1 && currentBond != null) {
				// bond stereo with MDL's convention
				String conv = attrs.getValue("convention");
				if ("MDL".equals(conv)) {
					convval = attrs.getValue("conventionValue");
					if (convval != null) {
						Stereo stereo = SDFParser.translateMdlStereo(Integer.parseInt(convval));
						currentBond.setStereo(stereo);
					}
				}

				// Bond stereo with CDK/JChemPaint convention
				String dictRef = attrs.getValue("dictRef");
				if (dictRef != null) {
					currentBond.setStereo(translateCMLBondStereo(dictRef.replace("cml:", "")));
				}
			}

		}
	}

	/**
	 * @see DefaultHandler#endElement(String, String, String)
	 */
	public void endElement(String namespaceURI, String sName, String qName) throws SAXException {
		if ("molecule".equals(sName)) {
			moleculeTagLevel--;
			if (moleculeTagLevel == 0) {
				if (sc != null) {
					sc.consumeStructure(currentMolecule);
				}
			}

			if (null != system) {
				system.addStructure(currentMolecule);
			}
		} else if ("bond".equals(sName)) {
			bondTagLevel--;
			if (currentBond != null)
				currentMolecule.addBond(currentBond);
		} else if ("scalar".equals(sName)) {
			String data = getCharacterData();
			if (data != null) {
				if (needToAssignStereo) {
					currentBond.setStereo(translateCMLBondStereo(data));
					needToAssignStereo = false;
				}

				if (needToAssignCharge) {
					currentAtom.setPartialCharge(Double.valueOf(data));
					needToAssignCharge = false;
				}
			}
		} else if ("bondStereo".equals(sName)) {
			String data = getCharacterData().trim();
			if (currentBond.getStereo().equals(Stereo.NONE) && data.length() > 0) {
				currentBond.setStereo(translateCMLBondStereo(data));
			}
		}

	}

	private Stereo translateCMLBondStereo(String s) {
		if ("W".equals(s))
			return Stereo.UP;
		else if ("H".equals(s))
			return Stereo.DOWN;
		return Stereo.NONE;
	}

	private int translateCMLBondOrder(String s) {
		if ("1".equals(s) || "S".equals(s))
			return 1;
		else if ("2".equals(s) || "D".equals(s))
			return 2;
		else if ("3".equals(s) || "T".equals(s))
			return 3;
		else if ("A".equals(s))
			return 4;
		return 1;
	}

	/**
	 * @see DefaultHandler#resolveEntity(String, String)
	 */
	public InputSource resolveEntity(String publicId, String systemId) throws SAXException {
		if (systemId != null) {
			if (systemId.lastIndexOf(".dtd") > -1 || systemId.lastIndexOf(".DTD") > -1) {
				InputStream ins = this.getClass().getClassLoader().getResourceAsStream(
						"org/openmolgrid/common/data/cml1_0_1.dtd");
				return new InputSource(new BufferedReader(new InputStreamReader(ins)));
			}
		}

		try {
			return super.resolveEntity(publicId, systemId);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	/**
	 * @see DefaultHandler#error(SAXParseException)
	 */
	public void error(SAXParseException e) throws SAXException {
		System.out.println("error: " + e.getMessage());
	}

	/**
	 * @see DefaultHandler#fatalError(SAXParseException)
	 */
	public void fatalError(SAXParseException e) throws SAXException {
		System.out.println("fatalError: " + e.getMessage());
	}

	/**
	 * @see DefaultHandler#warning(SAXParseException)
	 */
	public void warning(SAXParseException e) throws SAXException {
		System.out.println("warning: " + e.getMessage());
	}
}
