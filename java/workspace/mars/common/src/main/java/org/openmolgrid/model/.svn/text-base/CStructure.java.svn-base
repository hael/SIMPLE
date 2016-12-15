/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import org.openmolgrid.format.cml.CMLReader;
import org.openmolgrid.format.cml.CMLWriter;
import org.openmolgrid.format.mopac.MopacParser;
import org.openmolgrid.format.pdb.PdbParser;
import org.openmolgrid.format.sdf.SDFParser;
import org.openmolgrid.format.sdf.SDFWriter;
import org.openmolgrid.qsar.CAttribute;
import org.openmolgrid.qsar.CDescriptorList;
import org.openmolgrid.qsar.CPropertyList;
import org.openmolgrid.qsar.StructureConsumer;
import org.xml.sax.InputSource;

/**
 * This class represents chemical structures.
 * 
 * @author Andre Lomaka
 * @author Sulev Sild
 * @author Stefan Bozic
 */
public class CStructure implements StructureConsumer, Serializable {

	/**
	 * generated serialVersionUID
	 */
	private static final long serialVersionUID = -3834464327434328019L;

	/** Constant */
	public static final int CATEGORY_2D = 0;

	/** Constant */
	public static final int CATEGORY_3D = 1;

	/** Constant */
	public static final int CATEGORY_3DOPT = 2;

	/** Constant */
	private static final String[] categoryList = new String[] { "2D", "3D", "3DOpt" };

	/** The atoms of this structure */
	private ArrayList<CAtom> atoms = new ArrayList<CAtom>();

	/** The bonds of this structure */
	private ArrayList<CBond> bonds = new ArrayList<CBond>();

	/** Additional attributes for this structure */
	private LinkedHashMap<Class<? extends CAttribute>, CAttribute> attributes = new LinkedHashMap<Class<? extends CAttribute>, CAttribute>();

	/** An unique id for this structure, if used in a system that contains several structures */
	private String id = "";

	/** The category of the structure */
	private String category = null;

	/** the name of the structure */
	private String name = "";

	/** International Chemical Identifier e.g. for ethanol(CH3CH2OH) InChI=1/C2H6O/c1-2-3/h3H,2H2,1H3 */
	private String inchi = null;

	/**
	 * Creates a new instance of CStructure
	 */
	public CStructure() {
	}

	/**
	 * Creates a copy of a given CStructure
	 * 
	 * @param s the structure to copy
	 */
	public CStructure(CStructure s) {
		makeDeepCopy(s);
	}

	/**
	 * Makes a deep copy of a given structure
	 * 
	 * @param s the structure to copy
	 */
	private void makeDeepCopy(CStructure s) {
		id = s.id;
		category = s.category;
		name = s.name;
		inchi = s.inchi;

		copyStructure(s);

		removeAttributes();
		for (CAttribute at : s.attributes.values()) {
			attributes.put(at.getClass(), at.copy());
		}
	}

	/**
	 * Sets the name for this structure
	 * 
	 * @param name the name to set
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * Returns the name for this structure
	 * 
	 * @return the name for this structure
	 */
	public String getName() {
		return name;
	}

	/**
	 * @return the atoms
	 */
	public ArrayList<CAtom> getAtoms() {
		return atoms;
	}

	/**
	 * @return the bonds
	 */
	public ArrayList<CBond> getBonds() {
		return bonds;
	}

	/**
	 * Returns the number of atoms for this structure
	 * 
	 * @return the number of atoms for this structure
	 */
	public int getAtomCount() {
		return atoms.size();
	}

	/**
	 * Returns the atom of a certain position
	 * 
	 * @param i the atom index
	 * 
	 * @return the atom of a certain position
	 */
	public CAtom getAtom(int i) {
		return (CAtom) (atoms.get(i - 1));
	}

	/**
	 * Returns the number of bonds for this structure
	 * 
	 * @return the number of bonds for this structure
	 */
	public int getBondCount() {
		return bonds.size();
	}

	/**
	 * Returns the bond of a certain position
	 * 
	 * @param i the bond index
	 * 
	 * @return the bond of a certain position
	 */
	public CBond getBond(int i) {
		return (CBond) (bonds.get(i - 1));
	}

	/**
	 * Sets a id for the structure
	 * 
	 * @param id the id to set
	 */
	public void setId(String id) {
		if (id != null)
			this.id = id;
		else
			this.id = "m1";
	}

	/**
	 * Returns the id for the strcuture
	 * 
	 * @return the id
	 */
	public String getId() {
		return id;
	}

	/**
	 * Sets the chemical identifier
	 * 
	 * @param inchi the identifier to set
	 */
	public void setInChi(String inchi) {
		this.inchi = inchi;
	}

	/**
	 * Returns InChi code.
	 * 
	 * @return String with InChi code (or null if not available)
	 */
	public String getInChi() {
		return inchi;
	}

	/**
	 * sets the category for this structure
	 * 
	 * @param categoryId the category to set
	 */
	public void setCategory(int categoryId) // TODO: use enum instead of int
	{
		if (categoryId < 0 || categoryId > 2)
			category = null;
		category = categoryList[categoryId];
	}

	/**
	 * Adds an atom to the structure
	 * 
	 * @param atom the atom to add
	 */
	public void addAtom(CAtom atom) {
		atoms.add(atom);
	}

	/**
	 * Adds an bond to the structure
	 * 
	 * @param bond the bond to add
	 */
	public void addBond(CBond bond) {
		bonds.add(bond);
	}

	/**
	 * Remove all entries from the attribute list.
	 */
	public void removeAttributes() {
		attributes.clear();
	}

	/**
	 * Returns a List of atoms of a given element 
	 * @param element the element to search for
	 * 
	 * @return a List of atoms of a given element
	 */
	public List<CAtom> getAtomsOfElement(CConstants element) {
		List<CAtom> retVal = new ArrayList<CAtom>();
		
		for (CAtom atom : atoms) {
			if (atom.isElement(element)) {
				retVal.add(atom);
			}
		}
		
		return retVal;
	}
	
	
	/**
	 * Adds an entry to the attribute list.
	 * 
	 * @param <C> the class of the attribute
	 * 
	 * @param at the attribute to set
	 */
	public <C extends CAttribute> void addAttribute(C at) {
		attributes.put(at.getClass(), at);
	}

	/**
	 * Returns the attribute of a given class
	 * 
	 * @param <C> the class of the attribute
	 * @param cls the class instance
	 * @return the attribute of a given class
	 */
	public <C extends CAttribute> C getAttribute(Class<C> cls) {
		return cls.cast(attributes.get(cls));
	}

	/**
	 * Returns the attribute of a given class
	 * 
	 * @param <C> the class of the attribute
	 * @param cls the class instance
	 * @param create if <code>true</code> creates a new attribute instance if it noct exists
	 * @return the attribute of a given class
	 */
	public <C extends CAttribute> C getAttribute(Class<C> cls, boolean create) {
		C at = getAttribute(cls);
		if (at == null & create) {
			try {
				at = cls.newInstance();
			} catch (Exception e) {
				throw new RuntimeException("Can't instantiate " + cls.getName());
			}
		}
		return cls.cast(at);
	}

	/**
	 * Returns a list with all descriptor attributes for this structure
	 * 
	 * @return a list with all descriptor attributes
	 */
	public CDescriptorList getDescriptorList() {
		return getAttribute(CDescriptorList.class);
	}

	/**
	 * Returns a list with all descriptor attributes for this structure. Creates a new descriptor list if not existing
	 * 
	 * @return a list with all descriptor attributes
	 */
	public CDescriptorList getOrCreateDescriptorList() {
		return getAttribute(CDescriptorList.class, true);
	}

	/**
	 * Returns a list with all property attributes
	 * 
	 * @return a list with all property attributes
	 */
	public CPropertyList getPropertyList() {
		return getAttribute(CPropertyList.class);
	}

	/**
	 * Returns a list with all property attributes for this structure. Creates a new property list if not existing
	 * 
	 * @return a list with all property attributes
	 */
	public CPropertyList getOrCreatePropertyList() {
		return getAttribute(CPropertyList.class, true);
	}

	/**
	 * Copy chemical structure (atom and bonding data) from another structure.
	 * 
	 * @param str - source structure
	 */
	public void copyStructure(CStructure str) {
		atoms.clear();
		for (int i = 0; i < str.atoms.size(); i++) {
			atoms.add(new CAtom(str.atoms.get(i)));
		}

		bonds.clear();
		for (int i = 0; i < str.bonds.size(); i++) {
			CBond b = str.bonds.get(i);
			CAtom a1 = atoms.get(b.getAtom(1).getId() - 1);
			CAtom a2 = atoms.get(b.getAtom(2).getId() - 1);
			bonds.add(new CBond(a1, a2, b.getType(), b.getStereo()));
		}
	}

	/**
	 * Copy atomic coordinates from another structure. Both structures must have identical order of atoms and exactly
	 * the same connectivity, because this method overwrites existing XYZ coordinates.
	 * 
	 * @param str - source structure for atomic coordinates
	 */
	public void copyCoordinatesFrom(CStructure str) {
		if (getAtomCount() != str.getAtomCount()) {
			throw new RuntimeException("number of atoms doesn't match");
		}
		for (int i = 1; i <= getAtomCount(); ++i) {
			CAtom a1 = getAtom(i);
			CAtom a2 = str.getAtom(i);
			if (!a1.getElement().equals(a2.getElement())) {
				throw new RuntimeException("Copy coordinates from different element");
			}
			atoms.set(i - 1, new CAtom(a1.getElement(), a2.getX(), a2.getY(), a2.getZ(), a1.getId()));
		}
	}

	/**
	 * Tests the presence of 3D coordinates.
	 * 
	 * @return true if visualised structure has 3D coordinates
	 */
	public boolean has3DCoordinates() {
		for (int i = 1; i <= getAtomCount(); ++i) {
			double z = getAtom(i).getZ();
			if (Math.abs(z) > 0.00001)
				return true;
		}
		return false;
	}

	/**
	 * @deprecated move this method to another class
	 * 
	 * Return a mol file representation of this structure
	 * 
	 * @return a mol file representation of this structure
	 */
	public String molFormula() {
		HashMap<String, Integer> counts = new HashMap<String, Integer>();
		for (CAtom a : atoms) {
			String element = a.getElement();
			Integer n = counts.get(element);
			if (n != null) {
				counts.put(element, n + 1);
			} else {
				counts.put(element, 1);
			}
		}

		StringBuilder sb = new StringBuilder();
		Integer nC = counts.remove("C");
		if (nC != null) {
			sb.append(nC > 1 ? "C" + nC : "C");
		}
		Integer nH = counts.remove("H");
		if (nH != null) {
			sb.append(nH > 1 ? "H" + nH : "H");
		}

		ArrayList<String> elements = new ArrayList<String>(counts.keySet());
		Collections.sort(elements);
		for (String e : elements) {
			Integer n = counts.get(e);
			sb.append(n > 1 ? e + n : e);
		}

		return sb.toString();
	}

	/**
	 * @deprecated move this method to another class
	 * 
	 * Writes the structure for the slf
	 * 
	 * @param pw the writer instance
	 */
	public void writeCDR(PrintWriter pw) {
		pw.println("  <structure id=\"" + id + "\">");
		if (getInChi() != null && !getInChi().equals("")) {
			pw.println("    <inchi>" + getInChi() + "</inchi>");
		}
		if (!getName().equals("")) {
			pw.println("    <name>" + getName() + "</name>");
		}
		writeStructureDataCDR(pw);
		pw.println("  </structure>");
	}

	/**
	 * @deprecated move this method to another class
	 * Writes a cml representation of this structure for the slf
	 * 
	 * @param pw the writer instance
	 */
	protected void writeStructureDataCDR(PrintWriter pw) {
		if (atoms.size() > 0) {
			pw.print("    <coordinates format=\"chemical/x-cml\"");
			if (category != null)
				pw.print(" category=\"" + category + "\"");
			pw.println("><![CDATA[" + CMLWriter.getCmlMolecule(this) + "]]></coordinates>");
		}

		for (CAttribute at : attributes.values()) {
			at.writeCDR(pw);
		}
	}

	/**
	 * @deprecated move this method to another class
	 * 
	 * Returns the structure as a CML molecule tag.
	 * 
	 * @return the structure as a CML molecule tag.
	 */
	public String getCStructureAsCmlMolecule() {
		return CMLWriter.getCmlMolecule(this);
	}

	/**
	 * @deprecated move this method to another class
	 * 
	 * Writes a full cml document representation of this structure as string
	 * 
	 * @return a full cml document representation of this structure as string
	 */
	public String writeCML2String() {
		return CMLWriter.getCmlDocument(this);
	}



	/**
	 * @deprecated move this method to another class
	 * 
	 * Writes the structure as a mol file
	 * 
	 * @param fileName the path of the mol file
	 */
	public void writeMolFile(String fileName) {
		PrintWriter outWriter = null;
		File f = new File(fileName);
		try {
			outWriter = new PrintWriter(new FileWriter(f));
		} catch (IOException e) {
			throw new RuntimeException(e.getMessage());
		}
		writeMolFile(outWriter);
		outWriter.close();
	}

	/**
	 * @deprecated move this method to another class
	 * 
	 * Writes the structure as a mol file to the {@link PrintWriter}
	 * 
	 * @param pw the instance for writing the data
	 */
	public void writeMolFile(PrintWriter pw) {
		SDFWriter sdf = new SDFWriter(pw);
		sdf.addStructure(this);
	}

	/**
	 * @deprecated move this method to another class
	 * 
	 * Loads mopac file data that is provided by a {@link BufferedReader}
	 * 
	 * @param in the reader that provides the mopac data
	 * 
	 * @return <code>true</code> if the mopac file has been read successfully.
	 */
	public boolean loadMopOutFile(BufferedReader in) {
		CStructure s = new MopacParser(in).getFinalGeometry();
		if (s == null)
			return false;
		makeDeepCopy(s);
		return true;
	}

	/**
	 * @deprecated move this method to another class
	 * 
	 * Loads mol file data that is provided by a {@link BufferedReader}
	 * 
	 * @param in the reader that provides the mol data
	 * 
	 * @return <code>true</code> if the mol file has been read successfully.
	 */
	public boolean loadMOLFile(BufferedReader in) {
		SDFParser sdf = new SDFParser(in);
		CStructure s = sdf.getNextStructure();
		sdf.close();
		if (s != null) {
			makeDeepCopy(s);
			return true;
		}
		return false;
	}

	/**
	 * @deprecated move this method to another class
	 * 
	 * Loads a CML file and stores it in this structure instance
	 * 
	 * @param source the path to the cml file
	 * @return <code>true</code> if the cml file has been read successfully.
	 * 
	 * @throws ChemLibException is thrown when an error occurs during file reading
	 */
	public boolean loadCMLFile(String source) throws ChemLibException {
		CMLReader cr = new CMLReader(new InputSource(new StringReader(source)), null, this);
		try {
			cr.provideStructures();
		} catch (ChemLibException e) {
			throw e;
		}
		return true;
	}

	/**
	 * @deprecated move this method to another class
	 * @deprecated move this method to another class
	 * Loads a PdbFile
	 * 
	 * @param in an {@link BufferedReader}
	 * 
	 * @return <code>true</code> if the Pdb file has been loaded sucessfully
	 */
	public boolean loadPdbFile(BufferedReader in) {
		CStructure structure = new PdbParser(in).getFinalGeometry();
		if (structure == null)
			return false;
		makeDeepCopy(structure);

		return true;
	}

	// StructureConsumer implementation

	/**
	 * @see StructureConsumer#startHandling()
	 */
	public void startHandling() {
	}

	/**
	 * @see StructureConsumer#endHandling()
	 */
	public void endHandling() {
	}

	/**
	 * @see StructureConsumer#consumeStructure(CStructure)
	 */
	public void consumeStructure(CStructure str) {
		if (atoms.size() == 0) {
			makeDeepCopy(str);
		}
	}
}
