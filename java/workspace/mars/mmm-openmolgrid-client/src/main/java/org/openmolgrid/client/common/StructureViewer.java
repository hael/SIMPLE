/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.StringReader;

import javax.swing.AbstractAction;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.KeyStroke;
import javax.swing.WindowConstants;

import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.ChemLibException;

/**
 * This class is used by modelling plugins and data request plugin to visualise chemical structures. Both 2D and 3D representations of
 * molecular structures are supported. A renderer from CDK library is used for the visualisation of 2D structures (StructureViewer2D class)
 * and a renderer from Jmol is used for the visualisation of 3D structures (StructureViewer3D class).
 * 
 * @author Sulev Sild
 * @author Stefan Bozic
 */
public class StructureViewer extends JPanel {
	
	/** generated serialVersionUID */
	private static final long serialVersionUID = -4134874356163598248L;

	/** Structure to show */
	private CStructure str;

	/** file to show */
	private File file;

	/**
	 * Constructs a StructureViewer object.
	 * 
	 * @param structure is a CStructure object to be visualised
	 */
	public StructureViewer(CStructure structure) throws ChemLibException {
		this.str = structure;
		buildComponents();
	}

	/**
	 * Constructs a StructureViewer object.
	 * 
	 * @param molString is a String with a structure in MDL's MOL or CML format
	 */
	public StructureViewer(String molString) throws ChemLibException {
		str = parseMolString(molString);
		buildComponents();
	}

	/**
	 * Constructs a StructureViewer object.
	 * 
	 * @param molString is a String with a structure in MDL's MOL or CML format
	 */
	public StructureViewer(File file) throws ChemLibException {
		this.file = file;
		buildComponents();
	}

	/**
	 * Creates a {@link JFrame} for the viewer
	 * 
	 * @return a {@link JFrame} for the viewer
	 */
	public JFrame makeJFrame() {
		final JFrame f = new JFrame();
		f.getContentPane().add(this);
		f.setSize(400, 400);
		f.setTitle("Structure Viewer");
		f.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);

		// kill frame when escape is pressed
		getInputMap().put(KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), "hide StructureViewer");
		getActionMap().put("hide StructureViewer", new AbstractAction() {
			public void actionPerformed(ActionEvent e) {
				f.dispose();
			}
		});

		return f;
	}

	/**
	 * Builds the GUI
	 */
	private void buildComponents() throws ChemLibException {
		setLayout(new BorderLayout());

		JPanel viewer = null;

		if (str != null) {
			if (str.has3DCoordinates()) {
				viewer = new StructureViewer3D(str);
			} else {
				viewer = new StructureViewer2D(str);
			}
		}
		
		if (file != null) {
			viewer = new StructureViewer3D(file);
		}
		
		add(viewer, BorderLayout.CENTER);
	}

	/**
	 * Parses molecular structure that may be in MDL's MOL or CML format.
	 * 
	 * @param is a String with moleculat structure
	 * @return a CStructure object
	 * @throws ChemLibException on parse error
	 */
	private CStructure parseMolString(String molString) throws ChemLibException {
		CStructure mol = new CStructure();

		if (-1 != molString.lastIndexOf("</molecule>")) {
			mol.loadCMLFile(molString);
		} else {
			BufferedReader r = new BufferedReader(new StringReader(molString));
			mol.loadMOLFile(r);
		}

		return mol;
	}

}
