/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Rectangle;
import java.io.File;
import java.io.PrintWriter;
import java.io.StringWriter;

import javax.swing.JPanel;

import org.jmol.adapter.smarter.SmarterJmolAdapter;
import org.jmol.api.JmolAdapter;
import org.jmol.api.JmolViewer;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.ChemLibException;

/**
 * This class visualises structures in 3D representation by using Jmol renderer.
 * 
 * @author Sulev Sild
 * @author Stefan Bozic
 */
public class StructureViewer3D extends JPanel {

	/** generated serialVersionUID */
	private static final long serialVersionUID = 843300950789302322L;

	/** JMol viewer */
	private JmolViewer viewer;

	/** JMol adapter */
	private JmolAdapter adapter;

	/** Clipping size for the JPanel */
	final Dimension currentSize = new Dimension();
	
	/** Clipping topology for the JPanel */
	final Rectangle rectClip = new Rectangle();
	
	/**
	 * Constructs a StructureViewer3D object.
	 * 
	 * @param structure is a CStructure object to be visualised
	 * @throws ChemLibException if renderer can't handle the structure
	 */
	public StructureViewer3D(CStructure structure) throws ChemLibException {
		adapter = new SmarterJmolAdapter();
		viewer = JmolViewer.allocateViewer(this, adapter);

		// The following is quite ugly but works. The molecule is
		// converted to a string in the MOL file format and then parsed by
		// Jmol. There should be a more direct way for doing that.
		StringWriter output = new StringWriter();
		structure.writeMolFile(new PrintWriter(output));
		viewer.openStringInline(output.toString());
		String strError = viewer.getErrorMessage();
		if (strError != null) {
			throw new ChemLibException(strError);
		}
	}

	/**
	 * Constructs a StructureViewer3D object.
	 * 
	 * @param file to be visualised
	 * @throws ChemLibException if renderer can't handle the structure
	 */
	public StructureViewer3D(File file) throws ChemLibException {
		if (file == null) {
			throw new ChemLibException("File is null!");			
		} else {
			
			adapter = new SmarterJmolAdapter();
			viewer = JmolViewer.allocateViewer(this, adapter);
			viewer.openFileAsynchronously(file.getAbsolutePath());
			
			String strError = viewer.getErrorMessage();
			if (strError != null) {
				throw new ChemLibException(strError);		
			}
		}
	}
	
	/**
	 * paint method of the JPanel
	 */
	public void paint(Graphics g) {
		getSize(currentSize);
		g.getClipBounds(rectClip);
		viewer.renderScreenImage(g, currentSize, rectClip);
	}

}
