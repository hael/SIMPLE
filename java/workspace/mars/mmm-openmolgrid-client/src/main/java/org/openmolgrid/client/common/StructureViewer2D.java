/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;

import javax.swing.JPanel;

import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.Convertor;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.geometry.GeometryToolsInternalCoordinates;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.renderer.Renderer2D;
import org.openscience.cdk.renderer.Renderer2DModel;



/**
 * This class visualises structures in the 2D representation by using renderer
 * from the CDK library. 
 *
 * @author Sulev Sild
 */
public class StructureViewer2D extends JPanel
{
	private Molecule mol;

	// Renderer2D model controls how molecule is displayed
	private Renderer2DModel r2dm;

	// renderer object paints the molecule
	private Renderer2D renderer;

	/**
	 * Constructs a StructureViever2D object.
	 *
	 * @param  str is a CStructure object to be visualised
	 */
	public StructureViewer2D()
	{
		setPreferredSize(new Dimension(200, 200));
		r2dm = new Renderer2DModel();
		renderer = new Renderer2D(r2dm);
		setBackground(r2dm.getBackColor());
	}


	/**
	 * Constructs a StructureViever2D object.
	 *
	 * @param  str is a CStructure object to be visualised
	 */
	public StructureViewer2D(CStructure str)
	{
		this();
		adaptCStructure(str);
	}

	public void setStructure(CStructure str)
	{
		adaptCStructure(str);
	}
	

	private void adaptCStructure(CStructure str)
	{
		if (str != null)
		{
			mol = new Molecule(Convertor.convert(str));

			// TODO: skip generation of 2D coordinates if they are present 
			StructureDiagramGenerator sdg = new StructureDiagramGenerator(mol);
			try
			{
				sdg.generateCoordinates();
			}
			catch (Exception e)
			{
				throw new RuntimeException(e.getMessage());
			}

			GeometryToolsInternalCoordinates.translateAllPositive(mol);
		} else
		{
			mol = null;
		}
		repaint();
	}
	

	/**
	 * Paints structure.
	 *
	 * @param  g is a graphics context
	 */
	public void paint(Graphics g)
	{
		super.paint(g);
		if (mol == null) return;

		r2dm.setBackgroundDimension(getSize());
		GeometryToolsInternalCoordinates.scaleMolecule(mol, getSize(), 0.8);
		GeometryToolsInternalCoordinates.center(mol, getSize());
		renderer.paintMolecule(mol, (Graphics2D)g);
	}

	/**
	 * Sets background colour. This method is overriden in order to set
	 * background colour for the structure rendering model.
	 *
	 * @param bg is a new background colour
	 */
	public void setBackground(Color bg)
	{
		super.setBackground(bg);
		if (r2dm != null)
			r2dm.setBackColor(bg);
	}

}
