/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common.table;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

import org.openmolgrid.client.common.StructureViewer;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.ChemLibException;

public class ColumnViewStructure extends Column {

	public ColumnViewStructure(StructureList sl) {
		super(sl, "Structure");
	}
	
	public CStructure get(int row) {
		return sl.get(row);
	}

	/**
	 * Displays selected CStructure in the visualisation window.
	 *
	 * @param model - table model
	 * @param row - row number
	 * @param col - column number 
	 */
	@Override
	public void openViewer(StructureListModel model, int row, int col) {
		CStructure str = model.getStructure(row);
		
		if (str.getAtomCount() == 0) return;

		try
		{
			StructureViewer view = new StructureViewer(str);
			JFrame frame = view.makeJFrame();
			frame.setTitle("Structure (ID=" + str.getId() + ")");
			frame.setVisible(true);
		} catch (ChemLibException e)
		{
			String message = "There was an error starting the viewer: "+e.getMessage();
			JOptionPane.showMessageDialog(null, message, "Error", JOptionPane.ERROR_MESSAGE);
			return;
		}
	}

}
