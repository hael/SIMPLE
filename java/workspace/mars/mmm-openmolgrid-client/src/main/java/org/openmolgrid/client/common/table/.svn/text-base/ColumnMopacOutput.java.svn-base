/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common.table;

import javax.swing.JFrame;

import org.openmolgrid.client.common.TextViewer;
import org.openmolgrid.format.mopac.MopacOutput;
import org.openmolgrid.model.CStructure;

class ColumnMopacOutput extends Column {

	public ColumnMopacOutput(StructureList sl) {
		super(sl, "MOPAC output");
	}

	public MopacOutput get(int row) {
		return sl.get(row).getAttribute(MopacOutput.class);
	}
	
	/**
	 * Display a MOPAC output file for the selected structure
	 *
	 * @param model - table model
	 * @param row - row number
	 * @param col - column number 
	 */
	@Override
	public void openViewer(StructureListModel model, int row, int col) {
		CStructure str = model.getStructure(row);
		String mopacOutput = str.getAttribute(MopacOutput.class).getContent();
		
		JFrame frame = TextViewer.makeJFrame(mopacOutput);
		frame.setTitle("MOPAC Output (ID=" + str.getId() + ")");
		frame.setVisible(true);
		frame.toFront();
	}

}
