/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common.table;

import javax.swing.JFrame;

import org.openmolgrid.client.common.DescriptorViewer;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.qsar.CDescriptorList;
import org.openmolgrid.qsar.CDescriptorMetaMap;

class ColumnMDCOutput extends Column {

	public ColumnMDCOutput(StructureList sl) {
		super(sl, "Descriptors");
	}

	@Override
	public CDescriptorList get(int row) {
		return sl.get(row).getAttribute(CDescriptorList.class);
	}
	
	/** 
	 * Displays viewer for molecular descriptors.
	 * 
	 * @param model - table model
	 * @param row - row number
	 * @param col - column number 
	 */
	@Override
	public void openViewer(StructureListModel model, int row, int col) {
		CDescriptorMetaMap meta = model.getStructureList().getDescriptorMetaData();
		CStructure str = model.getStructure(row);
		DescriptorViewer view = new DescriptorViewer(str, meta);
		JFrame frame = view.makeJFrame();
		frame.setTitle("Descriptor viewer (ID="+ str.getId() + ")");
		frame.setVisible(true);
	}

}

