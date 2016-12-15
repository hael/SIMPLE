/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common.table;

import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JTable;

import org.openmolgrid.client.common.table.ColumnMDCOutput;
import org.openmolgrid.client.common.table.ColumnMopacOutput;
import org.openmolgrid.client.common.table.ColumnProperty;
import org.openmolgrid.client.common.table.ColumnViewStructure;

public class StructureListTable extends JTable {

	public StructureListTable(final StructureListModel model) {
		super(model);
		setDefaultRenderer(ColumnViewStructure.class, new ViewCellRenderer());
		setDefaultRenderer(ColumnMopacOutput.class, new ViewCellRenderer());
		setDefaultRenderer(ColumnMDCOutput.class, new ViewCellRenderer());
		
		addMouseListener(new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				int row = rowAtPoint(e.getPoint());
				int col = columnAtPoint(e.getPoint());
				if (row == -1 || col == -1) {
					return;
				}

				if (getValueAt(row, col) == null) {
					return;
				}
				
				model.openViewer(row, col);
			}
		});
	}

	public StructureList getStructureList() {
		return ((StructureListModel)getModel()).getStructureList();
	}

}
