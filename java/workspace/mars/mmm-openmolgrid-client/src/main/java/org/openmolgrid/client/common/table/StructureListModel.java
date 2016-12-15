/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common.table;

import java.util.ArrayList;
import java.util.HashSet;

import javax.swing.SwingUtilities;
import javax.swing.table.AbstractTableModel;

import org.openmolgrid.format.mopac.MopacOutput;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.qsar.CProperty;
import org.openmolgrid.qsar.CPropertyList;

public class StructureListModel extends AbstractTableModel implements StructureList.Listener {

	private StructureList strList;
	private ArrayList<Column> columns = new ArrayList<Column>();
	private HashSet<String> properties = new HashSet<String>();

	private boolean hasStructures;
	private boolean hasMopacOutput;
	private boolean hasDescriptors;
	private boolean hasName;
	

	public StructureListModel(StructureList strList) {
		this.strList = strList;
		strList.registerListener(this);
		resetColumns();
	}
	
	public void resetColumns() {
		columns.clear();
		columns.add(new ColumnStructureId(strList));
		hasStructures = false;
		hasMopacOutput = false;
		hasDescriptors = false;
		hasName = false;
		fireTableStructureChanged();
	}

	public int getColumnCount() {
		return columns.size();
	}
	
	@Override
	public String getColumnName(int col) {
		return columns.get(col).getName();
	}
	
	@Override
	public Class<?> getColumnClass(int col) {
		return columns.get(col).getClass();
	}

	public int getRowCount() {
		return strList.size();
	}

	public Object getValueAt(int row, int col) {
		return columns.get(col).get(row);
		
	}
	
	protected void updateColumns(CStructure s) {
		int oldSize = columns.size();
		
		if (hasStructures == false) { 
			hasStructures = s.getAtomCount() > 0;
			if (hasStructures) columns.add(new ColumnViewStructure(strList));
		}
		
		if (hasName == false) {
			hasName = s.getName().length() > 0;
			if (hasName ) columns.add(new ColumnName(strList));
		}
		
		CPropertyList pl = s.getAttribute(CPropertyList.class);
		if (pl != null) {
			for (int i=0; i<pl.getNumberOfProperties(); ++i)
			{
				CProperty p = pl.getProperty(i);
				if (!properties.contains(p.getId())) {
					String id = p.getId();
					properties.add(id);
					String name = strList.getPropertyMetaData().getName(id);
					columns.add(new ColumnProperty(strList, id, name));
				}
			}
		}
		
		if (hasMopacOutput == false) {
			hasMopacOutput = s.getAttribute(MopacOutput.class) != null;
			if (hasMopacOutput) columns.add(new ColumnMopacOutput(strList));
		}

		if (hasDescriptors == false) {
			hasDescriptors = s.getDescriptorList() != null;
			if (hasDescriptors) columns.add(new ColumnMDCOutput(strList));
		}
		
		if (oldSize != columns.size()) fireTableStructureChanged();
	}

	public void handleStructureAdded(final int row, final CStructure s) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				updateColumns(s);
				fireTableRowsInserted(row, row);
			}
		});
	}
	
	public void handleStructureRemoved(final int row, CStructure s) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				fireTableRowsDeleted(row, row);
			}
		});
	}
	
	public void handleStructuresCleared() {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				fireTableStructureChanged();
			}
		});
	}
	
	public void handleStructureChanged(final int row, final CStructure s) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				updateColumns(s);
				fireTableRowsUpdated(row, row);
			}
		});
	}
	
	public StructureList getStructureList() {
		return strList;
	}
	
	public boolean hasStructures() {
		return hasStructures;
	}
	
	public boolean hasDescriptors() {
		return hasDescriptors;
	}

	public CStructure getStructure(int row) {
		return strList.get(row);
	}

	public void openViewer(int row, int col) {
		columns.get(col).openViewer(this, row, col);
	}

}
