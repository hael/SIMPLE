/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common.table;

import java.awt.Component;

import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;

public abstract class Column {

	private String name;
	protected StructureList sl;

	public Column(StructureList sl, String name) {
		this.name = name;
		this.sl = sl;
	}

	public abstract Object get(int row);

	public String getName() {
		return name;
	}
	
	public static TableCellRenderer getRenderer() {
		final DefaultTableCellRenderer renderer = new DefaultTableCellRenderer();	
		renderer.setHorizontalAlignment(JLabel.CENTER);
		
		return new TableCellRenderer() {
			
			public Component getTableCellRendererComponent(JTable table, Object value,
					boolean isSelected, boolean hasFocus, int row, int column) {
				value = (value != null) ? "View" : "";
				return renderer.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
			}
		};
	}
	
	public void openViewer(StructureListModel model, int row, int col) {
	}

}
