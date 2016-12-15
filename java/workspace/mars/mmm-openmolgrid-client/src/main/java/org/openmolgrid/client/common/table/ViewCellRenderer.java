/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common.table;

import java.awt.Component;

import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;

import org.openmolgrid.model.CStructure;

/**
 * Table cell renderer that displays a view button. 
 * 
 * @author Sulev Sild
 */
public class ViewCellRenderer implements TableCellRenderer {

	private DefaultTableCellRenderer renderer = new DefaultTableCellRenderer();

	public ViewCellRenderer() {
		renderer.setHorizontalAlignment(JLabel.CENTER);
	}

	public Component getTableCellRendererComponent(JTable table, Object value,
			boolean isSelected, boolean hasFocus, int row, int column) {
		if (value instanceof CStructure && ((CStructure)value).getAtomCount() == 0) {
			value = "";
		}
		value = (value != null) ? "View" : "";
		return renderer.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
	}
		
}
