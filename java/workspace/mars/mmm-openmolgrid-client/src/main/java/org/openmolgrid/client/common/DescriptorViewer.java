/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.util.ArrayList;

import javax.swing.AbstractAction;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.KeyStroke;
import javax.swing.table.AbstractTableModel;

import org.openmolgrid.model.CStructure;
import org.openmolgrid.qsar.CDescriptor;
import org.openmolgrid.qsar.CDescriptorList;
import org.openmolgrid.qsar.CDescriptorMetaMap;


/**
 * This class implements panel for displaying descriptor values that ara
 * attached to a given chemical structure.
 * 
 * @author Sulev Sild
 */
public class DescriptorViewer extends JPanel 
{
	private CDescriptorMetaMap descMeta;
	private TableModel model;
	
	public DescriptorViewer(CStructure structure, CDescriptorMetaMap descMeta)
	{
		this.descMeta = descMeta;
		model = new TableModel(structure.getDescriptorList());
		buildComponents();
	}


	public JFrame makeJFrame()
	{
		final JFrame frame = new JFrame();
		frame.getContentPane().add(this);
		frame.setSize(600, 400);

		// kill frame when escape is pressed
		getInputMap().put(
			KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0),
			"hide DescriptorViewer");
		getActionMap().put("hide DescriptorViewer", new AbstractAction()
		{
			public void actionPerformed(ActionEvent e)
			{
				frame.dispose();
			}
		});


		return frame;
	}
	
	/**
	 * Build a GUI for the descriptor viewer.
	 */
	private void buildComponents()
	{
		setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		
		c.insets = new Insets(3, 3, 3, 3);

		// JTable for descriptors
		c.fill = GridBagConstraints.BOTH;
		c.weightx = 0.75;
		c.weighty = 0.75;
		JTable table = new JTable(model);
		JScrollPane scrollPane = new JScrollPane(table);
		add(scrollPane, c);
	}

	private class TableModel extends AbstractTableModel
	{
		private static final int COL_VALUE = 0;
		private static final int COL_TYPE = 1;
		private static final int COL_NAME = 2;
		
		private final String[] header = {"Value", "Type", "Name"};

		private ArrayList<CDescriptor> data;

		public TableModel(CDescriptorList dl)
		{
			data = new ArrayList<CDescriptor>();
			for (int i=0; i<dl.getNumberOfDescriptors(); ++i)
			{
				CDescriptor d = dl.getDescriptor(i);
				data.add(d);
			}
		}
		
		public int getColumnCount()
		{
			return header.length;
		}

		public String getColumnName(int col)
		{
			return header[col];
		}
		
		public int getRowCount()
		{
			return data.size();
		}
		
		public Object getValueAt(int row, int col)
		{
			CDescriptor desc = data.get(row);
			if (col == COL_VALUE)
				return new Double(desc.getValue());
			else if (col == COL_TYPE)
				return descMeta.getType(desc.getId());
			else if (col == COL_NAME) {
				String name = descMeta.getName(desc.getId());
				if (name == null) {
					return desc.getId();
				} else {
					return name;
				}
			}
			
			return "";
		}
		
	} // end of TableModel class
	
}
