/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common;

import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.table.DefaultTableCellRenderer;

import org.openmolgrid.client.common.table.StructureList;
import org.openmolgrid.client.common.table.StructureListModel;
import org.openmolgrid.client.common.table.StructureListTable;
import org.openmolgrid.format.slf.StructureWriter;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.ChemLibException;


/**
 * This class implements panel that selects molecular structures for importing
 * to SLF format.
 *
 * @author Sulev Sild
 */
public class StructureImportPanel extends JPanel
{
	// property change events
	public static final String FILENAME = "filename";
	public static final String MODIFIED = "modified";
	
	private volatile boolean isModified = false;
	private volatile boolean readOnly   = false;
	private volatile boolean isReading  = false;
	
	private StructureList strList = new StructureList();
	private StructureListModel model;
	private StructureListTable structureTable;
	private JFileChooser fileChooser;
	private JButton removeButton;
	private JButton saveButton;
	private JButton clearButton;
	private JButton addButton;
	private JLabel strCounterLabel;

	/**
	 * Creates a GUI for the new UploadPanel instance.
	 */
	public StructureImportPanel()
	{
		setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.insets = new Insets(3, 3, 3, 3);
		
		// table
		c.gridy = 1;
		c.gridheight = 4;
		c.fill = GridBagConstraints.BOTH;
		c.weightx = 0.75;
		c.weighty = 0.75;

		model = new StructureListModel(strList);
		
		structureTable = new StructureListTable(model);			
		structureTable.setRowSelectionAllowed(true);
		
		DefaultTableCellRenderer rend = new DefaultTableCellRenderer() {
			@Override
			public String getText() {
				return "view";
			}
		};
		rend.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
		structureTable.setDefaultRenderer(CStructure.class, rend);
		
		add(new JScrollPane(structureTable), c);
		
		// start column for save, add, etc. buttons
		c.gridy = 1;
		c.gridheight = 1;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.NORTH;
		c.weightx = 0;
		c.weighty = 0;
		
		// make save button
		clearButton = new JButton("Clear");
		clearButton.setEnabled(false);
		clearButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				clearStructures();
				updateButtons();
				firePropertyChange(FILENAME, null, null);
			}
		});
		add(clearButton, c);
		
		c.gridy += 1;
		saveButton = new JButton("Save");
		saveButton.setEnabled(false);
		saveButton.addActionListener(new ActionListener()
		{
			public void actionPerformed(ActionEvent arg0) {
				handleSave();
			}	
		});
		add(saveButton, c);
		
		// make add button
		c.gridy += 1;
		addButton = new JButton("Add");
		addButton.addActionListener(new ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				handleAddFile();
			}
		});
		add(addButton, c);
		
		// make remove button
		c.gridy += 1;
		removeButton = new JButton("Remove");
		removeButton.setEnabled(false);
		removeButton.addActionListener(new ActionListener()
		{
			public void actionPerformed(ActionEvent e)
			{
				handleRemove(structureTable.getSelectedRows());
			}
		});
		add(removeButton, c);

		// adds a number of structures
 		JPanel counterPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 0,0));
 		counterPanel.add(new JLabel("Number of structures: "));
 		strCounterLabel = new JLabel("0     ");
 		counterPanel.add(strCounterLabel);
 		c.gridy += 1;
 		c.gridwidth = 2;
 		add(counterPanel, c);
		
		// add listener for table model that enables/disables the remove button
		// and updates label with number of structures
 		model.addTableModelListener(new TableModelListener()
 		{
 			public void tableChanged(TableModelEvent e)
 			{
 				isModified = true;
 				updateButtons();
				firePropertyChange(MODIFIED, null, true);
 			}
 		});
 		
 		updateButtons();
	}
	
	private void updateButtons() {
		int nStruct = countStructures();
		strCounterLabel.setText(String.valueOf(nStruct));
		
		if (readOnly || isReading) {
			clearButton.setEnabled(false);
			saveButton.setEnabled(false);
			removeButton.setEnabled(false);
			addButton.setEnabled(false);
		} else {
			clearButton.setEnabled(nStruct > 0);
			saveButton.setEnabled(isModified);
			removeButton.setEnabled(nStruct > 0);
			addButton.setEnabled(true);
		}
	}
	
	/**
	 * Returns the number of structures in the panel.
	 *
	 * @return the number of structures
	 */
	public int countStructures() {
		return model.getRowCount();
	}
	
	@Override
	public void setEnabled(boolean enabled) {
		readOnly = !enabled;
		updateButtons();
	}

	/**
	 * Saves structures to a given SLF file.
	 *
	 * <p>Fires FILENAME event with the new file value.
	 *
	 * @param file - destination file
	 * @throws IOException - on I/O error
	 */
	public void saveStructures(File file) throws IOException {
		StructureWriter sw = new StructureWriter(file);
		if (strList.getPropertyMetaData().getNumberOfProperties() > 0) {
			sw.writePropertyMeta(strList.getPropertyMetaData());
		}
		if (strList.getDescriptorMetaData().getNumberOfDescriptors() > 0) {
			sw.writeDescriptorMeta(strList.getDescriptorMetaData());
		}
		for (int i=0; i<strList.size(); i++) {
			sw.addStructure(strList.get(i));
		}
		sw.close();
		
		isModified = false;
		
		updateButtons();
		firePropertyChange(FILENAME, null, file.getAbsoluteFile());
	}

	/**
	 * Loads structures from a given SLF file.
	 *
	 * @param file - source file
	 * @throws IOException - on I/O error
	 */
	public void loadStructures(File file) throws IOException {

		try {
			if (file.getPath().toLowerCase().endsWith("slf")) {
				clearStructures();
				addFile(file);
				isModified = false;
			}
		} catch (ChemLibException e) {
			throw new IOException("Could not load file: "+file, e);
		}
		
		firePropertyChange(FILENAME, null, file.getAbsoluteFile());
	}
	
	/**
	 * Resets structure data in the panel.
	 */
	public void clearStructures() {
		strList.clear();
		model.resetColumns();
		isModified = false;
		firePropertyChange(FILENAME, null, null);
	}

	private void handleSave()
	{
		JFileChooser saveChooser = new JFileChooser();
		int ret = saveChooser.showSaveDialog(this);
		if (ret == JFileChooser.APPROVE_OPTION) {
			try {
				saveStructures(saveChooser.getSelectedFile());
			} catch (IOException e) {
				JOptionPane.showMessageDialog(this, e.toString(), "Error", 
					JOptionPane.ERROR_MESSAGE);
			}
		}
	}
	
	/**
	 * Allows user to select file and to add it to imports table.
	 */
	private void handleAddFile() {
		if (fileChooser == null){
			fileChooser = new JFileChooser();
			fileChooser.setMultiSelectionEnabled(true);
			fileChooser.addChoosableFileFilter(
				new DialogExtensionFilter("slf", "Structure List File (*.slf)"));
			fileChooser.addChoosableFileFilter(
				new DialogExtensionFilter("mol", "MDL's Mol file (*.mol)"));
			fileChooser.addChoosableFileFilter(
				new DialogExtensionFilter("mdl", "MDL's Mol file (*.mdl)"));
			fileChooser.addChoosableFileFilter(
				new DialogExtensionFilter("sdf", "MDL's Structure-Data File (*.sdf)"));
			fileChooser.addChoosableFileFilter(
				new DialogExtensionFilter("cml", "Chemical Markup Language (*.cml)"));
			fileChooser.addChoosableFileFilter(
					new DialogExtensionFilter("tsv", "Tab Separated Values (*.tsv)"));
			fileChooser.addChoosableFileFilter(
					new DialogExtensionFilter("log", "GAMESS-US output file (*.log)"));
			fileChooser.setAcceptAllFileFilterUsed(true);
		}

		int ret = fileChooser.showDialog(this, null);
		if (ret == JFileChooser.APPROVE_OPTION) {
			new Reader().start();
		}
	}

	private class Reader extends Thread {
		
		@Override
		public void run() {
			File selection[] = fileChooser.getSelectedFiles();
			final StringBuilder errors = new StringBuilder();
			isReading = true;
			try {
				for (File f : selection) {
					try {
						addFile(f);
					} catch (ChemLibException e) {
						errors.append(f.getName()).append(": ").append(e.toString());
						errors.append("\n");
					}
				}
			} finally {
				isReading = false;
			}

			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					updateButtons();
					if (errors.length() > 0) {
						JOptionPane.showMessageDialog(StructureImportPanel.this,
							errors.toString(), "Error", JOptionPane.ERROR_MESSAGE);
					}
				}
			});
		}
	}

	private void handleRemove(int[] rows) {
		for (int i=rows.length-1; i >= 0; --i) {
			strList.remove(rows[i]);
		}
		updateButtons();
	}
	
	private void addFile(File file) throws ChemLibException
	{
		StructureImporter imp = new StructureImporter(file, model.getStructureList());
		imp.provideStructures();
	}
	
	public static void main(String[] args) {
		StructureImportPanel panel = new StructureImportPanel();
		
//		panel.addPropertyChangeListener(new PropertyChangeListener() {
//			public void propertyChange(PropertyChangeEvent e) {
//				if (FILENAME.equals(e.getPropertyName())) {
//					System.out.println("Filename: "+e.getNewValue());
//				} else if (MODIFIED.equals(e.getPropertyName())) {
//					System.out.println("Modified");
//				}
//			}
//		});
		
		JFrame frame = new JFrame();
		frame.getContentPane().add(panel);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

		if (args.length > 0) {
			try {
				panel.addFile(new File(args[0]));
			} catch (ChemLibException e) {
			}
		}
	}

}
