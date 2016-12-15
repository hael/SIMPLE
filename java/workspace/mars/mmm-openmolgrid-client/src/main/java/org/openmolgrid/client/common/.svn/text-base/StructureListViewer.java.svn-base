/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;

import org.openmolgrid.client.common.table.NavigationToolbar;
import org.openmolgrid.client.common.table.PaginatedStructureListModel;
import org.openmolgrid.client.common.table.StructureList;
import org.openmolgrid.client.common.table.StructureListTable;
import org.openmolgrid.format.slf.StructureReader;
import org.openmolgrid.format.slf.StructureWriter;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.ChemLibException;
import org.openmolgrid.qsar.StructureConsumer;

/**
 * This class is intended for visualising the output from
 * 2Dto3Dconversion, DescriptorCalculation and SemiempiricalCalculation
 * tasks. The StructureListViewer can parse output files in StructureList,
 * DescriptorFile and SemiempiricalOutput formats.
 *
 * <p>During the reading temporary files are created to the location specified 
 * by "openmolgrid.tmpdir" system property, or "java.io.tmpdir" if it is not
 * defined.
 *
 * @author Sulev Sild
 * @author Andre Lomaka
 */
public class StructureListViewer extends JPanel implements StructureList.Listener {
	// container for holding structures that are visualised
	private PaginatedStructureListModel model;

	// table component for structures
	StructureListTable table;

	// input stream to be visualised
	private BufferedInputStream inputStream;

	// Swing components for saving and export buttons
	private JPanel exportPanel;
	private JButton saveAsButton;
	private JButton exportStructButton;
	private JButton exportDescButton;

	private NavigationToolbar nt;

	private StructureList strList = new StructureList();

	/**
	 * Constructs StructureListViewer.
	 */
	public StructureListViewer() {
		init();
	}

	/**
	 * Constructs StructureListViewer.
	 * 
	 * @param is - an input stream from where the structure list will be loaded
	 */
	public StructureListViewer(InputStream is) {
		inputStream = new BufferedInputStream(is);
		init();
	}

	/**
	 * Constructs StructureListViewer.
	 */
	public StructureListViewer(String fileName) {
		try {
			inputStream = new BufferedInputStream(new FileInputStream(fileName));
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
		init();
	}

	/**
	 * Sets input stream for data.
	 * 
	 * @param is - an input stream from where the structure list will be loaded
	 */
	public void setInput(InputStream is) {
		if (is instanceof BufferedInputStream)
			inputStream = (BufferedInputStream) is;
		else
			inputStream = new BufferedInputStream(is);
	}

	private void init()
	{
		model = new PaginatedStructureListModel(strList);
		strList.registerListener(this);
		
		exportPanel = new JPanel();
		saveAsButton = new JButton("Save as");
		saveAsButton.setEnabled(false);
		saveAsButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveAs();
			}
		});
		exportPanel.add(saveAsButton);

		// make table for structures
		table = new StructureListTable(model);
 		JScrollPane scrollPane = new JScrollPane(table);

		// lay out these components in this panel
		setLayout(new GridBagLayout());

		GridBagConstraints c = new GridBagConstraints();
		c.insets = new Insets(3, 3, 3, 3);
		c.gridy = 0;
		c.anchor = GridBagConstraints.WEST;
		add(exportPanel, c);

		c.gridy += 1;
		c.fill = GridBagConstraints.BOTH;
		nt = new NavigationToolbar(model);
		add(nt, c);
		
		c.gridy += 1;
		c.weightx = 0.75;
		c.weighty = 0.75;
 		add(scrollPane, c);
	}

	private void updateExportPanel() {
		// make export structures button if necessary
		if (exportStructButton == null && model.hasStructures()) {
			exportStructButton = new JButton("Export structures");
			exportStructButton.setEnabled(false);
			exportStructButton.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					new StructureListExporter(table).handleExport();
				}
			});
			exportPanel.add(exportStructButton);
		}

		// make export descriptors button if necessary
		if (exportDescButton == null && model.hasDescriptors()) {
			exportDescButton = new JButton("Export descriptors");
			exportDescButton.setEnabled(false);
			exportDescButton.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					new DescriptorMatrixExporter(table).handleExport();
				}
			});
			exportPanel.add(exportDescButton);
		}
		
		saveAsButton.setEnabled(model.getTotalSize() > 0);
		if (exportDescButton != null) {
			exportDescButton.setEnabled(true);
		}
		if (exportStructButton != null) {
			exportStructButton.setEnabled(true);
		}
		
		validate();
	}

	public void handleStructureAdded(int row, CStructure s) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				if (model.getRowCount() > 0) nt.refresh();
			}
		});
	}
	
	public void handleStructureRemoved(int row, CStructure s) {
	}
	
	public void handleStructuresCleared() {
	}
	
	public void handleStructureChanged(int row, CStructure s) {
	}

	/**
	 * Loads data from a file.
	 */
	private synchronized void loadData() {
		
		class Reader extends Thread implements StructureConsumer {
			private StructureReader sr;
			@Override
			public void run() {
				sr = new StructureReader(inputStream, this);
				try {
					sr.provideStructures();
					updateExportPanel();
				} catch (ChemLibException e) {
					throw new RuntimeException("Structure loading error", e);
				}
			}

			public void consumeStructure(CStructure s) {
				strList.add(s);
			}

			public void endHandling() {}
			public void startHandling() {
				strList.setDescriptorMetaData(sr.getDescriptorMetaData());
				strList.setPropertyMetaData(sr.getPropertyMetaData());
			}
		};
		
		new Reader().start();
	}
	
	private synchronized void saveAs() {
		JFileChooser saveChooser = new JFileChooser();
		int ret = saveChooser.showSaveDialog(this);
		if (ret == JFileChooser.APPROVE_OPTION) {
			try {
				File file = saveChooser.getSelectedFile();
				
				StructureWriter sw = new StructureWriter(file);
				if (strList.getPropertyMetaData().getNumberOfProperties() > 0) {
					sw.writePropertyMeta(strList.getPropertyMetaData());
				}
				if (strList.getDescriptorMetaData().getNumberOfDescriptors() > 0) {
					sw.writeDescriptorMeta(strList.getDescriptorMetaData());
				}
				for (int i=0; i<strList.size(); ++i) {
					sw.addStructure(strList.get(i));
				}
				sw.close();
			} catch (Throwable e) {
				JOptionPane.showMessageDialog(this, e.toString(), "Error", 
					JOptionPane.ERROR_MESSAGE);
			}
		}
	}

	/**
	 * This method executes whenever the outcome panel becomes visible.
	 */
	public synchronized void updateValues()
	{
		if (inputStream != null) loadData();
	}

	public static void main(final String[] args) {
		if (args.length != 1) {
			System.out.println("Usage: " + StructureListViewer.class.getName() + " filename");
			System.exit(1);
		}

		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				JFrame frame = new JFrame();
				StructureListViewer viewer = new StructureListViewer(args[0]);
				frame.getContentPane().add(viewer);
				frame.pack();
				frame.setVisible(true);
				frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
				viewer.updateValues();
			}
		});
	}

}
