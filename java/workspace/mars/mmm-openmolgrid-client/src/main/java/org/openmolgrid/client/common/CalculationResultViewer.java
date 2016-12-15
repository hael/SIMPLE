/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common;

import java.awt.Component;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.logging.Logger;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.border.TitledBorder;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

import org.openmolgrid.model.ChemLibException;
import org.openmolgrid.qsar.CalculationResult;
import org.openmolgrid.qsar.CalculationResultConsumer;
import org.openmolgrid.qsar.CalculationResultReader;

/**
 * This class implements a status panel that gives overview about the 
 * completed task. The supported tasks include 2D to 3D conversion, semi-
 * empirical quantum calculation, and descriptor calculation. The displayed
 * information contains number of molecules processed, number of failures,
 * and detailed information each failure.  This information is read from a
 * calculation log file.
 *
 * @author Sulev Sild
 */
public class CalculationResultViewer 
	extends JPanel implements CalculationResultConsumer
{
	// logger
	private static final Logger logger =
		Logger.getLogger("omg.openmolgrid.common.client");

	// file or input stream for the calculation calculation log file
	private File dataFile;
	private BufferedInputStream inputStream;
	
	// table model for calculation results
	private TableModel model;

	// number of structures processed
	private int numStruct;
	private JLabel numStructLabel;

	// number of structures that failed
	private int numFailed;
	private JLabel numFailedLabel;

	// time when data file was loaded
	private long lastLoadingTime = 0;

	/**
	 * Constructs CalculationResultViewer panel.
	 */
	public CalculationResultViewer()
	{
		numStruct = 0;
		numFailed = 0;
		dataFile = new File("");
		buildComponents();
	}

	/**
	 * Constructs CalculationResultViewer panel.
	 */
	public CalculationResultViewer(String fileName)
	{
		numStruct = 0;
		numFailed = 0;
		dataFile = new File(fileName);
		buildComponents();
	}

	/**
	 * Constructs CalculationResultViewer panel.
	 */
	public CalculationResultViewer(InputStream is)
	{
		numStruct = 0;
		numFailed = 0;
		inputStream = new BufferedInputStream(is);
		dataFile = new File("");
		buildComponents();
	}

	/**
	 * Sets an input stream.
	 *
	 * @param is - InputStream
	 */
	public void setInput(InputStream is)
	{
		inputStream = new BufferedInputStream(is);
	}
	
	/**
	 * Builds the GUI.
	 */
	private void buildComponents()
	{
		setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();

		c.gridy = 1;
		c.insets = new Insets(3, 3, 3, 3);
		c.anchor = GridBagConstraints.WEST;
		add(new JLabel("Number of structures processed: "), c);
		numStructLabel = new JLabel(String.valueOf(numStruct));
		add(numStructLabel, c);

		c.gridy = 2;
		add(new JLabel("Number of structures failed: "), c);
		numFailedLabel = new JLabel(String.valueOf(numFailed));
		add(numFailedLabel, c);
		
		c.gridy = 3;
		c.gridwidth = 2;
		c.fill = GridBagConstraints.BOTH;
		c.weightx = 0.75;
		c.weighty = 0.75;
		c.insets = new Insets(10, 3, 3, 3);
		model = new TableModel();
		JTable table = new JTable(model);
		table.setRowSelectionAllowed(true);
		table.setColumnSelectionAllowed(true);
		TableColumn tc = table.getColumnModel().getColumn(TableModel.COL_ERROR);
		tc.setCellRenderer(new ErrorCellRenderer());

		JScrollPane scrollPane = new JScrollPane(table);
		scrollPane.setBorder(new TitledBorder("Failed structures: "));
		add(scrollPane, c);

		// when the error column is selected, display relevant error message
		table.addMouseListener(new MouseAdapter()
		{
			public void mouseClicked(MouseEvent e)
			{
				JTable tab = (JTable) e.getSource();
				int row = tab.getSelectedRow();
				int col = tab.getSelectedColumn();
				if (row != -1 && col == TableModel.COL_ERROR)
				{
					String error = (String)tab.getValueAt(row, col);
					JFrame frame = TextViewer.makeJFrame(error);
					frame.setTitle("Error message");
					frame.setVisible(true);
					frame.toFront();
				}
			}
		});
	}

	/**
	 * Processes CalculationResult object that was read by the 
	 * CalculationResultReader.
	 *
	 * @param result
	 */
	public void processResult(CalculationResult result)
	{
		++numStruct;
		if (! result.isSuccessful())
		{
			model.addResult(result);
			++numFailed;
		}
	}


	/**
	 * Loads calculation results.
	 *
	 * @param fname - the file name
	 */
	public synchronized void loadData()
	{
		if (dataFile.exists() && lastLoadingTime > dataFile.lastModified())
			return;

		// reset data in case data is reloaded from NJS
		numStruct = 0;
		numFailed = 0;
		model.resetData();

		CalculationResultReader dr;
		if (dataFile.exists()) 
			dr = new CalculationResultReader(dataFile.getAbsolutePath(), this);
		else
			dr = new CalculationResultReader(inputStream, this);
			
		try 
		{
			dr.produceData();
		} catch (ChemLibException e)
		{
			logger.severe("Can't read: " + e.getMessage());
		}
		numStructLabel.setText(String.valueOf(numStruct));
		numFailedLabel.setText(String.valueOf(numFailed));
		model.fireTableDataChanged();

		lastLoadingTime = System.currentTimeMillis();
	}

	// Implementation of com.pallas.unicore.client.panels.Applyable

	/**
	 * This method executes whenever the outcome panel becomes visible.
	 */
	public void updateValues()
	{
		if (dataFile.exists() || inputStream != null)
			loadData();
	}
                                                                                
	/**
	 * This method gets executed whenever the panel contents will be applied         * to an object.
	 */
	public void applyValues()
	{
	}
                                                                                
	/**
	 * This method executes whenever a new outcome arrives from the job.
	 */
	public void resetValues()
	{
		if (dataFile.exists() || inputStream != null)
			loadData();
	}
	


	/**
	 * This inner class represent a table model for the data that was read from
	 * the output files.
	 */
	class TableModel extends AbstractTableModel 
							 
	{
		public final static int COL_ID = 0;
		public final static int COL_STATUS = 1;
		public final static int COL_ERROR = 2;
		
		private String[] header = {"Structure Id", "Status", "Error"};
		private ArrayList data;

		public void addResult(CalculationResult result)
		{
			if (data == null)
				data = new ArrayList();
			data.add(result);
		}
		
		public void resetData()
		{
			data = null;
		}

		public int getColumnCount()
		{
			return header.length;
		}

		public int getRowCount()
		{
			if (data == null)
				return 0;
			else
				return data.size();
		}

		public String getColumnName(int col)
		{
			return header[col];
		}

		public Object getValueAt(int row, int col)
		{
			CalculationResult result = (CalculationResult)data.get(row);
			if (result == null)
				return "";

			if (col == COL_ID)
				return result.getId();
			else if (col == COL_STATUS)
				return result.getStatus();
			else if (col == COL_ERROR)
				return result.getError();
			else
				return "";
		}

		public Class getColumnClass(int c)
		{
			return getValueAt(0, c).getClass();
		}
	}

	/**
	 * This inner class defines renderer for error message cells.
	 */
	private class ErrorCellRenderer 
		extends JLabel
		implements TableCellRenderer
	{
		public ErrorCellRenderer()
		{
			super("View Error");
			setHorizontalAlignment(CENTER);
		}
		
		public Component getTableCellRendererComponent(JTable table, 
			Object value, boolean isSelected, boolean hasFocus, int rowIndex,
			int vColIndex)
		{
			return this;
		}
	}	

}
