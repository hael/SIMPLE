/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;

import javax.swing.AbstractListModel;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.openmolgrid.model.CStructure;
import org.openmolgrid.qsar.CDescriptorList;
import org.openmolgrid.qsar.CPropertyList;

/**
 * Maps data fields from SDF and CSV files.  
 * 
 * @author Sulev Sild
 */
public class DataFieldMapper extends JDialog {

	private static final long serialVersionUID = 1L;
	
	private LinkedHashMap<String,Mapper> mappers = new LinkedHashMap<String,Mapper>();
	
	private AvailableMappers availableMappers = new AvailableMappers();
	private JList targetList = new JList(availableMappers);

	private JList unmappedList = new JList();

	
	public DataFieldMapper(JFrame parent) {
		super(parent, "Data mapping", true);
		buildComponents();
		setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
	}
	
	public void apply(Map<String, String> data, CStructure str) {
		final ArrayList<String> undefined = new ArrayList<String>();
		for (String name: data.keySet()) {
			if (!mappers.containsKey(name)) {
				undefined.add(name);
				mappers.put(name, availableMappers.getDefaultMapper());
			}
		}
		
		if (undefined.size() > 0) {
			Runnable showMapper = new Runnable() {
				public void run() {
					unmappedList.setListData(undefined.toArray());
					unmappedList.setSelectedIndex(0);
					pack();
					setVisible(true);			
				}
			};
			try {
				SwingUtilities.invokeAndWait(showMapper);
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
			} catch (InvocationTargetException e) {
				throw new RuntimeException(e);
			}
		}
		
		for (Entry<String, String> e: data.entrySet()) {
			mappers.get(e.getKey()).apply(str, e.getKey(), e.getValue());
		}
	}
	
	void buildComponents() {
		JPanel panel = new JPanel(new GridBagLayout());

		GridBagConstraints gbc = new GridBagConstraints();
		gbc.fill = GridBagConstraints.BOTH;
		gbc.weightx = 1.0;
		gbc.insets = new Insets(3,10,3,3);
		gbc.gridy = 0;
		
		// 1st row
		panel.add(new JLabel("Data field:         "), gbc);
		panel.add(new JLabel("Use data as:        "), gbc);
		
		// 2nd row
		gbc.gridy += 1;
		
		unmappedList.addListSelectionListener(availableMappers);

		targetList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		targetList.addListSelectionListener(availableMappers);
		
		gbc.weighty = 1.0;
		panel.add(new JScrollPane(unmappedList), gbc);
		panel.add(new JScrollPane(targetList), gbc);

		JPanel content = new JPanel(new BorderLayout());
		content.add(panel, BorderLayout.CENTER);
		
		JButton ok = new JButton("OK");
		ok.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				dispose();
			}
		});
		content.add(ok, BorderLayout.SOUTH);
		
		setContentPane(content);
		pack();
	}
	
	/**
	 * Maintains a list of mappers available to selected data field(s). 
	 */
	class AvailableMappers extends AbstractListModel implements ListSelectionListener {
		
		private static final long serialVersionUID = 1L;
		
		private Mapper[] avail = {
				new IgnoreMapper(), new PropertyMapper(), new DescriptorMapper(), 
				new IdMapper(), new NameMapper()
		}; 

		public Mapper getDefaultMapper() {
			return avail[0];
		}

		public Object getElementAt(int i) {
			return avail[i];
		}

		public int getSize() {
			return avail.length;
		}

		public void valueChanged(ListSelectionEvent ev) {
			if (ev.getValueIsAdjusting()) {
				return;
			}
			if (ev.getSource() == unmappedList) {
				int nSelected = unmappedList.getSelectedValues().length;
				if (nSelected != 1) {
					targetList.clearSelection();
					return;
				}
				String sel = (String) unmappedList.getSelectedValue();
				Mapper m = mappers.get(sel);
				targetList.setSelectedValue(m, true);
			} else if (ev.getSource() == targetList) {
				Mapper sel = (Mapper) targetList.getSelectedValue();
				if (sel != null) {
					for (Object name: unmappedList.getSelectedValues()) {
						mappers.put((String)name, sel);
					}
				}
			}
		}
	}

	abstract class Mapper {
		abstract public void apply(CStructure str, String name, String value);
		abstract public String toString();
	}

	class IgnoreMapper extends Mapper {
		@Override
		public void apply(CStructure str, String name, String value) {
		}
		@Override
		public String toString() {
			return "Ignore";
		}
	}
	
	class PropertyMapper extends Mapper {
		@Override
		public void apply(CStructure str, String propid, String value) {
			CPropertyList pl = str.getOrCreatePropertyList();
			try {
				pl.addProperty(propid, Double.parseDouble(value));
			} catch (NumberFormatException e) {
				// ignore: caused by missing values
			}
		}
		@Override
		public String toString() {
			return "Property";
		}
	}
	
	class DescriptorMapper extends Mapper {
		@Override
		public void apply(CStructure str, String descid, String value) {
			CDescriptorList dl = str.getOrCreateDescriptorList();
			try {
				dl.addDescriptor(descid, Double.parseDouble(value));
			} catch (NumberFormatException e) {
				// ignore: caused by missing value
			}
		}
		@Override
		public String toString() {
			return "Descriptor";
		}
	}
	
	class NameMapper extends Mapper {
		@Override
		public void apply(CStructure str, String name, String value) {
			str.setName(value.trim());	
		}
		@Override
		public String toString() {
			return "Chemical name";
		}
	}
	
	class IdMapper extends Mapper {
		@Override
		public void apply(CStructure str, String name, String value) {
			str.setId(value);
		}
		@Override
		public String toString() {
			return "Structure ID";
		}
	}

}
