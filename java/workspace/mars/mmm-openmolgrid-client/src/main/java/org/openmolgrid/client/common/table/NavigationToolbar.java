/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common.table;

import javax.swing.*;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Toolbar for changing browsing pages.
 * 
 * @author Andre Lomaka
 * @author Sulev Sild
 */
public class NavigationToolbar extends JPanel implements ActionListener {
	private JButton firstBtn;
	private JButton lastBtn;
	private JButton forwardBtn;
	private JButton backwardBtn;
	private JButton sharpBtn;
	private JLabel infoLabel;
	private static String START = "org/openmolgrid/client/common/start.gif";
	private static String END = "org/openmolgrid/client/common/end.gif";
	private static String FF = "org/openmolgrid/client/common/ff.gif";
	private static String BACK = "org/openmolgrid/client/common/back.gif";
	private static String SHARP = "org/openmolgrid/client/common/sharp.gif";

	private PaginatedStructureListModel bn;

	public NavigationToolbar(PaginatedStructureListModel bn) {
		setLayout(new GridBagLayout());
		GridBagConstraints gc = new GridBagConstraints();
		firstBtn = makeButton(START, "Go to the first page");
		lastBtn = makeButton(END, "Go to the last page");
		forwardBtn = makeButton(FF, "next page");
		sharpBtn = makeButton(SHARP, "Go to record #");
		sharpBtn.setEnabled(true);
		backwardBtn = makeButton(BACK, "previous page");
		infoLabel = new JLabel("                ");
		gc.anchor = GridBagConstraints.CENTER;
		gc.gridwidth = 1;
		gc.insets = new java.awt.Insets(5, 1, 5, 1);
		add(firstBtn, gc);
		add(backwardBtn, gc);
		add(sharpBtn, gc);
		add(forwardBtn, gc);
		add(lastBtn, gc);
		gc.anchor = GridBagConstraints.WEST;
		gc.weightx = 0.75;
		gc.insets = new java.awt.Insets(1, 10, 1, 1);
		add(infoLabel, gc);
		this.bn = bn;
	}

	private JButton makeButton(String iconName, String tooltip) {
		JButton button = new JButton(loadIcon(iconName));
		button.setToolTipText(tooltip);
		button.addActionListener(this);
		button.setEnabled(false);
		return button;
	}

	public void actionPerformed(ActionEvent e) {
		Object src = e.getSource();
		if (firstBtn == src) {
			bn.gotoFirstPage();
		} else if (lastBtn == src) {
			bn.gotoLastPage();
		} else if (forwardBtn == src) {
			bn.gotoNextPage();
		} else if (sharpBtn == src) {
			gotoRecord();
		} else if (backwardBtn == src) {
			bn.gotoPreviousPage();
		}

		load();
	}

	private ImageIcon loadIcon(String fileName) {
		if (fileName.trim().equals(""))
			return null;
		ClassLoader loader = getClass().getClassLoader();
		java.net.URL url;
		ImageIcon icon = null;
		try {
			url = loader.getResource(fileName);
			icon = new ImageIcon(Toolkit.getDefaultToolkit().createImage(url));
		} catch (Exception e) {
			e.printStackTrace();
		}
		return icon;
	}

	private void gotoRecord() {
		int n;
		while (true) {
			String ns = JOptionPane.showInputDialog(this, "Record Number");
			if (ns != null) {
				try {
					n = Integer.parseInt(ns);
					int ts = bn.getTotalSize();
					if (n < 1 || n > ts) {
						JOptionPane.showMessageDialog(this,
								"That record is out of range. There are " + ts
										+ " total record(s).",
								"OpenMolGRID message",
								JOptionPane.PLAIN_MESSAGE);
					} else {
						bn.gotoRecord(n);
						load();
						break;
					}
				} catch (NumberFormatException e) {
					JOptionPane.showMessageDialog(this,
							"Invalid record number", "OpenMolGRID message",
							JOptionPane.PLAIN_MESSAGE);
				}
			} else
				break;
		}
	}

	private void load() {
		int s = bn.getPageIndex();
		int e = s + bn.getRowCount() - 1;
		int r = bn.getTotalSize();
		infoLabel.setText("Showing rows " + s + "-" + e + " of " + r);
		setButtonState();
	}

	private void setButtonState() {
		firstBtn.setEnabled(bn.hasPreviousPage());
		backwardBtn.setEnabled(bn.hasPreviousPage());
		lastBtn.setEnabled(bn.hasNextPage());
		forwardBtn.setEnabled(bn.hasNextPage());
	}

	public void refresh() {
		load();
	}

}
