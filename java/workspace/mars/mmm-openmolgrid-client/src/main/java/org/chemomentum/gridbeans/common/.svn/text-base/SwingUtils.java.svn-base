/*
 * (c) University of Tartu, Chemomentum project
 */

package org.chemomentum.gridbeans.common;

import java.awt.GridBagConstraints;
import javax.swing.*;

/**
 * Swing utilities for building dialogs.
 *
 * @author Sulev Sild
 */
public class SwingUtils 
{
	/**
	 * Utility function to insert JLabel and JTextField pair to the panel.
	 * 
	 * @param label a text entry field label
	 * @param c constraints for the layout manager
	 * @param p destination panel
	 * @return instance of the JTextField object
	 */
	public static JTextField addTextEntryField(
		String label,
		GridBagConstraints c,
		JPanel p)
	{
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 0.0;
		p.add(new JLabel(label), c);
		c.anchor = GridBagConstraints.EAST;
		c.weightx = 0.75;
		JTextField tf = new JTextField(30);
		p.add(tf, c);
		return tf;
	}

	/**
	 * Utility function to insert JLabel and JComboBox pair to the panel.
	 * 
	 * @param label is a name of the field
	 * @param choices string array with choices
	 * @param c constraints for the layout manager
	 * @param p destination panel
	 * @return instance of the JComboBox object
	 */
	public static JComboBox addComboBox(
		String label,
		String[] choices,
		GridBagConstraints c,
		JPanel p)
	{
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 0.0;
		p.add(new JLabel(label), c);
		c.anchor = GridBagConstraints.EAST;
		c.weightx = 0.75;
		JComboBox cb = new JComboBox(choices);
		cb.setEditable(false);
		p.add(cb, c);
		return cb;
	}

	/**
	 * Utility function to insert JLabel and JCheckBox pair to the panel.
	 * 
	 * @param label is a name of the field
	 * @param c constraints for the layout manager
	 * @param p destination panel
	 * @return instance of the JTextField object
	 */
	public static JCheckBox addCheckBox(String label, GridBagConstraints c, JPanel p)
	{
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.WEST;
		c.weightx = 0.0;
		c.gridwidth = 2;
		JCheckBox cb = new JCheckBox(label);
		p.add(cb, c);
		return cb;
	}	
}
