/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * GuiUtils.java
 *
 * Created on Mar 1, 2012, 5:24:07 PM
 */
package edu.kit.gridbeans.simpleeovolassemble.plugin;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionListener;
import java.lang.reflect.Method;
import java.sql.Wrapper;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JLabel;

/**
 * Class to be used for various utilities and function throughOut other classes.
 * 
 * {@link WelcomePanel EoRecVolPanel}
 *  
 * @author Frederic Bonnet
 * @Date 24th of November 2015
 *
 */
public class GuiUtils {

	private static final Dimension BUTTON_SIZE = new Dimension(100,25);
    private static final Dimension BUTTON_SIZE_2 = new Dimension(120,25);
    private static final Dimension BUTTON_SIZE_3 = new Dimension(200,25);

    /**
     * {@link Method} to create the ImageIcon for the images.
     * @param filename
     * @return new ImageIcon(filename):ImageIcon
     */
    public static ImageIcon makeImageIconFromJarredImage(String filename) {
//		return new ImageIcon(GuiUtils.class.getClassLoader().getResource(filename));
    	return new ImageIcon(filename);
    }
    /**
     * {@link Wrapper} that calls 
     * @param text
     * @param aa
     * @return makeJLabel(text,false,aa):JLabel
     */
	public static JLabel makeJLabel(String text,boolean aa)
	{
		return makeJLabel(text,false,aa);
	}
	/**
	 * {@link Method} that overloads the JLabel class
	 * @param text
	 * @param bold
	 * @param aa
	 * @return label:JLabel
	 */
	private static JLabel makeJLabel(String text,boolean bold,boolean aa)
	{
		JLabel label;
		
		if (aa)
		label = new JLabel(text);
		else
			label = new JLabel(text);
		
		if (bold)
			label.setFont(boldFont(label.getFont()));
		
		return label;
	}
	/**
	 * {@link Method} that creates a boldFont
	 * @param font
	 * @return font.deriveFont(Font.BOLD):Font
	 */
	public static Font boldFont(Font font)
	{
		return font.deriveFont(Font.BOLD);
	}
   /**
     * method to create the JButton of custom size
     * @param text
     * @return button
     */
	public static JButton makeJButton(String text) {
		return makeJButton(text,null);
	}
	/**
	 * method to create a button of custom size: Dimension(100,25)
	 * @param text
	 * @param listener
	 * @return button:JButton
	 */
	public static JButton makeJButton(String text,ActionListener listener) {
		JButton button;
		
		button = new JButton(text);
		button.setPreferredSize(BUTTON_SIZE);
		
		if (listener != null)
			button.addActionListener(listener);
		
		return button;
	}
	/**
	 * method to create a button of custom size: Dimension(120,25)
	 * @param text
	 * @param listener
	 * @return button:JButton
	 */
	public static JButton makeJButton_2(String text,ActionListener listener) {
		JButton button;
		
		button = new JButton(text);
		button.setPreferredSize(BUTTON_SIZE_2);
		
		if (listener != null)
			button.addActionListener(listener);
		
		return button;
	}
	/**
	 * method to create a button of custom size: Dimension(200,25)
	 * @param text
	 * @param listener
	 * @return button:JButton
	 */
	public static JButton makeJButton_3(String text,ActionListener listener) {
		JButton button;
		
		button = new JButton(text);
		button.setPreferredSize(BUTTON_SIZE_3);
		
		if (listener != null)
			button.addActionListener(listener);
		
		return button;
	}
    

}
