/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;


/**
 * Text viewer panel.
 * 
 * @author Sulev Sild
 */
public class TextViewer extends JPanel {
	
	private JTextArea textArea = new JTextArea(24, 81);
	private JScrollPane scrollPane = new JScrollPane(textArea);
	private JPanel toolbar = new JPanel(new FlowLayout(FlowLayout.LEFT));

	public TextViewer(String text) {
		setLayout(new BorderLayout());

		toolbar.add(makeSaveAsButton());
		add(toolbar, BorderLayout.NORTH);

		textArea.setFont(new Font("Monospaced", Font.PLAIN, textArea.getFont().getSize()));
		textArea.setText(text);
		textArea.setEditable(false);
		textArea.setCaretPosition(0);
		add(scrollPane, BorderLayout.CENTER);
	}
	
	private JButton makeCloseButton(final JFrame frame) {
		JButton close = new JButton("Close");
		close.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				frame.dispose();
			}
		});
		
		toolbar.add(close, 0);
		return close;
	}
	
	private JButton makeSaveAsButton() {
		JButton saveAs = new JButton("Save As");	
		saveAs.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveAs();
			}
		});
		return saveAs;
	}

	private void saveAs() {
		JFileChooser fileChooser = new JFileChooser();
		int ret = fileChooser.showSaveDialog(TextViewer.this);
		if (ret != JFileChooser.APPROVE_OPTION) {
			return;
		}

		File f = fileChooser.getSelectedFile();
		if (f.exists()) {
			int overwrite = JOptionPane.showConfirmDialog(this, "Overwrite: "+f, "Confirm", JOptionPane.OK_OPTION);
			if (overwrite != JOptionPane.OK_OPTION) return;
		}
		
		String text = textArea.getText();
		FileWriter writer = null;
		try {
			writer = new FileWriter(f);
			writer.write(text);
			writer.close();
		} catch (IOException e) {
			String msg = "Message: " + e.getMessage();
			JOptionPane.showMessageDialog(this, msg, "I/O Error", JOptionPane.ERROR_MESSAGE);
		}
	}

	public static JFrame makeJFrame(String text)
	{
		final JFrame frame = new JFrame();
		
		final TextViewer viewer = new TextViewer(text);
		viewer.makeCloseButton(frame);
		
		frame.getContentPane().add(viewer);
		frame.setSize(viewer.getPreferredSize());
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		
		// Key bindings: close with escape, scrolling with arrow keys
		final String DISPOSE = "Dispose";
		InputMap imap = viewer.textArea.getInputMap();
		imap.put(KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), DISPOSE);
		imap.put(KeyStroke.getKeyStroke(KeyEvent.VK_UP, 0), "negativeUnitIncrement");
		imap.put(KeyStroke.getKeyStroke(KeyEvent.VK_DOWN, 0), "positiveUnitIncrement");
		viewer.textArea.getActionMap().put(DISPOSE, new AbstractAction()
		{
			public void actionPerformed(ActionEvent e)
			{
				frame.dispose();
			}
		});
		
		frame.addWindowFocusListener(new WindowAdapter() {
		    public void windowGainedFocus(WindowEvent e) {
		        viewer.textArea.requestFocusInWindow();
		    }
		});

		frame.pack();
		return frame;
	}

}
