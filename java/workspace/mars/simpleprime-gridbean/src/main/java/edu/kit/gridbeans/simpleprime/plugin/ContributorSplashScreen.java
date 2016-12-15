/**
 * 
 */
package edu.kit.gridbeans.simpleprime.plugin;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SpringLayout;

/**
 * {@link WelcomePanel} for SimpleGridbean
 *   
 * @author Frederic Bonnet
 * 
 */
public class ContributorSplashScreen extends JDialog {

	/**
	 * Add UID
	 */
    private static final long serialVersionUID = 7377459880071114807L;

    private final static int WIDTH = 800;       /**width dimension of the window*/
    private final static int HEIGHT = 480;      /**height dimension of the window*/
	
    private JButton mButtonContinue = null;
	@SuppressWarnings("unused")
	private JButton	mClicked = null;
	private ActionListener mButtonListener = new ButtonListener();
	
    public ContributorSplashScreen() {
	super();
		
	setTitle("Contributors to the GridBean (java) "+VersionPrime.getVersionPrimeGridBean());
	setModal(true);
	setSize(WIDTH, HEIGHT);

	Container contentPane = getContentPane();

	contentPane.setLayout(new BorderLayout());

	JPanel imagePanel = new JPanel(new FlowLayout());

	imagePanel.setBackground(Color.WHITE);

	//TODO: need to fix the path name for the 
	//at home /home/frederic/Unicore/workspace/bigdft/src/main/java/edu/kit/gridbeans/bigdft/plugin/sample.gif
	//at work /home/unicore/workspace/bigdft/src/main/java/edu/kit/gridbeans/bigdft/plugin/sample.gif
	String filename = "/home/frederic/Unicore/workspace/bigdft/src/main/java/edu/kit/gridbeans/bigdft/plugin/sample.gif";
	JLabel image = new JLabel(GuiUtils.makeImageIconFromJarredImage(filename));

	imagePanel.add(image);
	contentPane.add(BorderLayout.NORTH,imagePanel);

	JPanel msgPanel = new JPanel(new SpringLayout());

	msgPanel.setBackground(Color.WHITE);

	String msg =
	    "<html>"
	    +"<table>"
	    +"<tr><td>"
	    
	    // The Simple GridBean developer
	    + "<h2 style='font-family: \"Luxi Sans\", \"Bitstream Vera Sans\", \"Lucida Grande\", \"Trebuchet MS\", helvetica, verdana, arial, sans-serif; color: #800000'>The Simple GridBean Developper</h2>"
	    + "<p style='font-family: \"Luxi Sans\", \"Bitstream Vera Sans\", \"Lucida Grande\", \"Trebuchet MS\", helvetica, verdana, arial, sans-serif; font-size: medium'>" +
	    "The Simple GridBean has been developped: Frederic Bonnet <a href=\"fbonnet08@gmail.com\">fbonnet08@gmail.com</a><br>, <a href=\"frederic.bonnet@monash.edu\">frederic.bonnet@monash.edu</a><br>" +
	    "<br>" +
	    "The current GridBean (version :"+VersionPrime.getVersionPrimeGridBean()+") is for Pre-processing Cryo-EM images using the SIMPLE library (version: "+VersionPrime.getVersionPrime()+").</p>"

	    // Simple Developer Team
	    + "<h2 style='font-family: \"Luxi Sans\", \"Bitstream Vera Sans\", \"Lucida Grande\", \"Trebuchet MS\", helvetica, verdana, arial, sans-serif; color: #800000'>Simple Developper Team</p>"
	    + "<p style='font-family: \"Luxi Sans\", \"Bitstream Vera Sans\", \"Lucida Grande\", \"Trebuchet MS\", helvetica, verdana, arial, sans-serif; font-size: medium'>"
	    + "SIMPLE is a library used to process images from an Electron Microscope in particular from Cryo-EM that is single source electron beam microscope."
	    + "The main developer in the library is composed of Hans Elmlund, Frederic Bonnet, Cyril Reboul."
	    + "</p>"
	    
	    + "</td></tr>"
	    + "</table>"
	    
	    + "</html>";
	
	// Add the tabs to the frame's content panel
	msgPanel.add(GuiUtils.makeJLabel(msg,true));
	SpringLayoutUtils.makeCompactGrid(msgPanel, 1, 1, 4, 4, 4, 4);
	contentPane.add(BorderLayout.CENTER,msgPanel);

	//add the button in the footer
	JPanel buttonPanel = new JPanel(new FlowLayout());
	
	buttonPanel.setBackground(Color.WHITE);
	
	mButtonContinue = GuiUtils.makeJButton("Continue", mButtonListener);
	buttonPanel.add(mButtonContinue);
	
	contentPane.add(BorderLayout.SOUTH, buttonPanel);

	Dimension me = getSize();
	Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

	int left = (screen.width - me.width) / 2;
	int top = (screen.height - me.height) / 2;

	setLocation(left,top);
	
    }
    
	class ButtonListener implements ActionListener {
		public void actionPerformed(ActionEvent e)
		{
			if (e.getSource() == mButtonContinue) doContinue();
		}
	}
	protected void doContinue() {
    	mClicked = mButtonContinue;
    	setVisible(false);
    	dispose();
		
	}
    
}
