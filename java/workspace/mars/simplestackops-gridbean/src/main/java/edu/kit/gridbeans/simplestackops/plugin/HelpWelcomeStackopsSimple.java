/**
 * 
 */
package edu.kit.gridbeans.simplestackops.plugin;

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
 * @author Frederic-PC.linux64
 *
 */
public class HelpWelcomeStackopsSimple extends JDialog {

	/**
	 * Add UID
	 */
	private static final long serialVersionUID = 5486361627846234269L;
    private final static int WIDTH = 800;       /**width dimension of the window*/
    private final static int HEIGHT = 580;      /**height dimension of the window*/
	
    private JButton mButtonContinue = null;
	@SuppressWarnings("unused")
	private JButton	mClicked = null;
	private ActionListener mButtonListener = new ButtonListener();

	public HelpWelcomeStackopsSimple() {
		super();
		
		setTitle("Help Using the Welcome Panel (java)"+VersionStackops.getVersionStackopsGridBean());
		setModal(true);
		setSize(WIDTH, HEIGHT);

		Container contentPane = getContentPane();

		contentPane.setLayout(new BorderLayout());

		JPanel imagePanel = new JPanel(new FlowLayout());

		imagePanel.setBackground(Color.WHITE);

		//TODO: need to fix the path name for the 
		JLabel image = new JLabel(GuiUtils.makeImageIconFromJarredImage("./sample.gif"));
		imagePanel.add(image);
		contentPane.add(imagePanel,BorderLayout.NORTH);

		JPanel msgPanel = new JPanel(new SpringLayout());

		msgPanel.setBackground(Color.WHITE);

		String msg =
		    "<html>"
		    +"<table>"
		    +"<tr><td>"
		    
		    // The select BigDFT GridBean developper
		    + "<h2 style='font-family: \"Luxi Sans\", \"Bitstream Vera Sans\", \"Lucida Grande\", \"Trebuchet MS\", helvetica, verdana, arial, sans-serif; color: #800000'>" +
		    "What to do first?</h2>"
		    + "<p style='font-family: \"Luxi Sans\", \"Bitstream Vera Sans\", \"Lucida Grande\", \"Trebuchet MS\", helvetica, verdana, arial, sans-serif; font-size: medium'>" +
		    "First select the start you wish to do either " +
		    "" +
		    "</p>" +

		    // What happens next?
		    "<h2 style='font-family: \"Luxi Sans\", \"Bitstream Vera Sans\", \"Lucida Grande\", \"Trebuchet MS\", helvetica, verdana, arial, sans-serif; color: #800000'>" +
		    "What happens next?</p>" +
		    "<p style='font-family: \"Luxi Sans\", \"Bitstream Vera Sans\", \"Lucida Grande\", \"Trebuchet MS\", helvetica, verdana, arial, sans-serif; font-size: medium'>" +
		    "Once that you have selected the starting configuration of the preprocessing, that is if one has the Boxfile files or not," +
		    " you must go the panel in quetions and modify the " +
		    "parameters correctly" +
		    "<stackops/begin> is a program that provides standard single-particle image processing routines that are applied to MRC or SPIDER stacks. " +
"<comment/begin> You can do many things with simple_stackops. Inputting two stacks of the same size results calculation of the " +
"joint Fourier Ring  Correlation (FRC) between the images. Inputting no stacks, but setting \texttt{nptcls}, results in production of " +
"\texttt{nptcls} pure noise images, unless \texttt{ctf=yes}, then CTF images are produced. Filtering is controlled by the \texttt{hp} " +
"and \texttt{lp} arguments. Two kinds of alignments are available: shift alignment and rotational alignment. If you input an alignment " +
"document (via \texttt{oritab}) \texttt{shalgn=yes} will produce a shift-aligned stack based on the inputted orientations, whereas if you " +
"do \textit{not} input an alignment document, the alignment will be done in a reference-free manner (remember to set \texttt{trs} to some " +
"nonzero value). \texttt{roalgn=yes} behaves in the same way, but note that if an alignment document is provided \texttt{roalgn=yes} will " +
"both shift and rotate the stack according to the inputted orientations. If you want to center the images based on their center of mass, " +
"set \texttt{masscen=yes}. If you want to extract a particular state, give an alignment document (\texttt{oritab}) and set \texttt{state} " +
"to the state that you want to extract. If you want to select the fraction of best particles (according to the goal function), input an " +
"alignment doc oritab.txt and set frac=[0,1]$. You can combine the \texttt{state} and \texttt{frac} options. If you " +
"want to apply noise to images, give the desired signal-to-noise ratio via \texttt{snr}. If you want to mask your images with a spherical " +
"mask with a soft falloff, set \texttt{msk} to the radius in pixels. If you want to binarize your images, set \texttt{bin=yes}. If " +
"\texttt{tres} is defined, the images are sigmoid normalised to in range of [0,1] and threshold binarized. If \texttt{tres} is not defined the " +
"foreground/background pixels are assigned by sort-means (a variant of the continuous k-means algorithm where the initial centers are " +
"obtained by sorting the real values). If you want to calculate the autocorrelation function of your images set \texttt{acf=yes}. If you " +
"want to randomise the phases of the Fourier transforms of your images, set \texttt{phrand=yes} and \texttt{lp} to the desired low-pass " +
"limit. If you want to extract a contiguous subset of particle images from the stack, set \texttt{fromp} and \texttt{top}. If you want " +
"to fish out a number of particle images from your stack at random, set \texttt{nran} to some nonzero integer number < \texttt{nptcls}. " +
"If you want to resize you images, set the desired box \texttt{newbox} < \texttt{box} or use the \texttt{scale} option. It is often " +
"convenient to use \texttt{scale} in combination with \texttt{clip} to resize images.  If you want to normalise your images, set " +
"\texttt{norm=yes}. \texttt{hfun} controls the normalisation function. With \texttt{avg=yes} the global average of the inputted stack " +
"is calculated. With \texttt{rankify=yes} the amplitudes and phases of the FT are replaced by their ranks and the images are reverse FTed. " +
"This removes the power-law dependency of the FT. With \texttt{ctf=flip} the contrast inversions due to the CTF are corrected by the " +
"infamous (but effective) phase-flipping heuristic. This requires additional input of CTF-related parameters (\texttt{kv; fraca; cs}) " +
"as well as defocus values and astigmatism angles, communicated either via \texttt{oritab} or via \texttt{deftab}. Even if you do initially " +
"phase-flip the images, which you should do for initial model production with PRIME, you can turn on the Wiener restoration later anyway, " +
"to accomplish correct weighting of information around the CTF zeroes and maximal noise reduction. \texttt{ft2img=yes} produces images of " +
"the square power spectrum of the images in \texttt{stk}. If you define \texttt{frameavg} to some integer number > 1 averages with chunk " +
"sizes of \texttt{frameavg} are produced, which may be useful for analysis of dose-fractionated image series. \texttt{clip} can be used " +
"to re-window or pad the images to a different box size. When \texttt{compare=yes}, the two inputted stacks are Fourier ring correlated. " +
"\texttt{neg} inverts the contrast of the images by multiplication with $-1$ in Fourier space. \texttt{ctfsq} applies the squared CTF to " +
"the inputted images. \texttt{inner} is for applying an inner mask with fall-off width \texttt{width}. Finally, \texttt{merge} is for " +
"merging stack \texttt{stk2} with stack \texttt{stk}, so that the \texttt{stk2} images occur last in the series.<comment/end> <stackops/end>" +
		    "<br>" +
		    "<br>" +
		    "Once all the parameters have been correctly adjusted you may click the run/start button in the main menu.</p>"
		    
		    // About the the GriDBean
		    + "<h2 style='font-family: \"Luxi Sans\", \"Bitstream Vera Sans\", \"Lucida Grande\", \"Trebuchet MS\", helvetica, verdana, arial, sans-serif; color: #800000'>" +
		    "General use of the GridBean</p>"
		    + "<p style='font-family: \"Luxi Sans\", \"Bitstream Vera Sans\", \"Lucida Grande\", \"Trebuchet MS\", helvetica, verdana, arial, sans-serif; font-size: medium'>" +
		    "The <Strong>Boxfile</Strong> " +
		    "<br>" +
		    "<br>" +
		    "</p>"
		    		    
		    + "</td></tr>"
		    + "</table>"
		    
		    + "</html>";
		
		// Add the tabs to the frame's content panel
		msgPanel.add(GuiUtils.makeJLabel(msg,true));
		SpringLayoutUtils.makeCompactGrid(msgPanel, 1, 1, 5, 5, 5, 5);		
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


