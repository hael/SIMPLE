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
 * @author Frederic-PC.linux64
 *
 */
public class HelpWelcomePrimeSimple extends JDialog {

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

    public HelpWelcomePrimeSimple() {
        super();
		
        setTitle("Help Using the Welcome Panel (java)"+VersionPrime.getVersionPrimeGridBean());
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
            "<prime2/begin> is an \textit{ab inito} reconstruction/low-resolution refinement program based on probabilistic projection matching. " +
            "PRIME is shorthand for PRobabilistic Initial 3D Model Generation for Single-Particle Cryo-Electron Microscopy. If you process " +
            "images of a molecule with a diameter of 200 angstrom it should be possible to obtain an 8 angstrom map using 1000 reference sections " +
            "(the default setting in PRIME) provided that the images are of sufficient quality, they are many enough, and they are sampled " +
            "finely enough. We do \textit{not} recommend using more than 1000 reference sections in attempt to push for high-resolution. " +
            "Then it is more efficient to use our new continuous probabilistic refinement code, implemented by the simple_oasis " +
            "program, described below. simple_oasis is at least ten times faster than PRIME for high-resolution refinement. We do " +
            "not recommend using PRIME for heterogeneity analysis (starting off with many random blobs), because there are more effective ways " +
            "to deal with the heterogeneity problem. We will address the heterogeneity problem in the next SIMPLE 2.2 release. However, if you " +
            "suspect that you have heterogeneity of the kind where totally different species co-exist it could be worth trying initialisation " +
            "with a few blobs. You should use phase-flipped images for initial model production with PRIME (phase flipping can be done with " +
            "simple_stackops). Do \textit{not} search the origin shifts initially, when the model is of very low quality. If your " +
            "images are far off centre, use simple_stackops with option \texttt{shalgn=yes} instead to shiftalign the images " +
            "beforehand (the algorithm implemented is the same as EMAN's \texttt{cenalignint} program). We recommend running the first round " +
            "of PRIME with dynamic resolution stepping \texttt{dynlp=yes}. The \texttt{dynlp} option implements a heuristic resolution " +
            "weighting/update scheme. The initial low-pass limit is set so that each image receives ten nonzero orientation weights. When " +
            "quasi-convergence has been reached, the limit is updated one Fourier index at the time until PRIME reaches the condition where " +
            "six nonzero orientation weights is assigned to each image. FSC-based filtering is unfortunately not possible to do in the " +
            "\textit{ab initio} reconstruction step, because when the orientations are mostly random, the FSC overestimates the resolution." +
            "Once the initial model has converged, we recommend start searching the shifts (by setting \texttt{trs} to some nonzero value), " +
            "apply the FSC for resolution-weighting (by setting \texttt{eo=yes}), and turn on the Wiener restoration by setting " +
            "\texttt{ctf=yes|flip|mul} where \texttt{yes} instructs PRIME to take care of all CTF correction, \texttt{flip} indicates that " +
            "the images have been phase-flipped beforehand and \texttt{mul} indicates that the images have been multiplied with the CTF " +
            "beforehand. To use Wiener restoration you also need to input CTF parameters, for example via \texttt{deftab=defocus_values.txt}. " +
            "Remember that the defocus values should be given in microns and the astigmatism angle in degrees (one row of the file " +
            "\texttt{defocus_values.txt} may look like: \texttt{dfx=3.5  dfy=3.3  angast=20.0}). There are many ways of using (and probably " +
            "also abusing) simple_prime. We will walk you through an examples (see section \ref{recsym}, below). Use " +
            "distr_simple.pl, described below, for distributed PRIME execution. <comment/begin> Since the search is probabilistic, " +
            "we figured that an elegant convergence criterion could be formulated based on the variance of the distribution of orientations " +
            "assigned to each image. This works well for asymmetrical reconstructions, but for symmetrical reconstructions the variance " +
            "increases towards the end of the run, when the shape most consistent with the point group is being established. Note that we do " +
            "not assume any point-group symmetry in the initial runs. However, the simple_symsrch program (described below) can " +
            "be used to align the reconstruction to its symmetry axis so that future searches can be restricted to the asymmetric unit in " +
            "refinement (see section \ref{recsym}, below). Less commonly used, and less obvious input parameters are \texttt{nspace}, which " +
            "controls the number of reference projections, \texttt{amsklp}, which controls the low-pass limit used in the automask routine, " +
            "\texttt{maxits}, which controls the maximum number of iterations executed, \texttt{pgrp}, which controls the point-group symmetry, " +
            "assuming that the starting volume is aligned to its principal symmetry axis, \texttt{edge}, which controls the size of the " +
            "softening edge in the automask routine. Setting \texttt{noise=yes} produces a random noise starting volume instead of a random " +
            "blob, if no input volume is given.<comment/end> <prime2/end>" +
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


