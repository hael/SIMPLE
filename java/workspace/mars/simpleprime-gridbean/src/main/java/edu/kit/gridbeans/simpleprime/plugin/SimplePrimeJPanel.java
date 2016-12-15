/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * SimplePrimeJPanel.java
 *
 * Created on 25/11/2015, 1:39:22 PM
 */

package edu.kit.gridbeans.simpleprime.plugin;

import javax.swing.JPanel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import com.intel.gpe.gridbeans.parameters.EnvironmentVariableParameterValue;
import com.intel.gpe.gridbeans.plugins.DataSetException;
import com.intel.gpe.gridbeans.plugins.swing.controls.ComboBoxDataControl;
import com.intel.gpe.gridbeans.plugins.swing.panels.GridBeanPanel;
import com.intel.gpe.gridbeans.plugins.translators.StringValueTranslator;
import com.intel.gpe.gridbeans.plugins.validators.DoubleValueValidator;
import com.intel.gpe.gridbeans.plugins.validators.IntegerValueValidator;
import com.intel.gpe.gridbeans.plugins.validators.NotNullValidator;

import edu.kit.gridbeans.simpleprime.CSymmetrySimplePrimeCalculationType;
import edu.kit.gridbeans.simpleprime.DSymmetrySimplePrimeCalculationType;
import edu.kit.gridbeans.simpleprime.SimplePrimeGridBean;

/**
 *
 * @author Frederic Bonnet
 * @Date 25of of November 2015
 * 
 */
public class SimplePrimeJPanel extends JPanel {
    private GridBeanPanel parentPanel;
    /**
     * Add UID
     */
    private static final long serialVersionUID = 3209227781643212004L;
    private int fracslidernumber;
    private int fracValidslidernumber;
    /** Creates new form SimplePrimeJPanel */
    public SimplePrimeJPanel(GridBeanPanel parentpanel) throws DataSetException {
        this.parentPanel = parentpanel;
        initComponents();
        initCustomComponents();
        bindComponents();
    }

    private void initCustomComponents() {
        //CTF button group
        ctfbuttonGroup.add(yesCTFPrimeRadioButton);
        ctfbuttonGroup.add(noCTFPrimeRadioButton);
        ctfbuttonGroup.add(flipCTFPrimeRadioButton);
        ctfbuttonGroup.add(mulCTFPrimeRadioButton);
        //refinement button group
        refinementbuttonGroup.add(noRefinementRadioButton);
        refinementbuttonGroup.add(exhaustRefinementRadioButton);
        refinementbuttonGroup.add(shiftRefinementRadioButton);
        //dynlp button group
        dynlpbuttonGroup.add(yesDynlpRadioButton);
        dynlpbuttonGroup.add(noDynlpRadioButton);
        //noise button group
        noisebuttonGroup.add(yesNoiseRadioButton);
        noisebuttonGroup.add(noNoiseRadioButton);
        //Diversify button group
        diversifybuttonGroup.add(yesDiversifyRadioButton);
        diversifybuttonGroup.add(noDiversifyRadioButton);
        //Eo button group
        eobuttonGroup.add(yesEoRadioButton);
        eobuttonGroup.add(noEoRadioButton);
        //Norec button group
        norecbuttonGroup.add(yesNorecRadioButton);
        norecbuttonGroup.add(noNorecRadioButton);
        //Symmetry button group
        symmetryGroupbuttonGroup.add(cSymmetryRadioButton);
        symmetryGroupbuttonGroup.add(dSymmetryRadioButton);
        symmetryGroupbuttonGroup.add(tSymmetryRadioButton);
        symmetryGroupbuttonGroup.add(iSymmetryRadioButton);
        symmetryGroupbuttonGroup.add(oSymmetryRadioButton);
    }

    private void bindComponents() throws DataSetException {
        //STK
        parentPanel.linkTextField(SimplePrimeGridBean.STK, stackNameTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.STK, StringValueTranslator.getInstance());
        //Vol1
        parentPanel.linkTextField(SimplePrimeGridBean.VOL1, vol1NameTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.VOL1, StringValueTranslator.getInstance());
        //Vol2
        parentPanel.linkTextField(SimplePrimeGridBean.VOL2, vol2NameTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.VOL2, StringValueTranslator.getInstance());
        //Orientation table
        parentPanel.linkTextField(SimplePrimeGridBean.ORITAB, oriTabTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.ORITAB, StringValueTranslator.getInstance());
        //Defocus values
        parentPanel.linkTextField(SimplePrimeGridBean.DEFTAB, defTabTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.DEFTAB, StringValueTranslator.getInstance());
        //Output file for the parallel output
        parentPanel.linkTextField(SimplePrimeGridBean.OUTFILE, outfileTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.OUTFILE, StringValueTranslator.getInstance());
        //Sampling distance
        parentPanel.linkTextField(SimplePrimeGridBean.SMPD, smpdTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.SMPD, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.SMPD, new DoubleValueValidator());
        //Masking value
        parentPanel.linkTextField(SimplePrimeGridBean.MSK, mskTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.MSK, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.MSK, new DoubleValueValidator());
        //Acceleration volatage in KeV
        parentPanel.linkTextField(SimplePrimeGridBean.KV, accelerationVoltageTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.KV, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.KV, new DoubleValueValidator());
        //Shperical abberation
        parentPanel.linkTextField(SimplePrimeGridBean.CS, csTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.CS, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.CS, new DoubleValueValidator());
        //frac amp contrast{0.07}
        parentPanel.linkTextField(SimplePrimeGridBean.FRACA, fracaTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.FRACA, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.FRACA, new DoubleValueValidator());
        //High-pass limit in A
        parentPanel.linkTextField(SimplePrimeGridBean.HP, hpTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.HP, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.HP, new DoubleValueValidator());
        //Low-pass limit {20}
        parentPanel.linkTextField(SimplePrimeGridBean.LP, lpTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.LP, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.LP, new DoubleValueValidator());
        //--mw=$MW \ #[mw=<molecular weight(in kD)>]
        parentPanel.linkTextField(SimplePrimeGridBean.MW, mwTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.MW, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.MW, new DoubleValueValidator());
        //--nstates=$NSTATES \               #[nstates=<nstates to reconstruct>]
        parentPanel.linkTextField(SimplePrimeGridBean.NSTATES, nstatesTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.NSTATES, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.NSTATES, new IntegerValueValidator());
        //--startit=$STARTIT \               #[startit=<start iteration>]
        parentPanel.linkTextField(SimplePrimeGridBean.STARTIT, startitTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.STARTIT, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.STARTIT, new IntegerValueValidator());
        //--trs=$TRS \                       #[trs=<origin shift(in pixels){0}>]
        parentPanel.linkTextField(SimplePrimeGridBean.TRS, trsTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.TRS, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.TRS, new DoubleValueValidator());
        //--maxits=$MAXITS \                 #[maxits=<max iterations{100}>]
        parentPanel.linkTextField(SimplePrimeGridBean.MAXITS, maxitsTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.MAXITS, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.MAXITS, new IntegerValueValidator());
        //--find=$FIND \                     #[find=<Fourier index>]
        parentPanel.linkTextField(SimplePrimeGridBean.FIND, findTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.FIND, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.FIND, new IntegerValueValidator());
        //--lpstop=$LPSTOP \                 #[lpstop=<stay at this low-pass limit (in A)>]
        parentPanel.linkTextField(SimplePrimeGridBean.LPSTOP, lpstopTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.LPSTOP, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.LPSTOP, new DoubleValueValidator());
        //--nspace=$NSPACE \                 #[nspace=<nr reference sections{1000}>]
        parentPanel.linkTextField(SimplePrimeGridBean.NSPACE, nspaceTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.NSPACE, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.NSPACE, new IntegerValueValidator());
        //--amsklp=$AMSKLP \                 #[amsklp=<automask low-pass limit(in A)>] 
        parentPanel.linkTextField(SimplePrimeGridBean.AMSKLP, amsklpTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.AMSKLP, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.AMSKLP, new DoubleValueValidator());
        //--nnn=$NNN \                       #[nnn=<nr nearest neighbors{300}>]
        parentPanel.linkTextField(SimplePrimeGridBean.NNN, nnnTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.NNN, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.NNN, new IntegerValueValidator());
        //--lpvalid=$LPVALID \               #[lpvalid=<low-pass limit validptcls{20}>]
        parentPanel.linkTextField(SimplePrimeGridBean.LPVALID, lpvalidTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.LPVALID, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.LPVALID, new DoubleValueValidator());
        //--fstep=$FSTEP \                   #[fstep=<Fourier step size>]
        parentPanel.linkTextField(SimplePrimeGridBean.FSTEP, fstepTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.FSTEP, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.FSTEP, new IntegerValueValidator());
        //--edge=$EDGE \                     #[edge=<edge size for softening molecular envelope(in pixels){3}>] 
        parentPanel.linkTextField(SimplePrimeGridBean.EDGE, edgeTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.EDGE, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.EDGE, new IntegerValueValidator());
        //--dens=$DENS \                     #[dens=<density(e.g. 9.368 Da/A3 4 gold clusters){0.}>]
        parentPanel.linkTextField(SimplePrimeGridBean.DENS, denseTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.DENS, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.DENS, new DoubleValueValidator());
        //--trsstep=$TRSSTEP \               #[trsstep=<origin shift stepsize{1}>] 
        parentPanel.linkTextField(SimplePrimeGridBean.TRSSTEP, trsstepTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.TRSSTEP, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.TRSSTEP, new DoubleValueValidator());
        //--nvox=$NVOX \                     #[nvox=<nr of voxels in mask{0}>] 
        parentPanel.linkTextField(SimplePrimeGridBean.NVOX, nvoxTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.NVOX, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.NVOX, new IntegerValueValidator());
        //--inner=$INNER \                   #[inner=<inner mask radius(in pixels)>]
        parentPanel.linkTextField(SimplePrimeGridBean.INNER, innerTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.INNER, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.INNER, new DoubleValueValidator());
        //--time_per_image=$TIME_PER_IMAGE \ #[time_per_image=<{100}>] 
        parentPanel.linkTextField(SimplePrimeGridBean.TIME_PER_IMAGE, timePerImageTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.TIME_PER_IMAGE, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.TIME_PER_IMAGE, new DoubleValueValidator());
        //--width=$WIDTH \                   #[width=<pixels falloff inner mask{10}>]
        parentPanel.linkTextField(SimplePrimeGridBean.WIDTH, widthTextField);
        parentPanel.setValueTranslator(SimplePrimeGridBean.WIDTH, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePrimeGridBean.WIDTH, new DoubleValueValidator());
        //link to the JSlider for the frac
        class FracSliderNumberListener implements ChangeListener {
            public void stateChanged(ChangeEvent event) {
                setFracSliderNumber();
            }
        }
        ChangeListener verbosityperfslidernumberlistener = new FracSliderNumberListener();
        fracSlider.addChangeListener(verbosityperfslidernumberlistener);
        //link to the JSlider for the frac Valid
        class FracValidSliderNumberListener implements ChangeListener {
            public void stateChanged(ChangeEvent event) {
                setFracValidSliderNumber();
            }
        }
        ChangeListener fracValidslidernumberlistener = new FracValidSliderNumberListener();
        fracValidSlider.addChangeListener(fracValidslidernumberlistener);
        
    }
    /**
     * method for the JSlider assignment method for frac
     */
    public void setFracSliderNumber() {
        fracslidernumber = fracSlider.getValue();

        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue(Integer.toString(fracslidernumber));
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.FRAC, value);
    }
    /**
     * method for the JSlider assignment method for frac Valid
     */
    public void setFracValidSliderNumber() {
        fracValidslidernumber = fracValidSlider.getValue();
        
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue(Integer.toString(fracValidslidernumber));
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.FRACVALID, value);
    }
    
    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        ctfbuttonGroup = new javax.swing.ButtonGroup();
        refinementbuttonGroup = new javax.swing.ButtonGroup();
        dynlpbuttonGroup = new javax.swing.ButtonGroup();
        noisebuttonGroup = new javax.swing.ButtonGroup();
        eobuttonGroup = new javax.swing.ButtonGroup();
        norecbuttonGroup = new javax.swing.ButtonGroup();
        diversifybuttonGroup = new javax.swing.ButtonGroup();
        symmetryGroupbuttonGroup = new javax.swing.ButtonGroup();
        mainSimplePrimeScrollPanel = new javax.swing.JScrollPane();
        mainSimplePrimePanel = new javax.swing.JPanel();
        stackDetailsPanel = new javax.swing.JPanel();
        stakNameLabel = new javax.swing.JLabel();
        stackNameTextField = new javax.swing.JTextField();
        browseStkButton = new javax.swing.JButton();
        oritabLabel = new javax.swing.JLabel();
        oriTabTextField = new javax.swing.JTextField();
        browseOriTabButton = new javax.swing.JButton();
        deftabLabel = new javax.swing.JLabel();
        defTabTextField = new javax.swing.JTextField();
        browseDefTabButton = new javax.swing.JButton();
        vol1NameLabel = new javax.swing.JLabel();
        vol1NameTextField = new javax.swing.JTextField();
        browseVol1Button = new javax.swing.JButton();
        outfileLabel = new javax.swing.JLabel();
        outfileTextField = new javax.swing.JTextField();
        browseOutFileButton = new javax.swing.JButton();
        vol2NameLabel = new javax.swing.JLabel();
        vol2NameTextField = new javax.swing.JTextField();
        browseVol2Button = new javax.swing.JButton();
        microscopeDetailsPanel = new javax.swing.JPanel();
        smpdLabel = new javax.swing.JLabel();
        smpdTextField = new javax.swing.JTextField();
        accelerationVoltageLabel = new javax.swing.JLabel();
        accelerationVoltageTextField = new javax.swing.JTextField();
        csLabel = new javax.swing.JLabel();
        csTextField = new javax.swing.JTextField();
        fracaLabel = new javax.swing.JLabel();
        fracaTextField = new javax.swing.JTextField();
        mskLabel = new javax.swing.JLabel();
        mskTextField = new javax.swing.JTextField();
        fracSlider = new javax.swing.JSlider();
        symmetryGroupPanel = new javax.swing.JPanel();
        cSymmetryGroupComboBox = new javax.swing.JComboBox();
        dSymmetryGroupComboBox = new javax.swing.JComboBox();
        tSymmetryRadioButton = new javax.swing.JRadioButton();
        oSymmetryRadioButton = new javax.swing.JRadioButton();
        iSymmetryRadioButton = new javax.swing.JRadioButton();
        cSymmetryRadioButton = new javax.swing.JRadioButton();
        dSymmetryRadioButton = new javax.swing.JRadioButton();
        ctfPrimePanel = new javax.swing.JPanel();
        yesCTFPrimeRadioButton = new javax.swing.JRadioButton();
        noCTFPrimeRadioButton = new javax.swing.JRadioButton();
        flipCTFPrimeRadioButton = new javax.swing.JRadioButton();
        mulCTFPrimeRadioButton = new javax.swing.JRadioButton();
        filterPanel = new javax.swing.JPanel();
        hpLabel = new javax.swing.JLabel();
        hpTextField = new javax.swing.JTextField();
        lpLabel = new javax.swing.JLabel();
        lpTextField = new javax.swing.JTextField();
        otherParametersScrollPanel = new javax.swing.JScrollPane();
        otherParametersPanel = new javax.swing.JPanel();
        refinePanel = new javax.swing.JPanel();
        noRefinementRadioButton = new javax.swing.JRadioButton();
        exhaustRefinementRadioButton = new javax.swing.JRadioButton();
        shiftRefinementRadioButton = new javax.swing.JRadioButton();
        dynlpNoiseDiversifyPanel = new javax.swing.JPanel();
        dynlpLabel = new javax.swing.JLabel();
        yesDynlpRadioButton = new javax.swing.JRadioButton();
        noDynlpRadioButton = new javax.swing.JRadioButton();
        noiseLabel = new javax.swing.JLabel();
        yesNoiseRadioButton = new javax.swing.JRadioButton();
        noNoiseRadioButton = new javax.swing.JRadioButton();
        diversifyLabel = new javax.swing.JLabel();
        yesDiversifyRadioButton = new javax.swing.JRadioButton();
        noDiversifyRadioButton = new javax.swing.JRadioButton();
        eoLabel = new javax.swing.JLabel();
        yesEoRadioButton = new javax.swing.JRadioButton();
        noEoRadioButton = new javax.swing.JRadioButton();
        norecLabel = new javax.swing.JLabel();
        yesNorecRadioButton = new javax.swing.JRadioButton();
        noNorecRadioButton = new javax.swing.JRadioButton();
        adjustmentInOptionalPanel = new javax.swing.JPanel();
        mwLabel = new javax.swing.JLabel();
        mwTextField = new javax.swing.JTextField();
        nstatesLabel = new javax.swing.JLabel();
        nstatesTextField = new javax.swing.JTextField();
        startitLabel = new javax.swing.JLabel();
        startitTextField = new javax.swing.JTextField();
        trsLabel = new javax.swing.JLabel();
        trsTextField = new javax.swing.JTextField();
        lpstopLabel = new javax.swing.JLabel();
        lpstopTextField = new javax.swing.JTextField();
        nspaceLabel = new javax.swing.JLabel();
        nspaceTextField = new javax.swing.JTextField();
        amsklpLabel = new javax.swing.JLabel();
        amsklpTextField = new javax.swing.JTextField();
        nnnTextField = new javax.swing.JTextField();
        nnnLabel = new javax.swing.JLabel();
        maxitsLabel = new javax.swing.JLabel();
        maxitsTextField = new javax.swing.JTextField();
        lpvalidLabel = new javax.swing.JLabel();
        lpvalidTextField = new javax.swing.JTextField();
        findLabel = new javax.swing.JLabel();
        findTextField = new javax.swing.JTextField();
        fstepLabel = new javax.swing.JLabel();
        fstepTextField = new javax.swing.JTextField();
        fracValidSlider = new javax.swing.JSlider();
        EdgeOtherPanel = new javax.swing.JPanel();
        edgeLabel = new javax.swing.JLabel();
        edgeTextField = new javax.swing.JTextField();
        denseLabel = new javax.swing.JLabel();
        denseTextField = new javax.swing.JTextField();
        trsstepLabel = new javax.swing.JLabel();
        trsstepTextField = new javax.swing.JTextField();
        nvoxLabel = new javax.swing.JLabel();
        nvoxTextField = new javax.swing.JTextField();
        innerLabel = new javax.swing.JLabel();
        innerTextField = new javax.swing.JTextField();
        timePerImageLabel = new javax.swing.JLabel();
        timePerImageTextField = new javax.swing.JTextField();
        widthLabel = new javax.swing.JLabel();
        widthTextField = new javax.swing.JTextField();

        setLayout(new java.awt.BorderLayout());

        stackDetailsPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Stack and input details..."));

        stakNameLabel.setText("Stk name:");

        browseStkButton.setText("Browse");

        oritabLabel.setText("Ori Table:");

        browseOriTabButton.setText("Browse");

        deftabLabel.setText("Def Table:");

        browseDefTabButton.setText("Browse");

        vol1NameLabel.setText("Vol1 name:");

        browseVol1Button.setText("Browse");

        outfileLabel.setText("Outfile:");

        browseOutFileButton.setText("Browse");

        vol2NameLabel.setText("Vol2 name:");

        browseVol2Button.setText("Browse");

        javax.swing.GroupLayout stackDetailsPanelLayout = new javax.swing.GroupLayout(stackDetailsPanel);
        stackDetailsPanel.setLayout(stackDetailsPanelLayout);
        stackDetailsPanelLayout.setHorizontalGroup(
            stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(stackDetailsPanelLayout.createSequentialGroup()
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(stackDetailsPanelLayout.createSequentialGroup()
                        .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(stakNameLabel)
                            .addComponent(vol1NameLabel))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(vol1NameTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 176, Short.MAX_VALUE)
                            .addComponent(stackNameTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 176, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(browseStkButton)
                            .addComponent(browseVol1Button)))
                    .addGroup(stackDetailsPanelLayout.createSequentialGroup()
                        .addComponent(vol2NameLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(vol2NameTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 176, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(browseVol2Button))
                    .addGroup(stackDetailsPanelLayout.createSequentialGroup()
                        .addComponent(oritabLabel)
                        .addGap(18, 18, 18)
                        .addComponent(oriTabTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 177, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(browseOriTabButton))
                    .addGroup(stackDetailsPanelLayout.createSequentialGroup()
                        .addComponent(deftabLabel)
                        .addGap(18, 18, 18)
                        .addComponent(defTabTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 175, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(browseDefTabButton))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, stackDetailsPanelLayout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(outfileLabel)
                        .addGap(18, 18, 18)
                        .addComponent(outfileTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 174, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(browseOutFileButton)))
                .addContainerGap())
        );
        stackDetailsPanelLayout.setVerticalGroup(
            stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(stackDetailsPanelLayout.createSequentialGroup()
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(stackNameTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(browseStkButton)
                    .addComponent(stakNameLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(vol1NameTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(vol1NameLabel)
                    .addComponent(browseVol1Button))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(vol2NameTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(vol2NameLabel)
                    .addComponent(browseVol2Button))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(oritabLabel)
                    .addComponent(browseOriTabButton)
                    .addComponent(oriTabTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(deftabLabel)
                    .addComponent(defTabTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(browseDefTabButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(outfileLabel)
                    .addComponent(outfileTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(browseOutFileButton))
                .addContainerGap(13, Short.MAX_VALUE))
        );

        oriTabTextField.getAccessibleContext().setAccessibleDescription("Previous alignment doc");
        defTabTextField.getAccessibleContext().setAccessibleDescription("Text file defocus values");
        outfileTextField.getAccessibleContext().setAccessibleDescription("Output alignment doc 4 parallell jobs");

        microscopeDetailsPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("MIcroscope details..."));

        smpdLabel.setText("Sampling distance smpd:");

        accelerationVoltageLabel.setText("acceleration voltage(in kV) kv:");

        csLabel.setText("Spherical Aberration cs:");

        fracaLabel.setText("Ampl. Contrast fraca:");

        mskLabel.setText("Mask msk:");

        javax.swing.GroupLayout microscopeDetailsPanelLayout = new javax.swing.GroupLayout(microscopeDetailsPanel);
        microscopeDetailsPanel.setLayout(microscopeDetailsPanelLayout);
        microscopeDetailsPanelLayout.setHorizontalGroup(
            microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(microscopeDetailsPanelLayout.createSequentialGroup()
                .addGroup(microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(microscopeDetailsPanelLayout.createSequentialGroup()
                        .addComponent(smpdLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(smpdTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 168, Short.MAX_VALUE))
                    .addGroup(microscopeDetailsPanelLayout.createSequentialGroup()
                        .addComponent(accelerationVoltageLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(accelerationVoltageTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 139, Short.MAX_VALUE))
                    .addGroup(microscopeDetailsPanelLayout.createSequentialGroup()
                        .addComponent(csLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(csTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 173, Short.MAX_VALUE))
                    .addGroup(microscopeDetailsPanelLayout.createSequentialGroup()
                        .addGroup(microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(fracaLabel)
                            .addComponent(mskLabel))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(mskTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 184, Short.MAX_VALUE)
                            .addComponent(fracaTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 184, Short.MAX_VALUE))))
                .addContainerGap())
        );
        microscopeDetailsPanelLayout.setVerticalGroup(
            microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(microscopeDetailsPanelLayout.createSequentialGroup()
                .addGroup(microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(smpdLabel)
                    .addComponent(smpdTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(accelerationVoltageLabel)
                    .addComponent(accelerationVoltageTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(csTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(csLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fracaTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(fracaLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(mskTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(mskLabel))
                .addContainerGap(55, Short.MAX_VALUE))
        );

        smpdTextField.getAccessibleContext().setAccessibleDescription("Sampling distance(in A)");
        accelerationVoltageTextField.getAccessibleContext().setAccessibleDescription("Acceleration voltage(in kV){300.}");
        csTextField.getAccessibleContext().setAccessibleDescription("Spherical aberration constant(in mm){2.7}");
        fracaTextField.getAccessibleContext().setAccessibleDescription("Frac amp contrast{0.07}");
        mskTextField.getAccessibleContext().setAccessibleDescription("Mask radius(in pixels)");

        fracSlider.setFont(new java.awt.Font("Dialog", 0, 8));
        fracSlider.setMajorTickSpacing(10);
        fracSlider.setMinorTickSpacing(2);
        fracSlider.setPaintLabels(true);
        fracSlider.setPaintTicks(true);
        fracSlider.setSnapToTicks(true);
        fracSlider.setToolTipText("Fraction of particles to include{1}");
        fracSlider.setValue(100);
        fracSlider.setBorder(javax.swing.BorderFactory.createTitledBorder("fraction of ptcls to include..."));
        fracSlider.setValueIsAdjusting(true);

        symmetryGroupPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("symmetry Group <cn|dn|t|o|i>..."));

        cSymmetryGroupComboBox.setMaximumRowCount(10);
        cSymmetryGroupComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9" }));
        cSymmetryGroupComboBox.setBorder(javax.swing.BorderFactory.createTitledBorder("C..."));
        cSymmetryGroupComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cSymmetryGroupComboBoxActionPerformed(evt);
            }
        });

        dSymmetryGroupComboBox.setMaximumRowCount(10);
        dSymmetryGroupComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9" }));
        dSymmetryGroupComboBox.setBorder(javax.swing.BorderFactory.createTitledBorder("D..."));
        dSymmetryGroupComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                dSymmetryGroupComboBoxActionPerformed(evt);
            }
        });

        tSymmetryRadioButton.setText("t");
        tSymmetryRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                tSymmetryRadioButtonActionPerformed(evt);
            }
        });

        oSymmetryRadioButton.setText("o");
        oSymmetryRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                oSymmetryRadioButtonActionPerformed(evt);
            }
        });

        iSymmetryRadioButton.setText("i");
        iSymmetryRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                iSymmetryRadioButtonActionPerformed(evt);
            }
        });

        cSymmetryRadioButton.setText("C Sym");
        cSymmetryRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cSymmetryRadioButtonActionPerformed(evt);
            }
        });

        dSymmetryRadioButton.setText("D");
        dSymmetryRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                dSymmetryRadioButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout symmetryGroupPanelLayout = new javax.swing.GroupLayout(symmetryGroupPanel);
        symmetryGroupPanel.setLayout(symmetryGroupPanelLayout);
        symmetryGroupPanelLayout.setHorizontalGroup(
            symmetryGroupPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(symmetryGroupPanelLayout.createSequentialGroup()
                .addGroup(symmetryGroupPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(cSymmetryRadioButton)
                    .addComponent(tSymmetryRadioButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(cSymmetryGroupComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, 78, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(symmetryGroupPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(dSymmetryRadioButton)
                    .addComponent(iSymmetryRadioButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(dSymmetryGroupComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, 78, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(oSymmetryRadioButton))
        );
        symmetryGroupPanelLayout.setVerticalGroup(
            symmetryGroupPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(symmetryGroupPanelLayout.createSequentialGroup()
                .addGroup(symmetryGroupPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(symmetryGroupPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(cSymmetryGroupComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, 46, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(dSymmetryGroupComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, 46, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(symmetryGroupPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                        .addGroup(javax.swing.GroupLayout.Alignment.LEADING, symmetryGroupPanelLayout.createSequentialGroup()
                            .addComponent(dSymmetryRadioButton)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(iSymmetryRadioButton))
                        .addGroup(javax.swing.GroupLayout.Alignment.LEADING, symmetryGroupPanelLayout.createSequentialGroup()
                            .addComponent(cSymmetryRadioButton)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(tSymmetryRadioButton)))
                    .addComponent(oSymmetryRadioButton))
                .addContainerGap(12, Short.MAX_VALUE))
        );

        ctfPrimePanel.setBorder(javax.swing.BorderFactory.createTitledBorder("CTF setup..."));

        yesCTFPrimeRadioButton.setText("Yes");
        yesCTFPrimeRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesCTFPrimeRadioButtonActionPerformed(evt);
            }
        });

        noCTFPrimeRadioButton.setText("No");
        noCTFPrimeRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noCTFPrimeRadioButtonActionPerformed(evt);
            }
        });

        flipCTFPrimeRadioButton.setText("Flip");
        flipCTFPrimeRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                flipCTFPrimeRadioButtonActionPerformed(evt);
            }
        });

        mulCTFPrimeRadioButton.setText("Mul");
        mulCTFPrimeRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                mulCTFPrimeRadioButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout ctfPrimePanelLayout = new javax.swing.GroupLayout(ctfPrimePanel);
        ctfPrimePanel.setLayout(ctfPrimePanelLayout);
        ctfPrimePanelLayout.setHorizontalGroup(
            ctfPrimePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(ctfPrimePanelLayout.createSequentialGroup()
                .addComponent(yesCTFPrimeRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(noCTFPrimeRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(flipCTFPrimeRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(mulCTFPrimeRadioButton)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        ctfPrimePanelLayout.setVerticalGroup(
            ctfPrimePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(ctfPrimePanelLayout.createSequentialGroup()
                .addGroup(ctfPrimePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(yesCTFPrimeRadioButton)
                    .addComponent(noCTFPrimeRadioButton)
                    .addComponent(flipCTFPrimeRadioButton)
                    .addComponent(mulCTFPrimeRadioButton))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        filterPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Filters..."));

        hpLabel.setText("High pass limit hp:");

        lpLabel.setText("Low pass limit lp:");

        javax.swing.GroupLayout filterPanelLayout = new javax.swing.GroupLayout(filterPanel);
        filterPanel.setLayout(filterPanelLayout);
        filterPanelLayout.setHorizontalGroup(
            filterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(filterPanelLayout.createSequentialGroup()
                .addComponent(hpLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(hpTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(lpLabel)
                .addGap(18, 18, 18)
                .addComponent(lpTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 49, Short.MAX_VALUE)
                .addGap(62, 62, 62))
        );
        filterPanelLayout.setVerticalGroup(
            filterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(filterPanelLayout.createSequentialGroup()
                .addGroup(filterPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(hpLabel)
                    .addComponent(hpTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(lpLabel)
                    .addComponent(lpTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(20, Short.MAX_VALUE))
        );

        otherParametersPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Optional parameters..."));

        refinePanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Refinement..."));

        noRefinementRadioButton.setText("No");
        noRefinementRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noRefinementRadioButtonActionPerformed(evt);
            }
        });

        exhaustRefinementRadioButton.setText("Exhaust");
        exhaustRefinementRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exhaustRefinementRadioButtonActionPerformed(evt);
            }
        });

        shiftRefinementRadioButton.setText("Shift");
        shiftRefinementRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                shiftRefinementRadioButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout refinePanelLayout = new javax.swing.GroupLayout(refinePanel);
        refinePanel.setLayout(refinePanelLayout);
        refinePanelLayout.setHorizontalGroup(
            refinePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(refinePanelLayout.createSequentialGroup()
                .addComponent(noRefinementRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(exhaustRefinementRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(shiftRefinementRadioButton)
                .addContainerGap(23, Short.MAX_VALUE))
        );
        refinePanelLayout.setVerticalGroup(
            refinePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(refinePanelLayout.createSequentialGroup()
                .addGroup(refinePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(noRefinementRadioButton)
                    .addComponent(exhaustRefinementRadioButton)
                    .addComponent(shiftRefinementRadioButton))
                .addContainerGap(42, Short.MAX_VALUE))
        );

        dynlpNoiseDiversifyPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Diversify Noise and dynlp panel..."));

        dynlpLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        dynlpLabel.setText("dynlp:");

        yesDynlpRadioButton.setText("Yes");
        yesDynlpRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesDynlpRadioButtonActionPerformed(evt);
            }
        });

        noDynlpRadioButton.setText("No");
        noDynlpRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noDynlpRadioButtonActionPerformed(evt);
            }
        });

        noiseLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        noiseLabel.setText("Noise:");

        yesNoiseRadioButton.setText("Yes");
        yesNoiseRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesNoiseRadioButtonActionPerformed(evt);
            }
        });

        noNoiseRadioButton.setText("No");
        noNoiseRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noNoiseRadioButtonActionPerformed(evt);
            }
        });

        diversifyLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        diversifyLabel.setText("Diversify:");

        yesDiversifyRadioButton.setText("Yes");
        yesDiversifyRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesDiversifyRadioButtonActionPerformed(evt);
            }
        });

        noDiversifyRadioButton.setText("No");
        noDiversifyRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noDiversifyRadioButtonActionPerformed(evt);
            }
        });

        eoLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        eoLabel.setText("eo:");

        yesEoRadioButton.setText("Yes");
        yesEoRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesEoRadioButtonActionPerformed(evt);
            }
        });

        noEoRadioButton.setText("No");
        noEoRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noEoRadioButtonActionPerformed(evt);
            }
        });

        norecLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        norecLabel.setText("norec:");

        yesNorecRadioButton.setText("Yes");
        yesNorecRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesNorecRadioButtonActionPerformed(evt);
            }
        });

        noNorecRadioButton.setText("No");
        noNorecRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noNorecRadioButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout dynlpNoiseDiversifyPanelLayout = new javax.swing.GroupLayout(dynlpNoiseDiversifyPanel);
        dynlpNoiseDiversifyPanel.setLayout(dynlpNoiseDiversifyPanelLayout);
        dynlpNoiseDiversifyPanelLayout.setHorizontalGroup(
            dynlpNoiseDiversifyPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(dynlpNoiseDiversifyPanelLayout.createSequentialGroup()
                .addGroup(dynlpNoiseDiversifyPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(dynlpLabel)
                    .addComponent(eoLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(dynlpNoiseDiversifyPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(dynlpNoiseDiversifyPanelLayout.createSequentialGroup()
                        .addComponent(yesDynlpRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noDynlpRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noiseLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(yesNoiseRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noNoiseRadioButton))
                    .addGroup(dynlpNoiseDiversifyPanelLayout.createSequentialGroup()
                        .addComponent(yesEoRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noEoRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(norecLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(yesNorecRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noNorecRadioButton)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(diversifyLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(yesDiversifyRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(noDiversifyRadioButton)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        dynlpNoiseDiversifyPanelLayout.setVerticalGroup(
            dynlpNoiseDiversifyPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(dynlpNoiseDiversifyPanelLayout.createSequentialGroup()
                .addGroup(dynlpNoiseDiversifyPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(dynlpLabel)
                    .addComponent(yesDynlpRadioButton)
                    .addComponent(noDynlpRadioButton)
                    .addComponent(noiseLabel)
                    .addComponent(yesNoiseRadioButton)
                    .addComponent(noNoiseRadioButton)
                    .addComponent(diversifyLabel)
                    .addComponent(yesDiversifyRadioButton)
                    .addComponent(noDiversifyRadioButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(dynlpNoiseDiversifyPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(eoLabel)
                    .addComponent(yesEoRadioButton)
                    .addComponent(noEoRadioButton)
                    .addComponent(norecLabel)
                    .addComponent(yesNorecRadioButton)
                    .addComponent(noNorecRadioButton))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        adjustmentInOptionalPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Adjustment parameters..."));

        mwLabel.setText("mw:");

        nstatesLabel.setText("nstates:");

        startitLabel.setText("startit:");

        trsLabel.setText("trs<origin>:");

        lpstopLabel.setText("lpstop:");

        nspaceLabel.setText("nspace:");

        amsklpLabel.setText("amsklp:");

        nnnLabel.setText("nnn:");

        maxitsLabel.setText("maxits:");

        lpvalidLabel.setText("lpvalid:");

        findLabel.setText("find:");

        fstepLabel.setText("fstep:");

        javax.swing.GroupLayout adjustmentInOptionalPanelLayout = new javax.swing.GroupLayout(adjustmentInOptionalPanel);
        adjustmentInOptionalPanel.setLayout(adjustmentInOptionalPanelLayout);
        adjustmentInOptionalPanelLayout.setHorizontalGroup(
            adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(adjustmentInOptionalPanelLayout.createSequentialGroup()
                .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, adjustmentInOptionalPanelLayout.createSequentialGroup()
                        .addComponent(lpstopLabel)
                        .addGap(18, 18, 18)
                        .addComponent(lpstopTextField))
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, adjustmentInOptionalPanelLayout.createSequentialGroup()
                        .addComponent(mwLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(mwTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(nstatesLabel)
                    .addComponent(nspaceLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(nspaceTextField)
                    .addComponent(nstatesTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 42, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(startitLabel)
                    .addComponent(amsklpLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(amsklpTextField)
                    .addComponent(startitTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(trsLabel)
                    .addComponent(nnnLabel))
                .addGap(2, 2, 2)
                .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(nnnTextField)
                    .addComponent(trsTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 39, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(adjustmentInOptionalPanelLayout.createSequentialGroup()
                        .addComponent(maxitsLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(maxitsTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(adjustmentInOptionalPanelLayout.createSequentialGroup()
                        .addComponent(lpvalidLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(lpvalidTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(adjustmentInOptionalPanelLayout.createSequentialGroup()
                        .addComponent(findLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(findTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(adjustmentInOptionalPanelLayout.createSequentialGroup()
                        .addComponent(fstepLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fstepTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(63, Short.MAX_VALUE))
        );
        adjustmentInOptionalPanelLayout.setVerticalGroup(
            adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(adjustmentInOptionalPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(adjustmentInOptionalPanelLayout.createSequentialGroup()
                        .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(findLabel)
                            .addComponent(findTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(fstepLabel)
                            .addComponent(fstepTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addGroup(adjustmentInOptionalPanelLayout.createSequentialGroup()
                        .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(mwLabel)
                            .addComponent(mwTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(nstatesLabel)
                            .addComponent(nstatesTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(startitLabel)
                            .addComponent(startitTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(trsLabel)
                            .addComponent(trsTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(maxitsLabel)
                            .addComponent(maxitsTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(adjustmentInOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(lpstopLabel)
                            .addComponent(lpstopTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(nspaceLabel)
                            .addComponent(nspaceTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(amsklpLabel)
                            .addComponent(amsklpTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(nnnLabel)
                            .addComponent(nnnTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(lpvalidLabel)
                            .addComponent(lpvalidTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))))
                .addContainerGap(17, Short.MAX_VALUE))
        );

        mwTextField.getAccessibleContext().setAccessibleDescription("molecular weight(in kD)");
        nstatesTextField.getAccessibleContext().setAccessibleDescription("nstates to reconstruct");
        startitTextField.getAccessibleContext().setAccessibleDescription("Start iteration");
        trsTextField.getAccessibleContext().setAccessibleDescription("Origin shift(in pixels){0}");
        lpstopTextField.getAccessibleContext().setAccessibleDescription("Stay at this low-pass limit (in A)");
        nspaceTextField.getAccessibleContext().setAccessibleDescription("nr reference sections{1000}");
        amsklpTextField.getAccessibleContext().setAccessibleDescription("Automask low-pass limit(in A)");
        nnnTextField.getAccessibleContext().setAccessibleDescription("nr nearest neighbors{300}");
        maxitsTextField.getAccessibleContext().setAccessibleDescription("Max iterations{100}");
        lpvalidTextField.getAccessibleContext().setAccessibleDescription("Low-pass limit validptcls{20}");
        findTextField.getAccessibleContext().setAccessibleDescription("Fourier index");
        fstepTextField.getAccessibleContext().setAccessibleDescription("Fourier step size");

        fracValidSlider.setFont(new java.awt.Font("Dialog", 0, 8));
        fracValidSlider.setMajorTickSpacing(10);
        fracValidSlider.setMinorTickSpacing(2);
        fracValidSlider.setPaintLabels(true);
        fracValidSlider.setPaintTicks(true);
        fracValidSlider.setSnapToTicks(true);
        fracValidSlider.setToolTipText("Fraction of particles to include{1}");
        fracValidSlider.setValue(0);
        fracValidSlider.setBorder(javax.swing.BorderFactory.createTitledBorder("Fraction of particles 4 validation..."));
        fracValidSlider.setValueIsAdjusting(true);

        EdgeOtherPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Edge, time/image, dens, trs step, nvox, inner and width panel... "));

        edgeLabel.setText("edge:");

        denseLabel.setText("dens:");

        trsstepLabel.setText("trsstep:");

        nvoxLabel.setText("nvox:");

        innerLabel.setText("inner:");

        timePerImageLabel.setText("time_per_image:");

        widthLabel.setText("width:");

        javax.swing.GroupLayout EdgeOtherPanelLayout = new javax.swing.GroupLayout(EdgeOtherPanel);
        EdgeOtherPanel.setLayout(EdgeOtherPanelLayout);
        EdgeOtherPanelLayout.setHorizontalGroup(
            EdgeOtherPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(EdgeOtherPanelLayout.createSequentialGroup()
                .addGroup(EdgeOtherPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(edgeLabel)
                    .addComponent(widthLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(EdgeOtherPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addComponent(widthTextField, javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(edgeTextField, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, 50, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(denseLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(denseTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 42, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(trsstepLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(trsstepTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(nvoxLabel)
                .addGap(2, 2, 2)
                .addComponent(nvoxTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 39, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(innerLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(innerTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(timePerImageLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(timePerImageTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 36, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(53, Short.MAX_VALUE))
        );
        EdgeOtherPanelLayout.setVerticalGroup(
            EdgeOtherPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(EdgeOtherPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(EdgeOtherPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(edgeLabel)
                    .addComponent(edgeTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(denseLabel)
                    .addComponent(denseTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(trsstepLabel)
                    .addComponent(trsstepTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(nvoxLabel)
                    .addComponent(nvoxTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(innerLabel)
                    .addComponent(innerTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(timePerImageLabel)
                    .addComponent(timePerImageTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(EdgeOtherPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(widthLabel)
                    .addComponent(widthTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(36, Short.MAX_VALUE))
        );

        denseTextField.getAccessibleContext().setAccessibleDescription("Density(e.g. 9.368 Da/A3 4 gold clusters){0.}");
        trsstepTextField.getAccessibleContext().setAccessibleDescription("Origin shift stepsize{1}");
        nvoxTextField.getAccessibleContext().setAccessibleDescription("Number of voxels in mask{0}");
        innerTextField.getAccessibleContext().setAccessibleDescription("Inner mask radius(in pixels)");
        timePerImageTextField.getAccessibleContext().setAccessibleDescription("Time_per_image=<{100}");
        widthTextField.getAccessibleContext().setAccessibleDescription("Pixels falloff inner mask{10}");

        javax.swing.GroupLayout otherParametersPanelLayout = new javax.swing.GroupLayout(otherParametersPanel);
        otherParametersPanel.setLayout(otherParametersPanelLayout);
        otherParametersPanelLayout.setHorizontalGroup(
            otherParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(otherParametersPanelLayout.createSequentialGroup()
                .addGroup(otherParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, otherParametersPanelLayout.createSequentialGroup()
                        .addComponent(refinePanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(dynlpNoiseDiversifyPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, otherParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                        .addComponent(fracValidSlider, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(EdgeOtherPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(adjustmentInOptionalPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .addContainerGap(133, Short.MAX_VALUE))
        );
        otherParametersPanelLayout.setVerticalGroup(
            otherParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(otherParametersPanelLayout.createSequentialGroup()
                .addGroup(otherParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(dynlpNoiseDiversifyPanel, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(refinePanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 12, Short.MAX_VALUE)
                .addComponent(adjustmentInOptionalPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(fracValidSlider, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(EdgeOtherPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(71, 71, 71))
        );

        fracValidSlider.getAccessibleContext().setAccessibleDescription("Fraction of particles 4 validation{0.}");

        otherParametersScrollPanel.setViewportView(otherParametersPanel);

        javax.swing.GroupLayout mainSimplePrimePanelLayout = new javax.swing.GroupLayout(mainSimplePrimePanel);
        mainSimplePrimePanel.setLayout(mainSimplePrimePanelLayout);
        mainSimplePrimePanelLayout.setHorizontalGroup(
            mainSimplePrimePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mainSimplePrimePanelLayout.createSequentialGroup()
                .addGroup(mainSimplePrimePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(mainSimplePrimePanelLayout.createSequentialGroup()
                        .addGroup(mainSimplePrimePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                            .addComponent(fracSlider, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(stackDetailsPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(mainSimplePrimePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(symmetryGroupPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(microscopeDetailsPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addGroup(mainSimplePrimePanelLayout.createSequentialGroup()
                        .addComponent(ctfPrimePanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(filterPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(otherParametersScrollPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 721, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(269, Short.MAX_VALUE))
        );
        mainSimplePrimePanelLayout.setVerticalGroup(
            mainSimplePrimePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mainSimplePrimePanelLayout.createSequentialGroup()
                .addGroup(mainSimplePrimePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addComponent(microscopeDetailsPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(stackDetailsPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addGroup(mainSimplePrimePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(mainSimplePrimePanelLayout.createSequentialGroup()
                        .addGap(13, 13, 13)
                        .addComponent(fracSlider, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(mainSimplePrimePanelLayout.createSequentialGroup()
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(symmetryGroupPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .addGroup(mainSimplePrimePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(mainSimplePrimePanelLayout.createSequentialGroup()
                        .addGap(8, 8, 8)
                        .addComponent(ctfPrimePanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(mainSimplePrimePanelLayout.createSequentialGroup()
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(filterPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(otherParametersScrollPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 431, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(180, Short.MAX_VALUE))
        );

        mainSimplePrimeScrollPanel.setViewportView(mainSimplePrimePanel);

        add(mainSimplePrimeScrollPanel, java.awt.BorderLayout.CENTER);
    }// </editor-fold>//GEN-END:initComponents

    private void cSymmetryGroupComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cSymmetryGroupComboBoxActionPerformed
        // TODO add your handling code here:
}//GEN-LAST:event_cSymmetryGroupComboBoxActionPerformed

    private void dSymmetryGroupComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_dSymmetryGroupComboBoxActionPerformed
        // TODO add your handling code here:
}//GEN-LAST:event_dSymmetryGroupComboBoxActionPerformed

    private void tSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_tSymmetryRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("t");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.PGRP, value);
}//GEN-LAST:event_tSymmetryRadioButtonActionPerformed
    private void oSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_oSymmetryRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("o");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.PGRP, value);
}//GEN-LAST:event_oSymmetryRadioButtonActionPerformed
    private void iSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_iSymmetryRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("i");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.PGRP, value);
}//GEN-LAST:event_iSymmetryRadioButtonActionPerformed
    private void yesCTFPrimeRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesCTFPrimeRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.CTF, value);
}//GEN-LAST:event_yesCTFPrimeRadioButtonActionPerformed
    private void noCTFPrimeRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noCTFPrimeRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.CTF, value);
}//GEN-LAST:event_noCTFPrimeRadioButtonActionPerformed
    private void flipCTFPrimeRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_flipCTFPrimeRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("flip");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.CTF, value);
}//GEN-LAST:event_flipCTFPrimeRadioButtonActionPerformed
    private void mulCTFPrimeRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_mulCTFPrimeRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("mul");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.CTF, value);
}//GEN-LAST:event_mulCTFPrimeRadioButtonActionPerformed
    private void noRefinementRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noRefinementRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.REFINE, value);
    }//GEN-LAST:event_noRefinementRadioButtonActionPerformed
    private void exhaustRefinementRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_exhaustRefinementRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("exhaust");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.REFINE, value);
    }//GEN-LAST:event_exhaustRefinementRadioButtonActionPerformed
    private void shiftRefinementRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_shiftRefinementRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("shift");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.REFINE, value);
    }//GEN-LAST:event_shiftRefinementRadioButtonActionPerformed
    private void yesDynlpRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesDynlpRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.DYNLP, value);
    }//GEN-LAST:event_yesDynlpRadioButtonActionPerformed
    private void noDynlpRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noDynlpRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.DYNLP, value);
    }//GEN-LAST:event_noDynlpRadioButtonActionPerformed
    private void yesNoiseRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesNoiseRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.NOISE, value);
    }//GEN-LAST:event_yesNoiseRadioButtonActionPerformed
    private void noNoiseRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noNoiseRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.NOISE, value);
    }//GEN-LAST:event_noNoiseRadioButtonActionPerformed
    private void yesDiversifyRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesDiversifyRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.DIVERSIFY, value);
    }//GEN-LAST:event_yesDiversifyRadioButtonActionPerformed
    private void noDiversifyRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noDiversifyRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.DIVERSIFY, value);
    }//GEN-LAST:event_noDiversifyRadioButtonActionPerformed
    private void yesEoRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesEoRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.EO, value);
    }//GEN-LAST:event_yesEoRadioButtonActionPerformed
    private void noEoRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noEoRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.EO, value);
    }//GEN-LAST:event_noEoRadioButtonActionPerformed
    private void yesNorecRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesNorecRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.NOREC, value);
    }//GEN-LAST:event_yesNorecRadioButtonActionPerformed
    private void noNorecRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noNorecRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimplePrimeGridBean.NOREC, value);
    }//GEN-LAST:event_noNorecRadioButtonActionPerformed
    private void cSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt){//GEN-FIRST:event_cSymmetryRadioButtonActionPerformed
        //C symmetry combo box
        try {
            ComboBoxDataControl cSymComboBoxDataControl = new ComboBoxDataControl(parentPanel.getClient(), SimplePrimeGridBean.PGRP, cSymmetryGroupComboBox);
            cSymComboBoxDataControl.setPossibleValues(CSymmetrySimplePrimeCalculationType.getAllTypes());
            parentPanel.linkDataControl(SimplePrimeGridBean.PGRP, cSymComboBoxDataControl);
            parentPanel.setValueTranslator(SimplePrimeGridBean.PGRP, StringValueTranslator.getInstance());
            parentPanel.setValueValidator(SimplePrimeGridBean.PGRP, NotNullValidator.getInstance());
        } catch (DataSetException e) {
            e.printStackTrace();
        }
    }//GEN-LAST:event_cSymmetryRadioButtonActionPerformed
    private void dSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_dSymmetryRadioButtonActionPerformed
        try {
            //D symmetry combo box
            ComboBoxDataControl dSymComboBoxDataControl = new ComboBoxDataControl(parentPanel.getClient(), SimplePrimeGridBean.PGRP, dSymmetryGroupComboBox);
            dSymComboBoxDataControl.setPossibleValues(DSymmetrySimplePrimeCalculationType.getAllTypes());
            parentPanel.linkDataControl(SimplePrimeGridBean.PGRP, dSymComboBoxDataControl);
            parentPanel.setValueTranslator(SimplePrimeGridBean.PGRP, StringValueTranslator.getInstance());
            parentPanel.setValueValidator(SimplePrimeGridBean.PGRP, NotNullValidator.getInstance());
        } catch (DataSetException e) {
            e.printStackTrace();
        }        
    }//GEN-LAST:event_dSymmetryRadioButtonActionPerformed


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JPanel EdgeOtherPanel;
    private javax.swing.JLabel accelerationVoltageLabel;
    private javax.swing.JTextField accelerationVoltageTextField;
    private javax.swing.JPanel adjustmentInOptionalPanel;
    private javax.swing.JLabel amsklpLabel;
    private javax.swing.JTextField amsklpTextField;
    private javax.swing.JButton browseDefTabButton;
    private javax.swing.JButton browseOriTabButton;
    private javax.swing.JButton browseOutFileButton;
    private javax.swing.JButton browseStkButton;
    private javax.swing.JButton browseVol1Button;
    private javax.swing.JButton browseVol2Button;
    private javax.swing.JComboBox cSymmetryGroupComboBox;
    private javax.swing.JRadioButton cSymmetryRadioButton;
    private javax.swing.JLabel csLabel;
    private javax.swing.JTextField csTextField;
    private javax.swing.JPanel ctfPrimePanel;
    private javax.swing.ButtonGroup ctfbuttonGroup;
    private javax.swing.JComboBox dSymmetryGroupComboBox;
    private javax.swing.JRadioButton dSymmetryRadioButton;
    private javax.swing.JTextField defTabTextField;
    private javax.swing.JLabel deftabLabel;
    private javax.swing.JLabel denseLabel;
    private javax.swing.JTextField denseTextField;
    private javax.swing.JLabel diversifyLabel;
    private javax.swing.ButtonGroup diversifybuttonGroup;
    private javax.swing.JLabel dynlpLabel;
    private javax.swing.JPanel dynlpNoiseDiversifyPanel;
    private javax.swing.ButtonGroup dynlpbuttonGroup;
    private javax.swing.JLabel edgeLabel;
    private javax.swing.JTextField edgeTextField;
    private javax.swing.JLabel eoLabel;
    private javax.swing.ButtonGroup eobuttonGroup;
    private javax.swing.JRadioButton exhaustRefinementRadioButton;
    private javax.swing.JPanel filterPanel;
    private javax.swing.JLabel findLabel;
    private javax.swing.JTextField findTextField;
    private javax.swing.JRadioButton flipCTFPrimeRadioButton;
    private javax.swing.JSlider fracSlider;
    private javax.swing.JSlider fracValidSlider;
    private javax.swing.JLabel fracaLabel;
    private javax.swing.JTextField fracaTextField;
    private javax.swing.JLabel fstepLabel;
    private javax.swing.JTextField fstepTextField;
    private javax.swing.JLabel hpLabel;
    private javax.swing.JTextField hpTextField;
    private javax.swing.JRadioButton iSymmetryRadioButton;
    private javax.swing.JLabel innerLabel;
    private javax.swing.JTextField innerTextField;
    private javax.swing.JLabel lpLabel;
    private javax.swing.JTextField lpTextField;
    private javax.swing.JLabel lpstopLabel;
    private javax.swing.JTextField lpstopTextField;
    private javax.swing.JLabel lpvalidLabel;
    private javax.swing.JTextField lpvalidTextField;
    private javax.swing.JPanel mainSimplePrimePanel;
    private javax.swing.JScrollPane mainSimplePrimeScrollPanel;
    private javax.swing.JLabel maxitsLabel;
    private javax.swing.JTextField maxitsTextField;
    private javax.swing.JPanel microscopeDetailsPanel;
    private javax.swing.JLabel mskLabel;
    private javax.swing.JTextField mskTextField;
    private javax.swing.JRadioButton mulCTFPrimeRadioButton;
    private javax.swing.JLabel mwLabel;
    private javax.swing.JTextField mwTextField;
    private javax.swing.JLabel nnnLabel;
    private javax.swing.JTextField nnnTextField;
    private javax.swing.JRadioButton noCTFPrimeRadioButton;
    private javax.swing.JRadioButton noDiversifyRadioButton;
    private javax.swing.JRadioButton noDynlpRadioButton;
    private javax.swing.JRadioButton noEoRadioButton;
    private javax.swing.JRadioButton noNoiseRadioButton;
    private javax.swing.JRadioButton noNorecRadioButton;
    private javax.swing.JRadioButton noRefinementRadioButton;
    private javax.swing.JLabel noiseLabel;
    private javax.swing.ButtonGroup noisebuttonGroup;
    private javax.swing.JLabel norecLabel;
    private javax.swing.ButtonGroup norecbuttonGroup;
    private javax.swing.JLabel nspaceLabel;
    private javax.swing.JTextField nspaceTextField;
    private javax.swing.JLabel nstatesLabel;
    private javax.swing.JTextField nstatesTextField;
    private javax.swing.JLabel nvoxLabel;
    private javax.swing.JTextField nvoxTextField;
    private javax.swing.JRadioButton oSymmetryRadioButton;
    private javax.swing.JTextField oriTabTextField;
    private javax.swing.JLabel oritabLabel;
    private javax.swing.JPanel otherParametersPanel;
    private javax.swing.JScrollPane otherParametersScrollPanel;
    private javax.swing.JLabel outfileLabel;
    private javax.swing.JTextField outfileTextField;
    private javax.swing.JPanel refinePanel;
    private javax.swing.ButtonGroup refinementbuttonGroup;
    private javax.swing.JRadioButton shiftRefinementRadioButton;
    private javax.swing.JLabel smpdLabel;
    private javax.swing.JTextField smpdTextField;
    private javax.swing.JPanel stackDetailsPanel;
    private javax.swing.JTextField stackNameTextField;
    private javax.swing.JLabel stakNameLabel;
    private javax.swing.JLabel startitLabel;
    private javax.swing.JTextField startitTextField;
    private javax.swing.JPanel symmetryGroupPanel;
    private javax.swing.ButtonGroup symmetryGroupbuttonGroup;
    private javax.swing.JRadioButton tSymmetryRadioButton;
    private javax.swing.JLabel timePerImageLabel;
    private javax.swing.JTextField timePerImageTextField;
    private javax.swing.JLabel trsLabel;
    private javax.swing.JTextField trsTextField;
    private javax.swing.JLabel trsstepLabel;
    private javax.swing.JTextField trsstepTextField;
    private javax.swing.JLabel vol1NameLabel;
    private javax.swing.JTextField vol1NameTextField;
    private javax.swing.JLabel vol2NameLabel;
    private javax.swing.JTextField vol2NameTextField;
    private javax.swing.JLabel widthLabel;
    private javax.swing.JTextField widthTextField;
    private javax.swing.JRadioButton yesCTFPrimeRadioButton;
    private javax.swing.JRadioButton yesDiversifyRadioButton;
    private javax.swing.JRadioButton yesDynlpRadioButton;
    private javax.swing.JRadioButton yesEoRadioButton;
    private javax.swing.JRadioButton yesNoiseRadioButton;
    private javax.swing.JRadioButton yesNorecRadioButton;
    // End of variables declaration//GEN-END:variables

}
