/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * SimpleStackopsJPanel.java
 *
 * Created on 26/11/2015, 9:49:41 PM
 */

package edu.kit.gridbeans.simplestackops.plugin;

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

import edu.kit.gridbeans.simplestackops.CSymmetrySimpleStackopsCalculationType;
import edu.kit.gridbeans.simplestackops.DSymmetrySimpleStackopsCalculationType;
import edu.kit.gridbeans.simplestackops.SimpleStackopsGridBean;

/**
 *
 * @author Frederic Bonnet
 * @Date 26th November 2015
 * 
 */
public class SimpleStackopsJPanel extends JPanel {
    /**
     * Add UID
     */
    private static final long serialVersionUID = -3217081134128067037L;
    private int fracslidernumber;
    private int fracZeroslidernumber;
    private int threshold4Binslidernumber;
    private GridBeanPanel parentPanel;
    /** Creates new form SimpleStackopsJPanel */
    public SimpleStackopsJPanel(GridBeanPanel parentpanel) throws DataSetException {
        this.parentPanel = parentpanel;
        initComponents();
        initCustomComponents();
        bindComponents();       
    }

    private void initCustomComponents() {
        //CTF button group
        ctfbuttonGroup.add(yesCTFRadioButton);
        ctfbuttonGroup.add(noCTFRadioButton);
        ctfbuttonGroup.add(flipCTFRadioButton);
        ctfbuttonGroup.add(mulCTFRadioButton);
        //CTF sq button group
        ctfsqbuttonGroup.add(yesCtfSqRadioButton);
        ctfsqbuttonGroup.add(noCtfSqRadioButton);
        //ACF button group
        acfbuttonGroup.add(yesAcfRadioButton);
        acfbuttonGroup.add(noAcfRadioButton);
        //masscen button group
        masscenbuttonGroup.add(yesMasscenRadioButton);
        masscenbuttonGroup.add(noMasscenRadioButton);
        //phrand button group
        phrandbuttonGroup.add(yesPhrandRadioButton);
        phrandbuttonGroup.add(noPhrandRadioButton);
        //Vis button group
        visbuttonGroup.add(yesVisRadioButton);
        visbuttonGroup.add(noVisRadioButton);
        //Shift alignment parameters
        shiftAlignmentbuttonGroup.add(yesShalgnRadioButton);
        shiftAlignmentbuttonGroup.add(noShalgnRadioButton);
        //Ro Alignment parameters button group
        roalgnbuttonGroup.add(yesRoalgnRadioButton);
        roalgnbuttonGroup.add(noRoalgnRadioButton);
        //Binarize button group
        binarizebuttonGroup.add(yesBinarizeRadioButton);
        binarizebuttonGroup.add(noBinarizeRadioButton);
        //Hfun button group
        hfunbuttonGroup.add(sigmHfunRadioButton);
        hfunbuttonGroup.add(tanhHfunRadioButton);
        hfunbuttonGroup.add(linHfunRadioButton);
        //Normalise button group
        normalisebuttonGroup.add(yesNormaliseRadioButton);
        normalisebuttonGroup.add(noNormaliseRadioButton);
        //Normalize noise button group
        normNoisebuttonGroup.add(yesNormNoiseRadioButton);
        normNoisebuttonGroup.add(noNormNoiseRadioButton);
        //Avg button group
        avgbuttonGroup.add(yesAvgRadioButton);
        avgbuttonGroup.add(noAvgRadioButton);
        //Rantify button group
        ranktifybuttonGroup.add(yesRanktifyRadioButton);
        ranktifybuttonGroup.add(noRanktifyRadioButton);
        //Stats button group
        statsbuttonGroup.add(yesStatsRadioButton1);
        statsbuttonGroup.add(noStatsRadioButton);
        //Compare button group
        comparebuttonGroup.add(yesCompareRadioButton);
        comparebuttonGroup.add(noCompareRadioButton);
        //Merge button group
        mergebuttonGroup.add(yesMergeRadioButton);
        mergebuttonGroup.add(noMergeRadioButton);
        //ft2img button group
        ft2imgbuttonGroup.add(yesFt2imgRadioButton);
        ft2imgbuttonGroup.add(noFt2imgRadioButton);
        
        //Symmetry button group
        symmetryGroupbuttonGroup.add(cSymmetryRadioButton);
        symmetryGroupbuttonGroup.add(dSymmetryRadioButton);
        symmetryGroupbuttonGroup.add(tSymmetryRadioButton);
        symmetryGroupbuttonGroup.add(iSymmetryRadioButton);
        symmetryGroupbuttonGroup.add(oSymmetryRadioButton);
       
    }

    private void bindComponents() throws DataSetException {
        //STK
        parentPanel.linkTextField(SimpleStackopsGridBean.STK, stackNameTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.STK, StringValueTranslator.getInstance());
        //STK2
        parentPanel.linkTextField(SimpleStackopsGridBean.STK2, stack2NameTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.STK2, StringValueTranslator.getInstance());
        //Orientation table
        parentPanel.linkTextField(SimpleStackopsGridBean.ORITAB, oriTabTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.ORITAB, StringValueTranslator.getInstance());
        //OUTSTK
        parentPanel.linkTextField(SimpleStackopsGridBean.OUTSTK, outstackTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.OUTSTK, StringValueTranslator.getInstance());
        //Defocus values
        parentPanel.linkTextField(SimpleStackopsGridBean.DEFTAB, defTabTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.DEFTAB, StringValueTranslator.getInstance());
        //Output file for the parallel output
        parentPanel.linkTextField(SimpleStackopsGridBean.OUTFILE, outputParameterFileTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.OUTFILE, StringValueTranslator.getInstance());
        //Filetab values
        parentPanel.linkTextField(SimpleStackopsGridBean.FILETAB, fileTabTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.FILETAB, StringValueTranslator.getInstance());
        //Sampling distance
        parentPanel.linkTextField(SimpleStackopsGridBean.SMPD, smpdTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.SMPD, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.SMPD, new DoubleValueValidator());
        //Acceleration volatage in KeV
        parentPanel.linkTextField(SimpleStackopsGridBean.KV, accelerationVoltageTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.KV, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.KV, new DoubleValueValidator());
        //Shperical abberation
        parentPanel.linkTextField(SimpleStackopsGridBean.CS, csTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.CS, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.CS, new DoubleValueValidator());
        //frac amp contrast{0.07}
        parentPanel.linkTextField(SimpleStackopsGridBean.FRACA, fracaTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.FRACA, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.FRACA, new DoubleValueValidator());
        //High-pass limit in A
        parentPanel.linkTextField(SimpleStackopsGridBean.HP, hpTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.HP, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.HP, new DoubleValueValidator());
        //Low-pass limit {20}
        parentPanel.linkTextField(SimpleStackopsGridBean.LP, lpTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.LP, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.LP, new DoubleValueValidator());
        //--fromp=$FROMP \           #[fromp=<start ptcl>]
        parentPanel.linkTextField(SimpleStackopsGridBean.FROMP, frompTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.FROMP, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.FROMP, new IntegerValueValidator());
        //--mul=$MUL \               #[mul=<shift multiplication factor{1}>]
        parentPanel.linkTextField(SimpleStackopsGridBean.MUL, mulTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.MUL, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.MUL, new DoubleValueValidator());
        //--trs=$TRS \                       #[trs=<origin shift(in pixels){0}>]
        parentPanel.linkTextField(SimpleStackopsGridBean.TRS, trsTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.TRS, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.TRS, new IntegerValueValidator());
        //--snr=$SNR \               #[snr=<signal2noise ratio>]
        parentPanel.linkTextField(SimpleStackopsGridBean.SNR, snrTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.SNR, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.SNR, new DoubleValueValidator());
        //--state=$STATE \           #[state=<state to extract>]
        parentPanel.linkTextField(SimpleStackopsGridBean.STATE, stateTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.STATE, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.STATE, new IntegerValueValidator());
        //--msk=$MSK \                      #msk=<mask radius(in pixels)> smpd=<sampling distance(in A)>
        parentPanel.linkTextField(SimpleStackopsGridBean.MSK, mskTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.MSK, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.MSK, new DoubleValueValidator());
        //--nran=$NRAN \             #[nran=<number of random images to select>]
        parentPanel.linkTextField(SimpleStackopsGridBean.NRAN, nranTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.NRAN, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.NRAN, new IntegerValueValidator());
        //--newbox=$NEWBOX \         #[newbox=<scaled box>]
        parentPanel.linkTextField(SimpleStackopsGridBean.NEWBOX, newBoxTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.NEWBOX, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.NEWBOX, new IntegerValueValidator());
        //--scale=$SCALE \           #[scale=<scale factor{1}>]
        parentPanel.linkTextField(SimpleStackopsGridBean.SCALE, scaleTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.SCALE, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.SCALE, new DoubleValueValidator());
        //--frameavg=$FRAMEAVG \     #[frameavg=<nr of frames to average{0}>] 
        parentPanel.linkTextField(SimpleStackopsGridBean.FRAMEAVG, frameAvgTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.FRAMEAVG, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.FRAMEAVG, new IntegerValueValidator());
        //--clip=$CLIP \             #[clip=<clipped box size{box}>]
        parentPanel.linkTextField(SimpleStackopsGridBean.CLIP, clipTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.CLIP, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.CLIP, new IntegerValueValidator());
        //--box=$BOX \               #[box=<image size(in pixels)>]
        parentPanel.linkTextField(SimpleStackopsGridBean.BOX, boxTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.BOX, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.BOX, new IntegerValueValidator());
        //--inner=$INNER \                   #[inner=<inner mask radius(in pixels)>]
        parentPanel.linkTextField(SimpleStackopsGridBean.INNER, innerTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.INNER, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.INNER, new DoubleValueValidator());
        //--width=$WIDTH \                   #[width=<pixels falloff inner mask{10}>]
        parentPanel.linkTextField(SimpleStackopsGridBean.WIDTH, widthTextField);
        parentPanel.setValueTranslator(SimpleStackopsGridBean.WIDTH, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleStackopsGridBean.WIDTH, new DoubleValueValidator());
        //link to the JSlider for the frac
        class FracSliderNumberListener implements ChangeListener {
            public void stateChanged(ChangeEvent event) {
                setFracSliderNumber();
            }
        }
        ChangeListener verbosityperfslidernumberlistener = new FracSliderNumberListener();
        fracSlider.addChangeListener(verbosityperfslidernumberlistener);
        //link to the JSlider for the frac Zero
        class FracZeroSliderNumberListener implements ChangeListener {
            public void stateChanged(ChangeEvent event) {
                setFracZeroSliderNumber();
            }
        }
        ChangeListener fracValidslidernumberlistener = new FracZeroSliderNumberListener();
        fracZeroSlider.addChangeListener(fracValidslidernumberlistener);
        //link to the JSlider for the tres ie tres4bin slider
        class Threshold4BinSliderNumberListener implements ChangeListener {
            public void stateChanged(ChangeEvent event) {
                setThreshold4BinSliderNUmer();
            }
        }
        ChangeListener threshold4Binslidernumberlistener = new Threshold4BinSliderNumberListener();
        threshold4BinSlider.addChangeListener(threshold4Binslidernumberlistener);
    }
    /**
     * method for the JSlider assignment method for tres ie tres4bin slider
     */
    private void setThreshold4BinSliderNUmer() {
        threshold4Binslidernumber = threshold4BinSlider.getValue();
        
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue(Integer.toString(threshold4Binslidernumber));
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.TRES, value);
    }
    /**
     * method for the JSlider assignment method for frac
     */
    public void setFracSliderNumber() {
        fracslidernumber = fracSlider.getValue();

        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue(Integer.toString(fracslidernumber));
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.FRAC, value);
    }
    /**
     * method for the JSlider assignment method for frac Valid
     */
    public void setFracZeroSliderNumber() {
        fracZeroslidernumber = fracZeroSlider.getValue();
        
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue(Integer.toString(fracZeroslidernumber));
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.FRACZERO, value);
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        shiftAlignmentbuttonGroup = new javax.swing.ButtonGroup();
        roalgnbuttonGroup = new javax.swing.ButtonGroup();
        acfbuttonGroup = new javax.swing.ButtonGroup();
        masscenbuttonGroup = new javax.swing.ButtonGroup();
        phrandbuttonGroup = new javax.swing.ButtonGroup();
        visbuttonGroup = new javax.swing.ButtonGroup();
        ctfsqbuttonGroup = new javax.swing.ButtonGroup();
        binarizebuttonGroup = new javax.swing.ButtonGroup();
        hfunbuttonGroup = new javax.swing.ButtonGroup();
        normalisebuttonGroup = new javax.swing.ButtonGroup();
        normNoisebuttonGroup = new javax.swing.ButtonGroup();
        avgbuttonGroup = new javax.swing.ButtonGroup();
        ranktifybuttonGroup = new javax.swing.ButtonGroup();
        statsbuttonGroup = new javax.swing.ButtonGroup();
        comparebuttonGroup = new javax.swing.ButtonGroup();
        negbuttonGroup = new javax.swing.ButtonGroup();
        mergebuttonGroup = new javax.swing.ButtonGroup();
        ft2imgbuttonGroup = new javax.swing.ButtonGroup();
        ctfbuttonGroup = new javax.swing.ButtonGroup();
        symmetryGroupbuttonGroup = new javax.swing.ButtonGroup();
        mainSimpleStackopsScrollPanel = new javax.swing.JScrollPane();
        mainSimpleStackopsPanel = new javax.swing.JPanel();
        microscopeDetailsPanel = new javax.swing.JPanel();
        smpdLabel = new javax.swing.JLabel();
        smpdTextField = new javax.swing.JTextField();
        accelerationVoltageLabel = new javax.swing.JLabel();
        accelerationVoltageTextField = new javax.swing.JTextField();
        csLabel = new javax.swing.JLabel();
        csTextField = new javax.swing.JTextField();
        fracaLabel = new javax.swing.JLabel();
        fracaTextField = new javax.swing.JTextField();
        nptclsLabel = new javax.swing.JLabel();
        nptclsTextField = new javax.swing.JTextField();
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
        stak2NameLabel = new javax.swing.JLabel();
        stack2NameTextField = new javax.swing.JTextField();
        browseStk2Button = new javax.swing.JButton();
        outstackLabel = new javax.swing.JLabel();
        outstackTextField = new javax.swing.JTextField();
        browseOutStackButton = new javax.swing.JButton();
        optionalPanel = new javax.swing.JPanel();
        shalgnLabel = new javax.swing.JLabel();
        yesShalgnRadioButton = new javax.swing.JRadioButton();
        noShalgnRadioButton = new javax.swing.JRadioButton();
        roalgnLabel = new javax.swing.JLabel();
        yesRoalgnRadioButton = new javax.swing.JRadioButton();
        noRoalgnRadioButton = new javax.swing.JRadioButton();
        acfOtherPanel = new javax.swing.JPanel();
        acfLabel = new javax.swing.JLabel();
        yesAcfRadioButton = new javax.swing.JRadioButton();
        noAcfRadioButton = new javax.swing.JRadioButton();
        masscenLabel = new javax.swing.JLabel();
        yesMasscenRadioButton = new javax.swing.JRadioButton();
        noMasscenRadioButton = new javax.swing.JRadioButton();
        phrandLabel = new javax.swing.JLabel();
        yesPhrandRadioButton = new javax.swing.JRadioButton();
        noPhrandRadioButton = new javax.swing.JRadioButton();
        visLabel = new javax.swing.JLabel();
        yesVisRadioButton = new javax.swing.JRadioButton();
        noVisRadioButton = new javax.swing.JRadioButton();
        optionalScrollPanel = new javax.swing.JScrollPane();
        optionalParametersPanel = new javax.swing.JPanel();
        hpLabel = new javax.swing.JLabel();
        hpTextField = new javax.swing.JTextField();
        frompLabel = new javax.swing.JLabel();
        frompTextField = new javax.swing.JTextField();
        mulLabel = new javax.swing.JLabel();
        mulTextField = new javax.swing.JTextField();
        trsLabel = new javax.swing.JLabel();
        trsTextField = new javax.swing.JTextField();
        lpLabel = new javax.swing.JLabel();
        lpTextField = new javax.swing.JTextField();
        snrLabel = new javax.swing.JLabel();
        snrTextField = new javax.swing.JTextField();
        stateLabel = new javax.swing.JLabel();
        stateTextField = new javax.swing.JTextField();
        mskLabel = new javax.swing.JLabel();
        mskTextField = new javax.swing.JTextField();
        binarizeLabel = new javax.swing.JLabel();
        yesBinarizeRadioButton = new javax.swing.JRadioButton();
        noBinarizeRadioButton = new javax.swing.JRadioButton();
        topLabel = new javax.swing.JLabel();
        topTextField = new javax.swing.JTextField();
        nranLabel = new javax.swing.JLabel();
        nranTextField = new javax.swing.JTextField();
        newBoxLabel = new javax.swing.JLabel();
        newBoxTextField = new javax.swing.JTextField();
        scaleLabel = new javax.swing.JLabel();
        hfunLabel = new javax.swing.JLabel();
        sigmHfunRadioButton = new javax.swing.JRadioButton();
        tanhHfunRadioButton = new javax.swing.JRadioButton();
        linHfunRadioButton = new javax.swing.JRadioButton();
        scaleTextField = new javax.swing.JTextField();
        normaliseLabel = new javax.swing.JLabel();
        yesNormaliseRadioButton = new javax.swing.JRadioButton();
        noNormaliseRadioButton = new javax.swing.JRadioButton();
        norm_noiseLabel = new javax.swing.JLabel();
        yesNormNoiseRadioButton = new javax.swing.JRadioButton();
        noNormNoiseRadioButton = new javax.swing.JRadioButton();
        avgLabel = new javax.swing.JLabel();
        yesAvgRadioButton = new javax.swing.JRadioButton();
        noAvgRadioButton = new javax.swing.JRadioButton();
        ranktifyLabel = new javax.swing.JLabel();
        yesRanktifyRadioButton = new javax.swing.JRadioButton();
        noRanktifyRadioButton = new javax.swing.JRadioButton();
        fileTabLabel = new javax.swing.JLabel();
        fileTabTextField = new javax.swing.JTextField();
        browseFileTabButton = new javax.swing.JButton();
        statsLabel = new javax.swing.JLabel();
        yesStatsRadioButton1 = new javax.swing.JRadioButton();
        noStatsRadioButton = new javax.swing.JRadioButton();
        compareLabel = new javax.swing.JLabel();
        yesCompareRadioButton = new javax.swing.JRadioButton();
        noCompareRadioButton = new javax.swing.JRadioButton();
        negLabel = new javax.swing.JLabel();
        yesNegRadioButton = new javax.swing.JRadioButton();
        noNegRadioButton = new javax.swing.JRadioButton();
        mergeLabel = new javax.swing.JLabel();
        yesMergeRadioButton = new javax.swing.JRadioButton();
        noMergeRadioButton = new javax.swing.JRadioButton();
        ft2imgLabel = new javax.swing.JLabel();
        yesFt2imgRadioButton = new javax.swing.JRadioButton();
        noFt2imgRadioButton = new javax.swing.JRadioButton();
        frameAvgLabel = new javax.swing.JLabel();
        frameAvgTextField = new javax.swing.JTextField();
        clipLabel = new javax.swing.JLabel();
        clipTextField = new javax.swing.JTextField();
        boxLabel = new javax.swing.JLabel();
        boxTextField = new javax.swing.JTextField();
        innerLabel = new javax.swing.JLabel();
        innerTextField = new javax.swing.JTextField();
        widthLabel = new javax.swing.JLabel();
        widthTextField = new javax.swing.JTextField();
        outputParameterFileLabel = new javax.swing.JLabel();
        outputParameterFileTextField = new javax.swing.JTextField();
        browseOutputParameterFileButton = new javax.swing.JButton();
        fracZeroSlider = new javax.swing.JSlider();
        threshold4BinSlider = new javax.swing.JSlider();
        symmetryCTFScrollPanel = new javax.swing.JScrollPane();
        symmetryCTFPanel = new javax.swing.JPanel();
        symmetryGroupPanel = new javax.swing.JPanel();
        cSymmetryGroupComboBox = new javax.swing.JComboBox();
        dSymmetryGroupComboBox = new javax.swing.JComboBox();
        tSymmetryRadioButton = new javax.swing.JRadioButton();
        oSymmetryRadioButton = new javax.swing.JRadioButton();
        iSymmetryRadioButton = new javax.swing.JRadioButton();
        cSymmetryRadioButton = new javax.swing.JRadioButton();
        dSymmetryRadioButton = new javax.swing.JRadioButton();
        ctfPanel = new javax.swing.JPanel();
        yesCTFRadioButton = new javax.swing.JRadioButton();
        noCTFRadioButton = new javax.swing.JRadioButton();
        flipCTFRadioButton = new javax.swing.JRadioButton();
        mulCTFRadioButton = new javax.swing.JRadioButton();
        ctfsqLabel = new javax.swing.JLabel();
        yesCtfSqRadioButton = new javax.swing.JRadioButton();
        noCtfSqRadioButton = new javax.swing.JRadioButton();
        fracSlider = new javax.swing.JSlider();

        setLayout(new javax.swing.BoxLayout(this, javax.swing.BoxLayout.LINE_AXIS));

        microscopeDetailsPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("MIcroscope details..."));

        smpdLabel.setText("Sampling distance smpd:");

        accelerationVoltageLabel.setText("acceleration voltage(in kV) kv:");

        csLabel.setText("Spherical Aberration cs:");

        fracaLabel.setText("Ampl. Contrast fraca:");

        nptclsLabel.setText("n particles nptcls:");

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
                            .addComponent(nptclsLabel))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(nptclsTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 184, Short.MAX_VALUE)
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
                    .addComponent(nptclsTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(nptclsLabel))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        stackDetailsPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Stack and input details..."));

        stakNameLabel.setText("Stk name:");

        browseStkButton.setText("Browse");

        oritabLabel.setText("Ori Table:");

        browseOriTabButton.setText("Browse");

        deftabLabel.setText("Def Table:");

        browseDefTabButton.setText("Browse");

        stak2NameLabel.setText("Stk2 name:");

        browseStk2Button.setText("Browse");

        outstackLabel.setText("Outstack:");

        browseOutStackButton.setText("Browse");

        javax.swing.GroupLayout stackDetailsPanelLayout = new javax.swing.GroupLayout(stackDetailsPanel);
        stackDetailsPanel.setLayout(stackDetailsPanelLayout);
        stackDetailsPanelLayout.setHorizontalGroup(
            stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, stackDetailsPanelLayout.createSequentialGroup()
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, stackDetailsPanelLayout.createSequentialGroup()
                        .addComponent(deftabLabel)
                        .addGap(18, 18, 18)
                        .addComponent(defTabTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 145, Short.MAX_VALUE))
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, stackDetailsPanelLayout.createSequentialGroup()
                        .addComponent(outstackLabel)
                        .addGap(18, 18, 18)
                        .addComponent(outstackTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 145, Short.MAX_VALUE))
                    .addGroup(stackDetailsPanelLayout.createSequentialGroup()
                        .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(stakNameLabel)
                            .addComponent(stak2NameLabel)
                            .addComponent(oritabLabel))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(oriTabTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 145, Short.MAX_VALUE)
                            .addComponent(stack2NameTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 145, Short.MAX_VALUE)
                            .addComponent(stackNameTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 145, Short.MAX_VALUE))))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(browseOriTabButton)
                    .addComponent(browseStkButton)
                    .addComponent(browseStk2Button)
                    .addComponent(browseOutStackButton)
                    .addComponent(browseDefTabButton))
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
                    .addComponent(stack2NameTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(stak2NameLabel)
                    .addComponent(browseStk2Button))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(oritabLabel)
                    .addComponent(browseOriTabButton)
                    .addComponent(oriTabTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(outstackLabel)
                    .addComponent(outstackTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(browseOutStackButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(deftabLabel)
                    .addComponent(defTabTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(browseDefTabButton)))
        );

        optionalPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Alignment parameter..."));

        shalgnLabel.setText("Shift Algn shalgn:");

        yesShalgnRadioButton.setText("Yes");
        yesShalgnRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesShalgnRadioButtonActionPerformed(evt);
            }
        });

        noShalgnRadioButton.setText("No");
        noShalgnRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noShalgnRadioButtonActionPerformed(evt);
            }
        });

        roalgnLabel.setText("Ro Algn roalgn:");

        yesRoalgnRadioButton.setText("Yes");
        yesRoalgnRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesRoalgnRadioButtonActionPerformed(evt);
            }
        });

        noRoalgnRadioButton.setText("No");
        noRoalgnRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noRoalgnRadioButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout optionalPanelLayout = new javax.swing.GroupLayout(optionalPanel);
        optionalPanel.setLayout(optionalPanelLayout);
        optionalPanelLayout.setHorizontalGroup(
            optionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(optionalPanelLayout.createSequentialGroup()
                .addGroup(optionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(optionalPanelLayout.createSequentialGroup()
                        .addComponent(shalgnLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(yesShalgnRadioButton))
                    .addGroup(optionalPanelLayout.createSequentialGroup()
                        .addComponent(roalgnLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(yesRoalgnRadioButton)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(noShalgnRadioButton)
                    .addComponent(noRoalgnRadioButton))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        optionalPanelLayout.setVerticalGroup(
            optionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(optionalPanelLayout.createSequentialGroup()
                .addGroup(optionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(shalgnLabel)
                    .addComponent(yesShalgnRadioButton)
                    .addComponent(noShalgnRadioButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(roalgnLabel)
                    .addComponent(yesRoalgnRadioButton)
                    .addComponent(noRoalgnRadioButton))
                .addContainerGap(21, Short.MAX_VALUE))
        );

        acfOtherPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Acf masscen phrand..."));

        acfLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        acfLabel.setText("Acf:");

        yesAcfRadioButton.setText("Yes");
        yesAcfRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesAcfRadioButtonActionPerformed(evt);
            }
        });

        noAcfRadioButton.setText("No");
        noAcfRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noAcfRadioButtonActionPerformed(evt);
            }
        });

        masscenLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        masscenLabel.setText("Masscen:");

        yesMasscenRadioButton.setText("Yes");
        yesMasscenRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesMasscenRadioButtonActionPerformed(evt);
            }
        });

        noMasscenRadioButton.setText("No");
        noMasscenRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noMasscenRadioButtonActionPerformed(evt);
            }
        });

        phrandLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        phrandLabel.setText("Phrand:");

        yesPhrandRadioButton.setText("Yes");
        yesPhrandRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesPhrandRadioButtonActionPerformed(evt);
            }
        });

        noPhrandRadioButton.setText("No");
        noPhrandRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noPhrandRadioButtonActionPerformed(evt);
            }
        });

        visLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        visLabel.setText("Vis:");

        yesVisRadioButton.setText("Yes");
        yesVisRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesVisRadioButtonActionPerformed(evt);
            }
        });

        noVisRadioButton.setText("No");
        noVisRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noVisRadioButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout acfOtherPanelLayout = new javax.swing.GroupLayout(acfOtherPanel);
        acfOtherPanel.setLayout(acfOtherPanelLayout);
        acfOtherPanelLayout.setHorizontalGroup(
            acfOtherPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(acfOtherPanelLayout.createSequentialGroup()
                .addGroup(acfOtherPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(acfOtherPanelLayout.createSequentialGroup()
                        .addComponent(acfLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 33, Short.MAX_VALUE)
                        .addComponent(yesAcfRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noAcfRadioButton))
                    .addGroup(acfOtherPanelLayout.createSequentialGroup()
                        .addComponent(phrandLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(yesPhrandRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noPhrandRadioButton)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(acfOtherPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(acfOtherPanelLayout.createSequentialGroup()
                        .addComponent(visLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(yesVisRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noVisRadioButton))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, acfOtherPanelLayout.createSequentialGroup()
                        .addComponent(masscenLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(yesMasscenRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noMasscenRadioButton)))
                .addContainerGap())
        );
        acfOtherPanelLayout.setVerticalGroup(
            acfOtherPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(acfOtherPanelLayout.createSequentialGroup()
                .addGroup(acfOtherPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(acfLabel)
                    .addComponent(masscenLabel)
                    .addComponent(yesMasscenRadioButton)
                    .addComponent(noMasscenRadioButton)
                    .addComponent(yesAcfRadioButton)
                    .addComponent(noAcfRadioButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(acfOtherPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(phrandLabel)
                    .addComponent(yesPhrandRadioButton)
                    .addComponent(noPhrandRadioButton)
                    .addComponent(visLabel)
                    .addComponent(yesVisRadioButton)
                    .addComponent(noVisRadioButton))
                .addContainerGap(21, Short.MAX_VALUE))
        );

        optionalParametersPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Optional parameters..."));

        hpLabel.setText("High pass limit hp:");

        frompLabel.setText("Start ptcl fromp:");

        mulLabel.setText("Shifting Mul mul:");

        trsLabel.setText("trs<origin>:");

        lpLabel.setText("Low pass limit lp:");

        snrLabel.setText("SNR:");

        stateLabel.setText("State to extract:");

        mskLabel.setText("Mask msk:");

        binarizeLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        binarizeLabel.setText("Binarize:");

        yesBinarizeRadioButton.setText("Yes");
        yesBinarizeRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesBinarizeRadioButtonActionPerformed(evt);
            }
        });

        noBinarizeRadioButton.setText("No");
        noBinarizeRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noBinarizeRadioButtonActionPerformed(evt);
            }
        });

        topLabel.setFont(new java.awt.Font("Dialog", 0, 12));
        topLabel.setText("Top:");

        nranLabel.setFont(new java.awt.Font("Dialog", 0, 12));
        nranLabel.setText("nran:");

        newBoxLabel.setText("NewBox:");

        scaleLabel.setText("Scale:");

        hfunLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        hfunLabel.setText("Hfun:");

        sigmHfunRadioButton.setText("sigm");
        sigmHfunRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                sigmHfunRadioButtonActionPerformed(evt);
            }
        });

        tanhHfunRadioButton.setText("tanh");
        tanhHfunRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                tanhHfunRadioButtonActionPerformed(evt);
            }
        });

        linHfunRadioButton.setText("lin");
        linHfunRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                linHfunRadioButtonActionPerformed(evt);
            }
        });

        normaliseLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        normaliseLabel.setText("Normalise:");

        yesNormaliseRadioButton.setText("Yes");
        yesNormaliseRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesNormaliseRadioButtonActionPerformed(evt);
            }
        });

        noNormaliseRadioButton.setText("No");
        noNormaliseRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noNormaliseRadioButtonActionPerformed(evt);
            }
        });

        norm_noiseLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        norm_noiseLabel.setText("Normalise Noise:");

        yesNormNoiseRadioButton.setText("Yes");
        yesNormNoiseRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesNormNoiseRadioButtonActionPerformed(evt);
            }
        });

        noNormNoiseRadioButton.setText("No");
        noNormNoiseRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noNormNoiseRadioButtonActionPerformed(evt);
            }
        });

        avgLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        avgLabel.setText("Avg:");

        yesAvgRadioButton.setText("Yes");
        yesAvgRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesAvgRadioButtonActionPerformed(evt);
            }
        });

        noAvgRadioButton.setText("No");
        noAvgRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noAvgRadioButtonActionPerformed(evt);
            }
        });

        ranktifyLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        ranktifyLabel.setText("Ranktify:");

        yesRanktifyRadioButton.setText("Yes");
        yesRanktifyRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesRanktifyRadioButtonActionPerformed(evt);
            }
        });

        noRanktifyRadioButton.setText("No");
        noRanktifyRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noRanktifyRadioButtonActionPerformed(evt);
            }
        });

        fileTabLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        fileTabLabel.setText("FileTab:");

        browseFileTabButton.setText("Browse");

        statsLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        statsLabel.setText("Stats:");

        yesStatsRadioButton1.setText("Yes");
        yesStatsRadioButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesStatsRadioButton1ActionPerformed(evt);
            }
        });

        noStatsRadioButton.setText("No");
        noStatsRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noStatsRadioButtonActionPerformed(evt);
            }
        });

        compareLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        compareLabel.setText("Compare:");

        yesCompareRadioButton.setText("Yes");
        yesCompareRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesCompareRadioButtonActionPerformed(evt);
            }
        });

        noCompareRadioButton.setText("No");
        noCompareRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noCompareRadioButtonActionPerformed(evt);
            }
        });

        negLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        negLabel.setText("Neg:");

        yesNegRadioButton.setText("Yes");
        yesNegRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesNegRadioButtonActionPerformed(evt);
            }
        });

        noNegRadioButton.setText("No");
        noNegRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noNegRadioButtonActionPerformed(evt);
            }
        });

        mergeLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        mergeLabel.setText("Merge:");

        yesMergeRadioButton.setText("Yes");
        yesMergeRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesMergeRadioButtonActionPerformed(evt);
            }
        });

        noMergeRadioButton.setText("No");
        noMergeRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noMergeRadioButtonActionPerformed(evt);
            }
        });

        ft2imgLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        ft2imgLabel.setText("ft2img:");

        yesFt2imgRadioButton.setText("Yes");
        yesFt2imgRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesFt2imgRadioButtonActionPerformed(evt);
            }
        });

        noFt2imgRadioButton.setText("No");
        noFt2imgRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noFt2imgRadioButtonActionPerformed(evt);
            }
        });

        frameAvgLabel.setText("Frames to avg:");

        clipLabel.setText("Clip:");

        boxLabel.setText("Box:");

        innerLabel.setText("Inner:");

        widthLabel.setText("Width:");

        outputParameterFileLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        outputParameterFileLabel.setText("Output parameter file:");

        browseOutputParameterFileButton.setText("Browse");

        fracZeroSlider.setFont(new java.awt.Font("Dialog", 0, 8));
        fracZeroSlider.setMajorTickSpacing(10);
        fracZeroSlider.setMinorTickSpacing(2);
        fracZeroSlider.setPaintLabels(true);
        fracZeroSlider.setPaintTicks(true);
        fracZeroSlider.setSnapToTicks(true);
        fracZeroSlider.setToolTipText("Fraction of zeros to include ");
        fracZeroSlider.setValue(80);
        fracZeroSlider.setBorder(javax.swing.BorderFactory.createTitledBorder("Fraction of zeros to include..."));
        fracZeroSlider.setValueIsAdjusting(true);

        threshold4BinSlider.setFont(new java.awt.Font("Dialog", 0, 8));
        threshold4BinSlider.setMajorTickSpacing(10);
        threshold4BinSlider.setMinorTickSpacing(2);
        threshold4BinSlider.setPaintLabels(true);
        threshold4BinSlider.setPaintTicks(true);
        threshold4BinSlider.setSnapToTicks(true);
        threshold4BinSlider.setToolTipText("Threshold for bin...");
        threshold4BinSlider.setValue(60);
        threshold4BinSlider.setBorder(javax.swing.BorderFactory.createTitledBorder("Threshold for bins..."));
        threshold4BinSlider.setValueIsAdjusting(true);

        javax.swing.GroupLayout optionalParametersPanelLayout = new javax.swing.GroupLayout(optionalParametersPanel);
        optionalParametersPanel.setLayout(optionalParametersPanelLayout);
        optionalParametersPanelLayout.setHorizontalGroup(
            optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                        .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                                .addComponent(binarizeLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(yesBinarizeRadioButton)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(noBinarizeRadioButton))
                            .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                .addGroup(javax.swing.GroupLayout.Alignment.LEADING, optionalParametersPanelLayout.createSequentialGroup()
                                    .addComponent(lpLabel)
                                    .addGap(18, 18, 18)
                                    .addComponent(lpTextField))
                                .addGroup(javax.swing.GroupLayout.Alignment.LEADING, optionalParametersPanelLayout.createSequentialGroup()
                                    .addComponent(hpLabel)
                                    .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                    .addComponent(hpTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE))))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(frompLabel)
                            .addComponent(snrLabel)
                            .addComponent(topLabel))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(snrTextField)
                            .addComponent(frompTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 42, Short.MAX_VALUE)
                            .addComponent(topTextField)))
                    .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                        .addComponent(hfunLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(sigmHfunRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(tanhHfunRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(linHfunRadioButton)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(mulLabel)
                    .addComponent(stateLabel)
                    .addComponent(nranLabel)
                    .addComponent(scaleLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(scaleTextField)
                    .addComponent(nranTextField)
                    .addComponent(stateTextField)
                    .addComponent(mulTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 36, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                        .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(trsLabel)
                            .addComponent(mskLabel))
                        .addGap(2, 2, 2)
                        .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(mskTextField)
                            .addComponent(trsTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 39, Short.MAX_VALUE)))
                    .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                        .addComponent(newBoxLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(newBoxTextField)))
                .addGap(783, 783, 783))
            .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(normaliseLabel)
                    .addComponent(ranktifyLabel)
                    .addComponent(statsLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                        .addComponent(yesNormaliseRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noNormaliseRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(norm_noiseLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(yesNormNoiseRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noNormNoiseRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(avgLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(yesAvgRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noAvgRadioButton))
                    .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                        .addComponent(yesRanktifyRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noRanktifyRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fileTabLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fileTabTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 227, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(browseFileTabButton))
                    .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                        .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                                .addComponent(yesMergeRadioButton)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(noMergeRadioButton))
                            .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                                .addComponent(yesStatsRadioButton1)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(noStatsRadioButton)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 53, Short.MAX_VALUE)
                        .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                                .addComponent(compareLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(yesCompareRadioButton)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(noCompareRadioButton))
                            .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                                .addComponent(ft2imgLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(yesFt2imgRadioButton)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(noFt2imgRadioButton)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(negLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(yesNegRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noNegRadioButton)))
                .addGap(824, 824, 824))
            .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                .addComponent(mergeLabel)
                .addGap(1298, 1298, 1298))
            .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                .addComponent(frameAvgLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(frameAvgTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(clipLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(clipTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(boxLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(boxTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(innerLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(innerTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 56, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(widthLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(widthTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(802, 802, 802))
            .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                .addComponent(outputParameterFileLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(outputParameterFileTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 227, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(browseOutputParameterFileButton)
                .addContainerGap())
            .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                .addComponent(fracZeroSlider, javax.swing.GroupLayout.PREFERRED_SIZE, 288, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(threshold4BinSlider, javax.swing.GroupLayout.PREFERRED_SIZE, 288, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(761, Short.MAX_VALUE))
        );
        optionalParametersPanelLayout.setVerticalGroup(
            optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(optionalParametersPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(hpLabel)
                    .addComponent(hpTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(frompLabel)
                    .addComponent(frompTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(mulLabel)
                    .addComponent(mulTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(trsLabel)
                    .addComponent(trsTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(lpLabel)
                    .addComponent(lpTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(snrLabel)
                    .addComponent(snrTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(stateLabel)
                    .addComponent(stateTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(mskLabel)
                    .addComponent(mskTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(binarizeLabel)
                    .addComponent(yesBinarizeRadioButton)
                    .addComponent(noBinarizeRadioButton)
                    .addComponent(topLabel)
                    .addComponent(topTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(nranLabel)
                    .addComponent(nranTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(newBoxLabel)
                    .addComponent(newBoxTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(hfunLabel)
                    .addComponent(sigmHfunRadioButton)
                    .addComponent(tanhHfunRadioButton)
                    .addComponent(linHfunRadioButton)
                    .addComponent(scaleLabel)
                    .addComponent(scaleTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(normaliseLabel)
                    .addComponent(yesNormaliseRadioButton)
                    .addComponent(noNormaliseRadioButton)
                    .addComponent(norm_noiseLabel)
                    .addComponent(yesNormNoiseRadioButton)
                    .addComponent(noNormNoiseRadioButton)
                    .addComponent(avgLabel)
                    .addComponent(yesAvgRadioButton)
                    .addComponent(noAvgRadioButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ranktifyLabel)
                    .addComponent(yesRanktifyRadioButton)
                    .addComponent(noRanktifyRadioButton)
                    .addComponent(fileTabLabel)
                    .addComponent(fileTabTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(browseFileTabButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(statsLabel)
                    .addComponent(yesStatsRadioButton1)
                    .addComponent(noStatsRadioButton)
                    .addComponent(noNegRadioButton)
                    .addComponent(yesNegRadioButton)
                    .addComponent(negLabel)
                    .addComponent(noCompareRadioButton)
                    .addComponent(yesCompareRadioButton)
                    .addComponent(compareLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(mergeLabel)
                    .addComponent(yesMergeRadioButton)
                    .addComponent(noMergeRadioButton)
                    .addComponent(noFt2imgRadioButton)
                    .addComponent(yesFt2imgRadioButton)
                    .addComponent(ft2imgLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(frameAvgLabel)
                    .addComponent(frameAvgTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(clipLabel)
                    .addComponent(clipTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(boxLabel)
                    .addComponent(boxTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(innerLabel)
                    .addComponent(innerTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(widthLabel)
                    .addComponent(widthTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(outputParameterFileLabel)
                    .addComponent(outputParameterFileTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(browseOutputParameterFileButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(optionalParametersPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fracZeroSlider, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(threshold4BinSlider, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        optionalScrollPanel.setViewportView(optionalParametersPanel);

        symmetryCTFScrollPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Symmetry and CTF..."));

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
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(dSymmetryGroupComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, 78, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(oSymmetryRadioButton)
                .addGap(72, 72, 72))
        );
        symmetryGroupPanelLayout.setVerticalGroup(
            symmetryGroupPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(cSymmetryGroupComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, 46, javax.swing.GroupLayout.PREFERRED_SIZE)
            .addGroup(symmetryGroupPanelLayout.createSequentialGroup()
                .addComponent(cSymmetryRadioButton)
                .addGap(18, 18, 18)
                .addComponent(tSymmetryRadioButton))
            .addComponent(dSymmetryGroupComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, 46, javax.swing.GroupLayout.PREFERRED_SIZE)
            .addGroup(symmetryGroupPanelLayout.createSequentialGroup()
                .addComponent(dSymmetryRadioButton)
                .addGap(18, 18, 18)
                .addComponent(iSymmetryRadioButton))
            .addComponent(oSymmetryRadioButton)
        );

        ctfPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("CTF setup..."));

        yesCTFRadioButton.setText("Yes");
        yesCTFRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesCTFRadioButtonActionPerformed(evt);
            }
        });

        noCTFRadioButton.setText("No");
        noCTFRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noCTFRadioButtonActionPerformed(evt);
            }
        });

        flipCTFRadioButton.setText("Flip");
        flipCTFRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                flipCTFRadioButtonActionPerformed(evt);
            }
        });

        mulCTFRadioButton.setText("Mul");
        mulCTFRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                mulCTFRadioButtonActionPerformed(evt);
            }
        });

        ctfsqLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        ctfsqLabel.setText("Ctf sq:");

        yesCtfSqRadioButton.setText("Yes");
        yesCtfSqRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesCtfSqRadioButtonActionPerformed(evt);
            }
        });

        noCtfSqRadioButton.setText("No");
        noCtfSqRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noCtfSqRadioButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout ctfPanelLayout = new javax.swing.GroupLayout(ctfPanel);
        ctfPanel.setLayout(ctfPanelLayout);
        ctfPanelLayout.setHorizontalGroup(
            ctfPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(ctfPanelLayout.createSequentialGroup()
                .addGroup(ctfPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(ctfPanelLayout.createSequentialGroup()
                        .addComponent(yesCTFRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noCTFRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(flipCTFRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(mulCTFRadioButton))
                    .addGroup(ctfPanelLayout.createSequentialGroup()
                        .addComponent(ctfsqLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(yesCtfSqRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(noCtfSqRadioButton)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        ctfPanelLayout.setVerticalGroup(
            ctfPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(ctfPanelLayout.createSequentialGroup()
                .addGroup(ctfPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(yesCTFRadioButton)
                    .addComponent(noCTFRadioButton)
                    .addComponent(flipCTFRadioButton)
                    .addComponent(mulCTFRadioButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(ctfPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ctfsqLabel)
                    .addComponent(yesCtfSqRadioButton)
                    .addComponent(noCtfSqRadioButton))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout symmetryCTFPanelLayout = new javax.swing.GroupLayout(symmetryCTFPanel);
        symmetryCTFPanel.setLayout(symmetryCTFPanelLayout);
        symmetryCTFPanelLayout.setHorizontalGroup(
            symmetryCTFPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(symmetryCTFPanelLayout.createSequentialGroup()
                .addComponent(symmetryGroupPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 322, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(ctfPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(243, Short.MAX_VALUE))
        );
        symmetryCTFPanelLayout.setVerticalGroup(
            symmetryCTFPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(symmetryCTFPanelLayout.createSequentialGroup()
                .addGroup(symmetryCTFPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(ctfPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(symmetryGroupPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.PREFERRED_SIZE, 94, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(43, Short.MAX_VALUE))
        );

        symmetryCTFScrollPanel.setViewportView(symmetryCTFPanel);

        fracSlider.setFont(new java.awt.Font("Dialog", 0, 8));
        fracSlider.setMajorTickSpacing(10);
        fracSlider.setMinorTickSpacing(2);
        fracSlider.setPaintLabels(true);
        fracSlider.setPaintTicks(true);
        fracSlider.setSnapToTicks(true);
        fracSlider.setToolTipText("Fraction of particles to include{0.7}");
        fracSlider.setValue(70);
        fracSlider.setBorder(javax.swing.BorderFactory.createTitledBorder("fraction of ptcls to include..."));
        fracSlider.setValueIsAdjusting(true);

        javax.swing.GroupLayout mainSimpleStackopsPanelLayout = new javax.swing.GroupLayout(mainSimpleStackopsPanel);
        mainSimpleStackopsPanel.setLayout(mainSimpleStackopsPanelLayout);
        mainSimpleStackopsPanelLayout.setHorizontalGroup(
            mainSimpleStackopsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mainSimpleStackopsPanelLayout.createSequentialGroup()
                .addGroup(mainSimpleStackopsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, mainSimpleStackopsPanelLayout.createSequentialGroup()
                        .addComponent(acfOtherPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(optionalPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(optionalScrollPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.PREFERRED_SIZE, 647, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, mainSimpleStackopsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                        .addComponent(symmetryCTFScrollPanel, javax.swing.GroupLayout.Alignment.LEADING, 0, 0, Short.MAX_VALUE)
                        .addGroup(javax.swing.GroupLayout.Alignment.LEADING, mainSimpleStackopsPanelLayout.createSequentialGroup()
                            .addComponent(stackDetailsPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                            .addComponent(microscopeDetailsPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addComponent(fracSlider, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.PREFERRED_SIZE, 609, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(886, 886, 886))
        );
        mainSimpleStackopsPanelLayout.setVerticalGroup(
            mainSimpleStackopsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mainSimpleStackopsPanelLayout.createSequentialGroup()
                .addGroup(mainSimpleStackopsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(microscopeDetailsPanel, 0, 177, Short.MAX_VALUE)
                    .addComponent(stackDetailsPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 177, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(fracSlider, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(symmetryCTFScrollPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 143, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addGroup(mainSimpleStackopsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(acfOtherPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(optionalPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(optionalScrollPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(466, 466, 466))
        );

        mainSimpleStackopsScrollPanel.setViewportView(mainSimpleStackopsPanel);

        add(mainSimpleStackopsScrollPanel);
    }// </editor-fold>//GEN-END:initComponents

    private void iSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_iSymmetryRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("i");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.SYM_CLASS, value);
    }//GEN-LAST:event_iSymmetryRadioButtonActionPerformed
    private void tSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_tSymmetryRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("t");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.SYM_CLASS, value);
    }//GEN-LAST:event_tSymmetryRadioButtonActionPerformed
    private void oSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_oSymmetryRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("o");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.SYM_CLASS, value);
    }//GEN-LAST:event_oSymmetryRadioButtonActionPerformed
    private void dSymmetryGroupComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_dSymmetryGroupComboBoxActionPerformed
        //TODO: DONE: Implemented below
    }//GEN-LAST:event_dSymmetryGroupComboBoxActionPerformed
    private void cSymmetryGroupComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cSymmetryGroupComboBoxActionPerformed
        //TODO: DONE implemented below
    }//GEN-LAST:event_cSymmetryGroupComboBoxActionPerformed
    private void yesCTFRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesCTFRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.CTF, value);
    }//GEN-LAST:event_yesCTFRadioButtonActionPerformed
    private void noCTFRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noCTFRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.CTF, value);
    }//GEN-LAST:event_noCTFRadioButtonActionPerformed
    private void flipCTFRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_flipCTFRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("flip");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.CTF, value);
    }//GEN-LAST:event_flipCTFRadioButtonActionPerformed
    private void mulCTFRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_mulCTFRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("mul");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.CTF, value);
    }//GEN-LAST:event_mulCTFRadioButtonActionPerformed
    private void yesShalgnRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesShalgnRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.SHALGN, value);
    }//GEN-LAST:event_yesShalgnRadioButtonActionPerformed
    private void noShalgnRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noShalgnRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.SHALGN, value);
    }//GEN-LAST:event_noShalgnRadioButtonActionPerformed
    private void yesRoalgnRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesRoalgnRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.ROALGN, value);
    }//GEN-LAST:event_yesRoalgnRadioButtonActionPerformed
    private void noRoalgnRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noRoalgnRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.ROALGN, value);
    }//GEN-LAST:event_noRoalgnRadioButtonActionPerformed
    private void yesAcfRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesAcfRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.ACF, value);
    }//GEN-LAST:event_yesAcfRadioButtonActionPerformed
    private void noAcfRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noAcfRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.ACF, value);
    }//GEN-LAST:event_noAcfRadioButtonActionPerformed
    private void yesMasscenRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesMasscenRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.MASSCEN, value);
    }//GEN-LAST:event_yesMasscenRadioButtonActionPerformed
    private void noMasscenRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noMasscenRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.MASSCEN, value);
    }//GEN-LAST:event_noMasscenRadioButtonActionPerformed
    private void yesPhrandRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesPhrandRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.PHRAND, value);
    }//GEN-LAST:event_yesPhrandRadioButtonActionPerformed
    private void noPhrandRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noPhrandRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.PHRAND, value);
    }//GEN-LAST:event_noPhrandRadioButtonActionPerformed
    private void yesVisRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesVisRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.VIS, value);
    }//GEN-LAST:event_yesVisRadioButtonActionPerformed
    private void noVisRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noVisRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.VIS, value);
    }//GEN-LAST:event_noVisRadioButtonActionPerformed
    private void noCtfSqRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noCtfSqRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.CTFSQ, value);
    }//GEN-LAST:event_noCtfSqRadioButtonActionPerformed
    private void yesCtfSqRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesCtfSqRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.CTFSQ, value);
    }//GEN-LAST:event_yesCtfSqRadioButtonActionPerformed
    private void yesBinarizeRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesBinarizeRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.BIN, value);
    }//GEN-LAST:event_yesBinarizeRadioButtonActionPerformed
    private void noBinarizeRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noBinarizeRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.BIN, value);
    }//GEN-LAST:event_noBinarizeRadioButtonActionPerformed
    private void sigmHfunRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_sigmHfunRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("sigm");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.HFUN, value);
    }//GEN-LAST:event_sigmHfunRadioButtonActionPerformed
    private void tanhHfunRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_tanhHfunRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("tanh");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.HFUN, value);
    }//GEN-LAST:event_tanhHfunRadioButtonActionPerformed
    private void linHfunRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_linHfunRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("lin");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.HFUN, value);
    }//GEN-LAST:event_linHfunRadioButtonActionPerformed
    private void yesNormNoiseRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesNormNoiseRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.NOISE_NORM, value);
    }//GEN-LAST:event_yesNormNoiseRadioButtonActionPerformed
    private void noNormNoiseRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noNormNoiseRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.NOISE_NORM, value);
    }//GEN-LAST:event_noNormNoiseRadioButtonActionPerformed
    private void yesCompareRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesCompareRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.COMPARE, value);
    }//GEN-LAST:event_yesCompareRadioButtonActionPerformed
    private void noCompareRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noCompareRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.COMPARE, value);
    }//GEN-LAST:event_noCompareRadioButtonActionPerformed
    private void yesNormaliseRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesNormaliseRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.NORM, value);
    }//GEN-LAST:event_yesNormaliseRadioButtonActionPerformed
    private void noNormaliseRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noNormaliseRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.NORM, value);
    }//GEN-LAST:event_noNormaliseRadioButtonActionPerformed
    private void yesAvgRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesAvgRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.AVG, value);
    }//GEN-LAST:event_yesAvgRadioButtonActionPerformed
    private void noAvgRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noAvgRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.AVG, value);
    }//GEN-LAST:event_noAvgRadioButtonActionPerformed
    private void yesRanktifyRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesRanktifyRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.RANKIFY, value);
    }//GEN-LAST:event_yesRanktifyRadioButtonActionPerformed
    private void noRanktifyRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noRanktifyRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.RANKIFY, value);
    }//GEN-LAST:event_noRanktifyRadioButtonActionPerformed
    private void yesStatsRadioButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesStatsRadioButton1ActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.STATS, value);
    }//GEN-LAST:event_yesStatsRadioButton1ActionPerformed
    private void noStatsRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noStatsRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.STATS, value);
    }//GEN-LAST:event_noStatsRadioButtonActionPerformed
    private void yesNegRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesNegRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.NEG, value);
    }//GEN-LAST:event_yesNegRadioButtonActionPerformed
    private void noNegRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noNegRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.NEG, value);
    }//GEN-LAST:event_noNegRadioButtonActionPerformed
    private void yesMergeRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesMergeRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.MERGE, value);
    }//GEN-LAST:event_yesMergeRadioButtonActionPerformed
    private void noMergeRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noMergeRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.MERGE, value);
    }//GEN-LAST:event_noMergeRadioButtonActionPerformed
    private void yesFt2imgRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesFt2imgRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.FT2IMG, value);
    }//GEN-LAST:event_yesFt2imgRadioButtonActionPerformed
    private void noFt2imgRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noFt2imgRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleStackopsGridBean.FT2IMG, value);
    }//GEN-LAST:event_noFt2imgRadioButtonActionPerformed
    private void cSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cSymmetryRadioButtonActionPerformed
        //C symmetry combo box
        try {
            ComboBoxDataControl cSymComboBoxDataControl = new ComboBoxDataControl(parentPanel.getClient(), SimpleStackopsGridBean.SYM_CLASS, cSymmetryGroupComboBox);
            cSymComboBoxDataControl.setPossibleValues(CSymmetrySimpleStackopsCalculationType.getAllTypes());
            parentPanel.linkDataControl(SimpleStackopsGridBean.SYM_CLASS, cSymComboBoxDataControl);
            parentPanel.setValueTranslator(SimpleStackopsGridBean.SYM_CLASS, StringValueTranslator.getInstance());
            parentPanel.setValueValidator(SimpleStackopsGridBean.SYM_CLASS, NotNullValidator.getInstance());
        } catch (DataSetException e) {
            e.printStackTrace();
        }
    }//GEN-LAST:event_cSymmetryRadioButtonActionPerformed
    private void dSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_dSymmetryRadioButtonActionPerformed
        try {
            //D symmetry combo box
            ComboBoxDataControl dSymComboBoxDataControl = new ComboBoxDataControl(parentPanel.getClient(), SimpleStackopsGridBean.SYM_CLASS, dSymmetryGroupComboBox);
            dSymComboBoxDataControl.setPossibleValues(DSymmetrySimpleStackopsCalculationType.getAllTypes());
            parentPanel.linkDataControl(SimpleStackopsGridBean.SYM_CLASS, dSymComboBoxDataControl);
            parentPanel.setValueTranslator(SimpleStackopsGridBean.SYM_CLASS, StringValueTranslator.getInstance());
            parentPanel.setValueValidator(SimpleStackopsGridBean.SYM_CLASS, NotNullValidator.getInstance());
        } catch (DataSetException e) {
            e.printStackTrace();
        }
    }//GEN-LAST:event_dSymmetryRadioButtonActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel accelerationVoltageLabel;
    private javax.swing.JTextField accelerationVoltageTextField;
    private javax.swing.JLabel acfLabel;
    private javax.swing.JPanel acfOtherPanel;
    private javax.swing.ButtonGroup acfbuttonGroup;
    private javax.swing.JLabel avgLabel;
    private javax.swing.ButtonGroup avgbuttonGroup;
    private javax.swing.JLabel binarizeLabel;
    private javax.swing.ButtonGroup binarizebuttonGroup;
    private javax.swing.JLabel boxLabel;
    private javax.swing.JTextField boxTextField;
    private javax.swing.JButton browseDefTabButton;
    private javax.swing.JButton browseFileTabButton;
    private javax.swing.JButton browseOriTabButton;
    private javax.swing.JButton browseOutStackButton;
    private javax.swing.JButton browseOutputParameterFileButton;
    private javax.swing.JButton browseStk2Button;
    private javax.swing.JButton browseStkButton;
    private javax.swing.JComboBox cSymmetryGroupComboBox;
    private javax.swing.JRadioButton cSymmetryRadioButton;
    private javax.swing.JLabel clipLabel;
    private javax.swing.JTextField clipTextField;
    private javax.swing.JLabel compareLabel;
    private javax.swing.ButtonGroup comparebuttonGroup;
    private javax.swing.JLabel csLabel;
    private javax.swing.JTextField csTextField;
    private javax.swing.JPanel ctfPanel;
    private javax.swing.ButtonGroup ctfbuttonGroup;
    private javax.swing.JLabel ctfsqLabel;
    private javax.swing.ButtonGroup ctfsqbuttonGroup;
    private javax.swing.JComboBox dSymmetryGroupComboBox;
    private javax.swing.JRadioButton dSymmetryRadioButton;
    private javax.swing.JTextField defTabTextField;
    private javax.swing.JLabel deftabLabel;
    private javax.swing.JLabel fileTabLabel;
    private javax.swing.JTextField fileTabTextField;
    private javax.swing.JRadioButton flipCTFRadioButton;
    private javax.swing.JSlider fracSlider;
    private javax.swing.JSlider fracZeroSlider;
    private javax.swing.JLabel fracaLabel;
    private javax.swing.JTextField fracaTextField;
    private javax.swing.JLabel frameAvgLabel;
    private javax.swing.JTextField frameAvgTextField;
    private javax.swing.JLabel frompLabel;
    private javax.swing.JTextField frompTextField;
    private javax.swing.JLabel ft2imgLabel;
    private javax.swing.ButtonGroup ft2imgbuttonGroup;
    private javax.swing.JLabel hfunLabel;
    private javax.swing.ButtonGroup hfunbuttonGroup;
    private javax.swing.JLabel hpLabel;
    private javax.swing.JTextField hpTextField;
    private javax.swing.JRadioButton iSymmetryRadioButton;
    private javax.swing.JLabel innerLabel;
    private javax.swing.JTextField innerTextField;
    private javax.swing.JRadioButton linHfunRadioButton;
    private javax.swing.JLabel lpLabel;
    private javax.swing.JTextField lpTextField;
    private javax.swing.JPanel mainSimpleStackopsPanel;
    private javax.swing.JScrollPane mainSimpleStackopsScrollPanel;
    private javax.swing.JLabel masscenLabel;
    private javax.swing.ButtonGroup masscenbuttonGroup;
    private javax.swing.JLabel mergeLabel;
    private javax.swing.ButtonGroup mergebuttonGroup;
    private javax.swing.JPanel microscopeDetailsPanel;
    private javax.swing.JLabel mskLabel;
    private javax.swing.JTextField mskTextField;
    private javax.swing.JRadioButton mulCTFRadioButton;
    private javax.swing.JLabel mulLabel;
    private javax.swing.JTextField mulTextField;
    private javax.swing.JLabel negLabel;
    private javax.swing.ButtonGroup negbuttonGroup;
    private javax.swing.JLabel newBoxLabel;
    private javax.swing.JTextField newBoxTextField;
    private javax.swing.JRadioButton noAcfRadioButton;
    private javax.swing.JRadioButton noAvgRadioButton;
    private javax.swing.JRadioButton noBinarizeRadioButton;
    private javax.swing.JRadioButton noCTFRadioButton;
    private javax.swing.JRadioButton noCompareRadioButton;
    private javax.swing.JRadioButton noCtfSqRadioButton;
    private javax.swing.JRadioButton noFt2imgRadioButton;
    private javax.swing.JRadioButton noMasscenRadioButton;
    private javax.swing.JRadioButton noMergeRadioButton;
    private javax.swing.JRadioButton noNegRadioButton;
    private javax.swing.JRadioButton noNormNoiseRadioButton;
    private javax.swing.JRadioButton noNormaliseRadioButton;
    private javax.swing.JRadioButton noPhrandRadioButton;
    private javax.swing.JRadioButton noRanktifyRadioButton;
    private javax.swing.JRadioButton noRoalgnRadioButton;
    private javax.swing.JRadioButton noShalgnRadioButton;
    private javax.swing.JRadioButton noStatsRadioButton;
    private javax.swing.JRadioButton noVisRadioButton;
    private javax.swing.ButtonGroup normNoisebuttonGroup;
    private javax.swing.JLabel norm_noiseLabel;
    private javax.swing.JLabel normaliseLabel;
    private javax.swing.ButtonGroup normalisebuttonGroup;
    private javax.swing.JLabel nptclsLabel;
    private javax.swing.JTextField nptclsTextField;
    private javax.swing.JLabel nranLabel;
    private javax.swing.JTextField nranTextField;
    private javax.swing.JRadioButton oSymmetryRadioButton;
    private javax.swing.JPanel optionalPanel;
    private javax.swing.JPanel optionalParametersPanel;
    private javax.swing.JScrollPane optionalScrollPanel;
    private javax.swing.JTextField oriTabTextField;
    private javax.swing.JLabel oritabLabel;
    private javax.swing.JLabel outputParameterFileLabel;
    private javax.swing.JTextField outputParameterFileTextField;
    private javax.swing.JLabel outstackLabel;
    private javax.swing.JTextField outstackTextField;
    private javax.swing.JLabel phrandLabel;
    private javax.swing.ButtonGroup phrandbuttonGroup;
    private javax.swing.JLabel ranktifyLabel;
    private javax.swing.ButtonGroup ranktifybuttonGroup;
    private javax.swing.JLabel roalgnLabel;
    private javax.swing.ButtonGroup roalgnbuttonGroup;
    private javax.swing.JLabel scaleLabel;
    private javax.swing.JTextField scaleTextField;
    private javax.swing.JLabel shalgnLabel;
    private javax.swing.ButtonGroup shiftAlignmentbuttonGroup;
    private javax.swing.JRadioButton sigmHfunRadioButton;
    private javax.swing.JLabel smpdLabel;
    private javax.swing.JTextField smpdTextField;
    private javax.swing.JLabel snrLabel;
    private javax.swing.JTextField snrTextField;
    private javax.swing.JTextField stack2NameTextField;
    private javax.swing.JPanel stackDetailsPanel;
    private javax.swing.JTextField stackNameTextField;
    private javax.swing.JLabel stak2NameLabel;
    private javax.swing.JLabel stakNameLabel;
    private javax.swing.JLabel stateLabel;
    private javax.swing.JTextField stateTextField;
    private javax.swing.JLabel statsLabel;
    private javax.swing.ButtonGroup statsbuttonGroup;
    private javax.swing.JPanel symmetryCTFPanel;
    private javax.swing.JScrollPane symmetryCTFScrollPanel;
    private javax.swing.JPanel symmetryGroupPanel;
    private javax.swing.ButtonGroup symmetryGroupbuttonGroup;
    private javax.swing.JRadioButton tSymmetryRadioButton;
    private javax.swing.JRadioButton tanhHfunRadioButton;
    private javax.swing.JSlider threshold4BinSlider;
    private javax.swing.JLabel topLabel;
    private javax.swing.JTextField topTextField;
    private javax.swing.JLabel trsLabel;
    private javax.swing.JTextField trsTextField;
    private javax.swing.JLabel visLabel;
    private javax.swing.ButtonGroup visbuttonGroup;
    private javax.swing.JLabel widthLabel;
    private javax.swing.JTextField widthTextField;
    private javax.swing.JRadioButton yesAcfRadioButton;
    private javax.swing.JRadioButton yesAvgRadioButton;
    private javax.swing.JRadioButton yesBinarizeRadioButton;
    private javax.swing.JRadioButton yesCTFRadioButton;
    private javax.swing.JRadioButton yesCompareRadioButton;
    private javax.swing.JRadioButton yesCtfSqRadioButton;
    private javax.swing.JRadioButton yesFt2imgRadioButton;
    private javax.swing.JRadioButton yesMasscenRadioButton;
    private javax.swing.JRadioButton yesMergeRadioButton;
    private javax.swing.JRadioButton yesNegRadioButton;
    private javax.swing.JRadioButton yesNormNoiseRadioButton;
    private javax.swing.JRadioButton yesNormaliseRadioButton;
    private javax.swing.JRadioButton yesPhrandRadioButton;
    private javax.swing.JRadioButton yesRanktifyRadioButton;
    private javax.swing.JRadioButton yesRoalgnRadioButton;
    private javax.swing.JRadioButton yesShalgnRadioButton;
    private javax.swing.JRadioButton yesStatsRadioButton1;
    private javax.swing.JRadioButton yesVisRadioButton;
    // End of variables declaration//GEN-END:variables

}
