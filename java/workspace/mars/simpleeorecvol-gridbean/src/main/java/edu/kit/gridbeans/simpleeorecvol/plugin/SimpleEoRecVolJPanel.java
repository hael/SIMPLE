/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * SimpleEoRecVolJPanel.java
 *
 * Created on 20/11/2015, 12:01:09 PM
 */

package edu.kit.gridbeans.simpleeorecvol.plugin;

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

import edu.kit.gridbeans.simpleeorecvol.CSymmetrySimpleEoRecVolCalculationType;
import edu.kit.gridbeans.simpleeorecvol.DSymmetrySimpleEoRecVolCalculationType;
import edu.kit.gridbeans.simpleeorecvol.SimpleEoRecVolGridBean;

/**
 *
 * @author frederic
 */
public class SimpleEoRecVolJPanel extends JPanel {
    private GridBeanPanel parentPanel;
    /**
     * Add UID
     */
    private static final long serialVersionUID = -3724146872530150564L;
    private int fracslidernumber;
    /** Creates new form SimpleEoRecVolJPanel */
    public SimpleEoRecVolJPanel(GridBeanPanel parentpanel) throws DataSetException {
        this.parentPanel = parentpanel;
        initComponents();
        initCustomComponents();// custom code by us
        bindComponents();
    }

    private void initCustomComponents() {
        //CTF button group
        ctfbuttonGroup.add(yesCTFRadioButton);
        ctfbuttonGroup.add(noCTFRadioButton);
        ctfbuttonGroup.add(flipCTFRadioButton);
        ctfbuttonGroup.add(mulCTFRadioButton);
        //Symmetry button group
        symmetryGroupbuttonGroup.add(cSymmetryRadioButton);
        symmetryGroupbuttonGroup.add(dSymmetryRadioButton);
        symmetryGroupbuttonGroup.add(tSymmetryRadioButton);
        symmetryGroupbuttonGroup.add(iSymmetryRadioButton);
        symmetryGroupbuttonGroup.add(oSymmetryRadioButton);
    }

    private void bindComponents() throws DataSetException {
        //STK
        parentPanel.linkTextField(SimpleEoRecVolGridBean.STK, stackNameTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.STK, StringValueTranslator.getInstance());
        //Orientation table
        parentPanel.linkTextField(SimpleEoRecVolGridBean.ORITAB, oriTabTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.ORITAB, StringValueTranslator.getInstance());
        //Sampling distance
        parentPanel.linkTextField(SimpleEoRecVolGridBean.SMPD, smpdTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.SMPD, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.SMPD, new DoubleValueValidator());
        //Masking value
        parentPanel.linkTextField(SimpleEoRecVolGridBean.MSK, mskTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.MSK, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.MSK, new DoubleValueValidator());
        //Acceleration volatage in KeV
        parentPanel.linkTextField(SimpleEoRecVolGridBean.KV, accelerationVoltageTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.KV, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.KV, new DoubleValueValidator());
        //--mul=$MUL \               #[mul=<shift multiplication factor{1}>]
        parentPanel.linkTextField(SimpleEoRecVolGridBean.MUL, mulTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.MUL, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.MUL, new DoubleValueValidator());
        //--fromp=$FROMP \                  # [fromp=<start ptcl{1)>]
        parentPanel.linkTextField(SimpleEoRecVolGridBean.FROMP, frompTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.FROMP, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.FROMP, new IntegerValueValidator());
        //--top=$TOP \                      # [top=<stop ptcl{nptcls}>]
        parentPanel.linkTextField(SimpleEoRecVolGridBean.TOP, topTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.TOP, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.TOP, new IntegerValueValidator());
        //frac amp contrast{0.07}
        parentPanel.linkTextField(SimpleEoRecVolGridBean.FRACA, fracaTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.FRACA, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.FRACA, new DoubleValueValidator());
        //Shperical abberation
        parentPanel.linkTextField(SimpleEoRecVolGridBean.CS, csTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.CS, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.CS, new DoubleValueValidator());
        //Defocus values
        parentPanel.linkTextField(SimpleEoRecVolGridBean.DEFTAB, defTabTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.DEFTAB, StringValueTranslator.getInstance());
        //--state=$STATE \           #[state=<state to extract>]
        parentPanel.linkTextField(SimpleEoRecVolGridBean.STATE, stateTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.STATE, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.STATE, new IntegerValueValidator());
        //--alpha=$ALPHA \                  # [alpha=<oversampling ratio or padding factor{2.}>]
        parentPanel.linkTextField(SimpleEoRecVolGridBean.ALPHA, alphaTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.ALPHA, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.ALPHA, new DoubleValueValidator());
        //--mw=$MW \ #[mw=<molecular weight(in kD)>]
        parentPanel.linkTextField(SimpleEoRecVolGridBean.MW, mwTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.MW, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.MW, new DoubleValueValidator());
        //--amsklp=$AMSKLP \                 #[amsklp=<automask low-pass limit(in A)>] 
        parentPanel.linkTextField(SimpleEoRecVolGridBean.AMSKLP, amsklpTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.AMSKLP, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.AMSKLP, new DoubleValueValidator());
        //--edge=$EDGE \                     #[edge=<edge size for softening molecular envelope(in pixels){3}>] 
        parentPanel.linkTextField(SimpleEoRecVolGridBean.EDGE, edgeTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.EDGE, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.EDGE, new IntegerValueValidator());
        //--dens=$DENS \                     #[dens=<density(e.g. 9.368 Da/A3 4 gold clusters){0.}>]
        parentPanel.linkTextField(SimpleEoRecVolGridBean.DENS, denseTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.DENS, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.DENS, new DoubleValueValidator());
        //--inner=$INNER \                   #[inner=<inner mask radius(in pixels)>]
        parentPanel.linkTextField(SimpleEoRecVolGridBean.INNER, innerTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.INNER, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.INNER, new DoubleValueValidator());
        //--box=$BOX \                      # [box=<image size(in pixels)>]
        parentPanel.linkTextField(SimpleEoRecVolGridBean.BOX, boxTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.BOX, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.BOX, new DoubleValueValidator());
        //--width=$WIDTH \                   #[width=<pixels falloff inner mask{10}>]
        parentPanel.linkTextField(SimpleEoRecVolGridBean.WIDTH, widthTextField);
        parentPanel.setValueTranslator(SimpleEoRecVolGridBean.WIDTH, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimpleEoRecVolGridBean.WIDTH, new DoubleValueValidator());
        //link to the JSlider for the frac
        class FracSliderNumberListener implements ChangeListener {
            public void stateChanged(ChangeEvent event) {
                setFracSliderNumber();
            }
        }
        ChangeListener verbosityperfslidernumberlistener = new FracSliderNumberListener();
        fracSlider.addChangeListener(verbosityperfslidernumberlistener);

    }
    /**
     * method for the JSlider assignment method for frac
     */
    public void setFracSliderNumber() {
        fracslidernumber = fracSlider.getValue();

        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue(Integer.toString(fracslidernumber));
        parentPanel.getGridBeanModel().set(SimpleEoRecVolGridBean.FRAC, value);
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        symmetryGroupbuttonGroup = new javax.swing.ButtonGroup();
        ctfbuttonGroup = new javax.swing.ButtonGroup();
        mainSimpleEoRecVolScrollPanel = new javax.swing.JScrollPane();
        mainSimpleEoRecVolPanel = new javax.swing.JPanel();
        microscopeDetailsPanel = new javax.swing.JPanel();
        smpdLabel = new javax.swing.JLabel();
        smpdTextField = new javax.swing.JTextField();
        accelerationVoltageLabel = new javax.swing.JLabel();
        accelerationVoltageTextField = new javax.swing.JTextField();
        csLabel = new javax.swing.JLabel();
        csTextField = new javax.swing.JTextField();
        fracaLabel = new javax.swing.JLabel();
        fracaTextField = new javax.swing.JTextField();
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
        boxLabel = new javax.swing.JLabel();
        boxTextField = new javax.swing.JTextField();
        parameterOptionalPanel = new javax.swing.JPanel();
        fracSlider = new javax.swing.JSlider();
        mskLabel = new javax.swing.JLabel();
        mskTextField = new javax.swing.JTextField();
        mulLabel = new javax.swing.JLabel();
        mulTextField = new javax.swing.JTextField();
        frompLabel = new javax.swing.JLabel();
        frompTextField = new javax.swing.JTextField();
        topLabel = new javax.swing.JLabel();
        topTextField = new javax.swing.JTextField();
        ctfPanel = new javax.swing.JPanel();
        yesCTFRadioButton = new javax.swing.JRadioButton();
        noCTFRadioButton = new javax.swing.JRadioButton();
        flipCTFRadioButton = new javax.swing.JRadioButton();
        mulCTFRadioButton = new javax.swing.JRadioButton();
        paramPanel = new javax.swing.JPanel();
        stateLabel = new javax.swing.JLabel();
        stateTextField = new javax.swing.JTextField();
        alphaLabel = new javax.swing.JLabel();
        alphaTextField = new javax.swing.JTextField();
        otherParamsPanel = new javax.swing.JPanel();
        mwLabel = new javax.swing.JLabel();
        mwTextField = new javax.swing.JTextField();
        amsklpLabel = new javax.swing.JLabel();
        amsklpTextField = new javax.swing.JTextField();
        edgeLabel = new javax.swing.JLabel();
        edgeTextField = new javax.swing.JTextField();
        densLabel = new javax.swing.JLabel();
        denseTextField = new javax.swing.JTextField();
        innerLabel = new javax.swing.JLabel();
        innerTextField = new javax.swing.JTextField();
        widthLabel = new javax.swing.JLabel();
        widthTextField = new javax.swing.JTextField();
        symmetryGroupPanel = new javax.swing.JPanel();
        cSymmetryGroupComboBox = new javax.swing.JComboBox();
        dSymmetryGroupComboBox = new javax.swing.JComboBox();
        tSymmetryRadioButton = new javax.swing.JRadioButton();
        oSymmetryRadioButton = new javax.swing.JRadioButton();
        iSymmetryRadioButton = new javax.swing.JRadioButton();
        cSymmetryRadioButton = new javax.swing.JRadioButton();
        dSymmetryRadioButton = new javax.swing.JRadioButton();

        setLayout(new java.awt.BorderLayout());

        microscopeDetailsPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("MIcroscope details..."));

        smpdLabel.setText("Sampling distance smpd:");

        accelerationVoltageLabel.setText("acceleration voltage(in kV) kv:");

        csLabel.setText("Spherical Aberration cs:");

        fracaLabel.setText("Ampl. Contrast fraca:");

        javax.swing.GroupLayout microscopeDetailsPanelLayout = new javax.swing.GroupLayout(microscopeDetailsPanel);
        microscopeDetailsPanel.setLayout(microscopeDetailsPanelLayout);
        microscopeDetailsPanelLayout.setHorizontalGroup(
            microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(microscopeDetailsPanelLayout.createSequentialGroup()
                .addGroup(microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(microscopeDetailsPanelLayout.createSequentialGroup()
                        .addComponent(smpdLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(smpdTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 106, Short.MAX_VALUE))
                    .addGroup(microscopeDetailsPanelLayout.createSequentialGroup()
                        .addComponent(accelerationVoltageLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(accelerationVoltageTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 77, Short.MAX_VALUE))
                    .addGroup(microscopeDetailsPanelLayout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(csLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(csTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 99, Short.MAX_VALUE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, microscopeDetailsPanelLayout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(fracaLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fracaTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 110, Short.MAX_VALUE)))
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
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 16, Short.MAX_VALUE)
                .addGroup(microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(csTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(csLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(microscopeDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(fracaTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(fracaLabel))
                .addContainerGap())
        );

        stackDetailsPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Stack and input details..."));

        stakNameLabel.setText("Stk name:");

        browseStkButton.setText("Browse");

        oritabLabel.setText("Ori Table:");

        browseOriTabButton.setText("Browse");

        deftabLabel.setText("Def Table:");

        browseDefTabButton.setText("Browse");

        boxLabel.setText("Box size:");

        javax.swing.GroupLayout stackDetailsPanelLayout = new javax.swing.GroupLayout(stackDetailsPanel);
        stackDetailsPanel.setLayout(stackDetailsPanelLayout);
        stackDetailsPanelLayout.setHorizontalGroup(
            stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(stackDetailsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(stackDetailsPanelLayout.createSequentialGroup()
                        .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(oritabLabel)
                            .addComponent(stakNameLabel))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(stackDetailsPanelLayout.createSequentialGroup()
                                .addComponent(stackNameTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 99, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(6, 6, 6))
                            .addGroup(stackDetailsPanelLayout.createSequentialGroup()
                                .addComponent(oriTabTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 99, Short.MAX_VALUE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)))
                        .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(stackDetailsPanelLayout.createSequentialGroup()
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(browseStkButton)
                                .addContainerGap(38, Short.MAX_VALUE))
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, stackDetailsPanelLayout.createSequentialGroup()
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(browseOriTabButton)
                                .addGap(151, 151, 151))))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, stackDetailsPanelLayout.createSequentialGroup()
                        .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(javax.swing.GroupLayout.Alignment.LEADING, stackDetailsPanelLayout.createSequentialGroup()
                                .addComponent(boxLabel)
                                .addGap(18, 18, 18)
                                .addComponent(boxTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 99, Short.MAX_VALUE))
                            .addGroup(stackDetailsPanelLayout.createSequentialGroup()
                                .addComponent(deftabLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(defTabTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 99, Short.MAX_VALUE)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(browseDefTabButton)
                        .addGap(151, 151, 151))))
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
                    .addComponent(oritabLabel)
                    .addComponent(browseOriTabButton)
                    .addComponent(oriTabTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(deftabLabel)
                    .addComponent(defTabTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(browseDefTabButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 10, Short.MAX_VALUE)
                .addGroup(stackDetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(boxLabel)
                    .addComponent(boxTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );

        parameterOptionalPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Parameter optional..."));

        fracSlider.setFont(new java.awt.Font("Dialog", 0, 8));
        fracSlider.setMajorTickSpacing(10);
        fracSlider.setMinorTickSpacing(2);
        fracSlider.setPaintLabels(true);
        fracSlider.setPaintTicks(true);
        fracSlider.setSnapToTicks(true);
        fracSlider.setToolTipText("value for the nsteps variable intermediate steps between reports");
        fracSlider.setValue(70);
        fracSlider.setBorder(javax.swing.BorderFactory.createTitledBorder("fraction of ptcls to include..."));
        fracSlider.setValueIsAdjusting(true);

        mskLabel.setText("Mask radius msk:");

        mulLabel.setText("Shift multipication mul:");

        frompLabel.setText("fromp:");

        topLabel.setText("Top:");

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

        javax.swing.GroupLayout ctfPanelLayout = new javax.swing.GroupLayout(ctfPanel);
        ctfPanel.setLayout(ctfPanelLayout);
        ctfPanelLayout.setHorizontalGroup(
            ctfPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(ctfPanelLayout.createSequentialGroup()
                .addComponent(yesCTFRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(noCTFRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(flipCTFRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(mulCTFRadioButton)
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
                .addContainerGap(23, Short.MAX_VALUE))
        );

        paramPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("More params..."));

        stateLabel.setText("State:");

        alphaLabel.setText("Sampling ratio alpha:");

        javax.swing.GroupLayout paramPanelLayout = new javax.swing.GroupLayout(paramPanel);
        paramPanel.setLayout(paramPanelLayout);
        paramPanelLayout.setHorizontalGroup(
            paramPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(paramPanelLayout.createSequentialGroup()
                .addComponent(stateLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(stateTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 49, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(alphaLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(alphaTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 41, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(22, Short.MAX_VALUE))
        );
        paramPanelLayout.setVerticalGroup(
            paramPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(paramPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(paramPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(stateLabel)
                    .addComponent(stateTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(alphaLabel)
                    .addComponent(alphaTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        otherParamsPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Other params..."));

        mwLabel.setText("Molecular weight mw:");

        amsklpLabel.setText("Low pass limit amsklp:");

        edgeLabel.setText("edge:");

        densLabel.setText("Density dens:");

        innerLabel.setText("Inner mask radius inner:");

        widthLabel.setText("pixels falloff width:");

        javax.swing.GroupLayout otherParamsPanelLayout = new javax.swing.GroupLayout(otherParamsPanel);
        otherParamsPanel.setLayout(otherParamsPanelLayout);
        otherParamsPanelLayout.setHorizontalGroup(
            otherParamsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(otherParamsPanelLayout.createSequentialGroup()
                .addGroup(otherParamsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(otherParamsPanelLayout.createSequentialGroup()
                        .addComponent(densLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(denseTextField))
                    .addComponent(mwLabel))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(otherParamsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(otherParamsPanelLayout.createSequentialGroup()
                        .addComponent(mwTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 34, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(amsklpLabel))
                    .addGroup(otherParamsPanelLayout.createSequentialGroup()
                        .addComponent(innerLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(innerTextField)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(otherParamsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(otherParamsPanelLayout.createSequentialGroup()
                        .addComponent(amsklpTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 35, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(edgeLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(edgeTextField))
                    .addGroup(otherParamsPanelLayout.createSequentialGroup()
                        .addComponent(widthLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(widthTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 42, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(112, Short.MAX_VALUE))
        );
        otherParamsPanelLayout.setVerticalGroup(
            otherParamsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(otherParamsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(otherParamsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(mwLabel)
                    .addComponent(mwTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(amsklpLabel)
                    .addComponent(amsklpTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(edgeLabel)
                    .addComponent(edgeTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(otherParamsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(densLabel)
                    .addComponent(denseTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(innerLabel)
                    .addComponent(innerTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(widthLabel)
                    .addComponent(widthTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(16, Short.MAX_VALUE))
        );

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

        javax.swing.GroupLayout parameterOptionalPanelLayout = new javax.swing.GroupLayout(parameterOptionalPanel);
        parameterOptionalPanel.setLayout(parameterOptionalPanelLayout);
        parameterOptionalPanelLayout.setHorizontalGroup(
            parameterOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(parameterOptionalPanelLayout.createSequentialGroup()
                .addGroup(parameterOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(otherParamsPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, parameterOptionalPanelLayout.createSequentialGroup()
                        .addComponent(fracSlider, javax.swing.GroupLayout.PREFERRED_SIZE, 240, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(symmetryGroupPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 322, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, parameterOptionalPanelLayout.createSequentialGroup()
                        .addComponent(mskLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(mskTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(mulLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(mulTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 44, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(frompLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(frompTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 42, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(topLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(topTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 42, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, parameterOptionalPanelLayout.createSequentialGroup()
                        .addComponent(ctfPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(paramPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap())
        );
        parameterOptionalPanelLayout.setVerticalGroup(
            parameterOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, parameterOptionalPanelLayout.createSequentialGroup()
                .addGroup(parameterOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(fracSlider, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(symmetryGroupPanel, javax.swing.GroupLayout.DEFAULT_SIZE, 94, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED, 24, Short.MAX_VALUE)
                .addGroup(parameterOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(mskLabel)
                    .addComponent(mskTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(mulLabel)
                    .addComponent(mulTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(frompLabel)
                    .addComponent(frompTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(topLabel)
                    .addComponent(topTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(parameterOptionalPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(parameterOptionalPanelLayout.createSequentialGroup()
                        .addGap(1, 1, 1)
                        .addComponent(paramPanel, 0, 69, Short.MAX_VALUE))
                    .addComponent(ctfPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(otherParamsPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );

        fracSlider.getAccessibleContext().setAccessibleDescription("Fraction of particles to includes 0.7 is the default");

        javax.swing.GroupLayout mainSimpleEoRecVolPanelLayout = new javax.swing.GroupLayout(mainSimpleEoRecVolPanel);
        mainSimpleEoRecVolPanel.setLayout(mainSimpleEoRecVolPanelLayout);
        mainSimpleEoRecVolPanelLayout.setHorizontalGroup(
            mainSimpleEoRecVolPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mainSimpleEoRecVolPanelLayout.createSequentialGroup()
                .addGroup(mainSimpleEoRecVolPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addComponent(parameterOptionalPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, mainSimpleEoRecVolPanelLayout.createSequentialGroup()
                        .addComponent(microscopeDetailsPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(stackDetailsPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 278, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(74, Short.MAX_VALUE))
        );
        mainSimpleEoRecVolPanelLayout.setVerticalGroup(
            mainSimpleEoRecVolPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mainSimpleEoRecVolPanelLayout.createSequentialGroup()
                .addGroup(mainSimpleEoRecVolPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(stackDetailsPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(microscopeDetailsPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(parameterOptionalPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(1308, 1308, 1308))
        );

        mainSimpleEoRecVolScrollPanel.setViewportView(mainSimpleEoRecVolPanel);

        add(mainSimpleEoRecVolScrollPanel, java.awt.BorderLayout.CENTER);
    }// </editor-fold>//GEN-END:initComponents

    private void yesCTFRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesCTFRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("yes");
        parentPanel.getGridBeanModel().set(SimpleEoRecVolGridBean.CTF, value);
    }//GEN-LAST:event_yesCTFRadioButtonActionPerformed
    private void noCTFRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noCTFRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimpleEoRecVolGridBean.CTF, value);
    }//GEN-LAST:event_noCTFRadioButtonActionPerformed
    private void flipCTFRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_flipCTFRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("flip");
        parentPanel.getGridBeanModel().set(SimpleEoRecVolGridBean.CTF, value);
    }//GEN-LAST:event_flipCTFRadioButtonActionPerformed
    private void mulCTFRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_mulCTFRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("mul");
        parentPanel.getGridBeanModel().set(SimpleEoRecVolGridBean.CTF, value);
    }//GEN-LAST:event_mulCTFRadioButtonActionPerformed
    private void cSymmetryGroupComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cSymmetryGroupComboBoxActionPerformed
        //TODO: DONE implemented below
    }//GEN-LAST:event_cSymmetryGroupComboBoxActionPerformed
    private void dSymmetryGroupComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_dSymmetryGroupComboBoxActionPerformed
        //TODO: DONE: Implemented below
    }//GEN-LAST:event_dSymmetryGroupComboBoxActionPerformed
    private void tSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_tSymmetryRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("t");
        parentPanel.getGridBeanModel().set(SimpleEoRecVolGridBean.PGRP, value);
    }//GEN-LAST:event_tSymmetryRadioButtonActionPerformed
    private void oSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_oSymmetryRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("o");
        parentPanel.getGridBeanModel().set(SimpleEoRecVolGridBean.PGRP, value);
    }//GEN-LAST:event_oSymmetryRadioButtonActionPerformed
    private void iSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_iSymmetryRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("i");
        parentPanel.getGridBeanModel().set(SimpleEoRecVolGridBean.PGRP, value);
    }//GEN-LAST:event_iSymmetryRadioButtonActionPerformed
    private void cSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cSymmetryRadioButtonActionPerformed
        //C symmetry combo box
        try {
            ComboBoxDataControl cSymComboBoxDataControl = new ComboBoxDataControl(parentPanel.getClient(), SimpleEoRecVolGridBean.PGRP, cSymmetryGroupComboBox);
            cSymComboBoxDataControl.setPossibleValues(CSymmetrySimpleEoRecVolCalculationType.getAllTypes());
            parentPanel.linkDataControl(SimpleEoRecVolGridBean.PGRP, cSymComboBoxDataControl);
            parentPanel.setValueTranslator(SimpleEoRecVolGridBean.PGRP, StringValueTranslator.getInstance());
            parentPanel.setValueValidator(SimpleEoRecVolGridBean.PGRP, NotNullValidator.getInstance());
        } catch (DataSetException e) {
            e.printStackTrace();
        }
    }//GEN-LAST:event_cSymmetryRadioButtonActionPerformed
    private void dSymmetryRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_dSymmetryRadioButtonActionPerformed
        try {
            //D symmetry combo box
            ComboBoxDataControl dSymComboBoxDataControl = new ComboBoxDataControl(parentPanel.getClient(), SimpleEoRecVolGridBean.PGRP, dSymmetryGroupComboBox);
            dSymComboBoxDataControl.setPossibleValues(DSymmetrySimpleEoRecVolCalculationType.getAllTypes());
            parentPanel.linkDataControl(SimpleEoRecVolGridBean.PGRP, dSymComboBoxDataControl);
            parentPanel.setValueTranslator(SimpleEoRecVolGridBean.PGRP, StringValueTranslator.getInstance());
            parentPanel.setValueValidator(SimpleEoRecVolGridBean.PGRP, NotNullValidator.getInstance());
        } catch (DataSetException e) {
            e.printStackTrace();
        }
    }//GEN-LAST:event_dSymmetryRadioButtonActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel accelerationVoltageLabel;
    private javax.swing.JTextField accelerationVoltageTextField;
    private javax.swing.JLabel alphaLabel;
    private javax.swing.JTextField alphaTextField;
    private javax.swing.JLabel amsklpLabel;
    private javax.swing.JTextField amsklpTextField;
    private javax.swing.JLabel boxLabel;
    private javax.swing.JTextField boxTextField;
    private javax.swing.JButton browseDefTabButton;
    private javax.swing.JButton browseOriTabButton;
    private javax.swing.JButton browseStkButton;
    private javax.swing.JComboBox cSymmetryGroupComboBox;
    private javax.swing.JRadioButton cSymmetryRadioButton;
    private javax.swing.JLabel csLabel;
    private javax.swing.JTextField csTextField;
    private javax.swing.JPanel ctfPanel;
    private javax.swing.ButtonGroup ctfbuttonGroup;
    private javax.swing.JComboBox dSymmetryGroupComboBox;
    private javax.swing.JRadioButton dSymmetryRadioButton;
    private javax.swing.JTextField defTabTextField;
    private javax.swing.JLabel deftabLabel;
    private javax.swing.JLabel densLabel;
    private javax.swing.JTextField denseTextField;
    private javax.swing.JLabel edgeLabel;
    private javax.swing.JTextField edgeTextField;
    private javax.swing.JRadioButton flipCTFRadioButton;
    private javax.swing.JSlider fracSlider;
    private javax.swing.JLabel fracaLabel;
    private javax.swing.JTextField fracaTextField;
    private javax.swing.JLabel frompLabel;
    private javax.swing.JTextField frompTextField;
    private javax.swing.JRadioButton iSymmetryRadioButton;
    private javax.swing.JLabel innerLabel;
    private javax.swing.JTextField innerTextField;
    private javax.swing.JPanel mainSimpleEoRecVolPanel;
    private javax.swing.JScrollPane mainSimpleEoRecVolScrollPanel;
    private javax.swing.JPanel microscopeDetailsPanel;
    private javax.swing.JLabel mskLabel;
    private javax.swing.JTextField mskTextField;
    private javax.swing.JRadioButton mulCTFRadioButton;
    private javax.swing.JLabel mulLabel;
    private javax.swing.JTextField mulTextField;
    private javax.swing.JLabel mwLabel;
    private javax.swing.JTextField mwTextField;
    private javax.swing.JRadioButton noCTFRadioButton;
    private javax.swing.JRadioButton oSymmetryRadioButton;
    private javax.swing.JTextField oriTabTextField;
    private javax.swing.JLabel oritabLabel;
    private javax.swing.JPanel otherParamsPanel;
    private javax.swing.JPanel paramPanel;
    private javax.swing.JPanel parameterOptionalPanel;
    private javax.swing.JLabel smpdLabel;
    private javax.swing.JTextField smpdTextField;
    private javax.swing.JPanel stackDetailsPanel;
    private javax.swing.JTextField stackNameTextField;
    private javax.swing.JLabel stakNameLabel;
    private javax.swing.JLabel stateLabel;
    private javax.swing.JTextField stateTextField;
    private javax.swing.JPanel symmetryGroupPanel;
    private javax.swing.ButtonGroup symmetryGroupbuttonGroup;
    private javax.swing.JRadioButton tSymmetryRadioButton;
    private javax.swing.JLabel topLabel;
    private javax.swing.JTextField topTextField;
    private javax.swing.JLabel widthLabel;
    private javax.swing.JTextField widthTextField;
    private javax.swing.JRadioButton yesCTFRadioButton;
    // End of variables declaration//GEN-END:variables

}
