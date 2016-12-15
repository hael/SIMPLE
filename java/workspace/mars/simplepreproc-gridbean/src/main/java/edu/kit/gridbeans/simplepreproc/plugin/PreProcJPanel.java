/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * PreProcJPanel.java
 *
 * Created on 19/11/2015, 8:18:58 PM
 */

package edu.kit.gridbeans.simplepreproc.plugin;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JPanel;

import com.intel.gpe.gridbeans.parameters.EnvironmentVariableParameterValue;
import com.intel.gpe.gridbeans.plugins.DataSetException;
import com.intel.gpe.gridbeans.plugins.swing.panels.GridBeanPanel;
import com.intel.gpe.gridbeans.plugins.translators.StringValueTranslator;
import com.intel.gpe.gridbeans.plugins.validators.DoubleValueValidator;
import com.intel.gpe.gridbeans.plugins.validators.IntegerValueValidator;

import edu.kit.gridbeans.simplepreproc.SimplePreProcGridBean;

/**
 *
 * @author Frederic Bonnet
 * @Date 19th of November 2015
 * 
 */
public class PreProcJPanel extends JPanel {
    private GridBeanPanel parentPanel;
    //ButtonModel unblurButtonModel;
    /**
     * Add UID
     */
    private static final long serialVersionUID = 2026487374939363876L;

    /** Creates new form PreProcJPanel */
    /**
     * @param parentpanel
     * @throws DataSetException
     */
    public PreProcJPanel(GridBeanPanel parentpanel) throws DataSetException {
        this.parentPanel = parentpanel;
        initComponents();
        initCustomComponents();
        bindComponents();
    }

    /**
     * Method to initialise the binding components 
     */
    private void initCustomComponents() {
        //Create Dir button group
        createDirbuttonGroup.add(yesCreateDirRadioButton);
        createDirbuttonGroup.add(noCreateDirRadioButton);
        //Data organization button group 
        filePreProcDirectivesbuttonGroup.add(lookupDataRadioButton);
        filePreProcDirectivesbuttonGroup.add(distributeRadioButton);
        filePreProcDirectivesbuttonGroup.add(regroupRadioButton);
        filePreProcDirectivesbuttonGroup.add(preProcessRadioButton);
        //Data organization directive button group
        organizeReorNothingbuttonGroup.add(organizeRadioButton);
        organizeReorNothingbuttonGroup.add(reOrganizeRadioButton);
        organizeReorNothingbuttonGroup.add(doNothingRadioButton);
        organizeReorNothingbuttonGroup.add(unorganizeRadioButton);
        //Data extraction directive button group
        dataExtractionbuttonGroup.add(extractDataRadioButton);
        dataExtractionbuttonGroup.add(noDataExtractionRadioButton);
        //Verbosity button group
        verbositybuttonGroup.add(yesVerbosityRadioButton);
        verbositybuttonGroup.add(noVerbosityRadioButton);
        //Delete button group
        deleteDirYesNobuttonGroup.add(yesDeleteDirRadioButton);
        deleteDirYesNobuttonGroup.add(noDeleteDirRadioButton);
    }

    /**
     * @throws DataSetException
     */
    private void bindComponents() throws DataSetException {
        //--simple_data_path=$SIMPLE_DATA_PATH \ #='path to data';
        parentPanel.linkTextField(SimplePreProcGridBean.SIMPLE_DATA_PATH, dataPathTextField);
        parentPanel.setValueTranslator(SimplePreProcGridBean.SIMPLE_DATA_PATH, StringValueTranslator.getInstance());
        //--target_file=$TARGET_FILE \   #= 'target.txt';#Name of output file list
        parentPanel.linkTextField(SimplePreProcGridBean.TARGET_FILE, targetFileTextField);
        parentPanel.setValueTranslator(SimplePreProcGridBean.TARGET_FILE, StringValueTranslator.getInstance());
        //--unblur_dir=$UNBLUR_DIR \     #= "/opt/Devel_tools/unblur_1.0.2/bin";
        parentPanel.linkTextField(SimplePreProcGridBean.UNBLUR_DIR, unblurPathTextField);
        parentPanel.setValueTranslator(SimplePreProcGridBean.UNBLUR_DIR, StringValueTranslator.getInstance());
        //--ctffind_dir=$CTFFIND_DIR \   #= "/opt/Devel_tools/ctffind-4.0.16/build";
        parentPanel.linkTextField(SimplePreProcGridBean.CTFFIND_DIR, ctffindPathTextField);
        parentPanel.setValueTranslator(SimplePreProcGridBean.CTFFIND_DIR, StringValueTranslator.getInstance());
        //--nframes=$NFRAMES \           #= 7;
        parentPanel.linkTextField(SimplePreProcGridBean.NFRAMES, nframesTextField);
        parentPanel.setValueTranslator(SimplePreProcGridBean.NFRAMES, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePreProcGridBean.NFRAMES, new IntegerValueValidator());
        //--nconfig=$NCONFIG \           #= 452; #not used for checking only, real nconfig extracted
        parentPanel.linkTextField(SimplePreProcGridBean.NCONFIG, nconfigsTextField);
        parentPanel.setValueTranslator(SimplePreProcGridBean.NCONFIG, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePreProcGridBean.NCONFIG, new IntegerValueValidator());
        //--sph_abe=$SPH_ABE \           #= 2.7; #sphererical aberration
        parentPanel.linkTextField(SimplePreProcGridBean.SPH_ABE, sphericalAberrationTextField);
        parentPanel.setValueTranslator(SimplePreProcGridBean.SPH_ABE, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePreProcGridBean.SPH_ABE, new DoubleValueValidator());
        //--amp_con=$AMP_CON \           #= 0.07;#Amplitude constrast
        parentPanel.linkTextField(SimplePreProcGridBean.AMP_CON, amplitudeContrastTextField);
        parentPanel.setValueTranslator(SimplePreProcGridBean.AMP_CON, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePreProcGridBean.AMP_CON, new DoubleValueValidator());
        //--sze_pwr_spc=$SZE_PWR_SPC \   #= 512; #Size of prower spectrum
        parentPanel.linkTextField(SimplePreProcGridBean.SZE_PWR_SPC, sizePowerSpectrumTextField);
        parentPanel.setValueTranslator(SimplePreProcGridBean.SZE_PWR_SPC, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePreProcGridBean.SZE_PWR_SPC, new IntegerValueValidator());
        //--data_rot=$DATA_ROT \         #= -90; #Rotation angle for the Size of prower spectrum
        parentPanel.linkTextField(SimplePreProcGridBean.DATA_ROT, rotationAngleTextField);
        parentPanel.setValueTranslator(SimplePreProcGridBean.DATA_ROT, StringValueTranslator.getInstance());
        parentPanel.setValueValidator(SimplePreProcGridBean.DATA_ROT, new IntegerValueValidator());
        //Listeners for the Radio buttons
        //Action listener for the unblur path radio button 
        unblurRadioButton.addActionListener(new ActionListener() {
            @Override public void actionPerformed(ActionEvent e) {
                if (unblurRadioButton.isSelected() ) {
                    EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("1");
                    parentPanel.getGridBeanModel().set(SimplePreProcGridBean.UNBLUR, value);
                } else {
                    EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("0");
                    parentPanel.getGridBeanModel().set(SimplePreProcGridBean.UNBLUR, value);            
                }
            }
        });
        //Action listener for the ctfFind path radio button 
        ctffindPathSetRadioButton.addActionListener(new ActionListener() {
            @Override public void actionPerformed(ActionEvent e) {
                if (ctffindPathSetRadioButton.isSelected() ) {
                    EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("1");
                    parentPanel.getGridBeanModel().set(SimplePreProcGridBean.CTFFIND, value);
                } else {
                    EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("0");
                    parentPanel.getGridBeanModel().set(SimplePreProcGridBean.CTFFIND, value);            
                }
            }
        });

        
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        filePreProcDirectivesbuttonGroup = new javax.swing.ButtonGroup();
        organizeReorNothingbuttonGroup = new javax.swing.ButtonGroup();
        dataExtractionbuttonGroup = new javax.swing.ButtonGroup();
        verbositybuttonGroup = new javax.swing.ButtonGroup();
        deleteDirYesNobuttonGroup = new javax.swing.ButtonGroup();
        createDirbuttonGroup = new javax.swing.ButtonGroup();
        mainSimplePreProcScrollPanel = new javax.swing.JScrollPane();
        mainSimplePreProcPanel = new javax.swing.JPanel();
        pathPanel = new javax.swing.JPanel();
        dataPathTextField = new javax.swing.JTextField();
        browsePathButton = new javax.swing.JButton();
        pathDataLabel = new javax.swing.JLabel();
        microscopedetailsPanel = new javax.swing.JPanel();
        nframesLabel = new javax.swing.JLabel();
        nframesTextField = new javax.swing.JTextField();
        nconfigsLabel = new javax.swing.JLabel();
        nconfigsTextField = new javax.swing.JTextField();
        sphericalAberrationLabel = new javax.swing.JLabel();
        sphericalAberrationTextField = new javax.swing.JTextField();
        amplitudeContrastLabel = new javax.swing.JLabel();
        amplitudeContrastTextField = new javax.swing.JTextField();
        sizePowerSpectrumLabel = new javax.swing.JLabel();
        sizePowerSpectrumTextField = new javax.swing.JTextField();
        rotationAngleLabel = new javax.swing.JLabel();
        rotationAngleTextField = new javax.swing.JTextField();
        targetFilePanel = new javax.swing.JPanel();
        targetFileLabel = new javax.swing.JLabel();
        targetFileTextField = new javax.swing.JTextField();
        viewTargetFileButton = new javax.swing.JButton();
        createDirPanel = new javax.swing.JPanel();
        noCreateDirRadioButton = new javax.swing.JRadioButton();
        yesCreateDirRadioButton = new javax.swing.JRadioButton();
        preProcessingPanel = new javax.swing.JPanel();
        lookupDataRadioButton = new javax.swing.JRadioButton();
        distributeRadioButton = new javax.swing.JRadioButton();
        regroupRadioButton = new javax.swing.JRadioButton();
        preProcessPanel = new javax.swing.JPanel();
        preProcessRadioButton = new javax.swing.JRadioButton();
        organizeRadioButton = new javax.swing.JRadioButton();
        reOrganizeRadioButton = new javax.swing.JRadioButton();
        doNothingRadioButton = new javax.swing.JRadioButton();
        unorganizeRadioButton = new javax.swing.JRadioButton();
        dataExtractionPanel = new javax.swing.JPanel();
        noDataExtractionRadioButton = new javax.swing.JRadioButton();
        extractDataRadioButton = new javax.swing.JRadioButton();
        movieCorrectionPanel = new javax.swing.JPanel();
        unblurRadioButton = new javax.swing.JRadioButton();
        unblurPathTextField = new javax.swing.JTextField();
        unblurPathSetButton = new javax.swing.JButton();
        ctffindPathSetRadioButton = new javax.swing.JRadioButton();
        ctffindPathTextField = new javax.swing.JTextField();
        ctffindPathSetButton = new javax.swing.JButton();
        verbosityPanel = new javax.swing.JPanel();
        yesVerbosityRadioButton = new javax.swing.JRadioButton();
        noVerbosityRadioButton = new javax.swing.JRadioButton();
        dataCleanUpPanel = new javax.swing.JPanel();
        yesDeleteDirRadioButton = new javax.swing.JRadioButton();
        noDeleteDirRadioButton = new javax.swing.JRadioButton();

        setLayout(new javax.swing.BoxLayout(this, javax.swing.BoxLayout.LINE_AXIS));

        mainSimplePreProcPanel.setToolTipText("");

        pathPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Path to the data panel..."));

        browsePathButton.setText("Browse");

        pathDataLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        pathDataLabel.setText("Data:");

        javax.swing.GroupLayout pathPanelLayout = new javax.swing.GroupLayout(pathPanel);
        pathPanel.setLayout(pathPanelLayout);
        pathPanelLayout.setHorizontalGroup(
            pathPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pathPanelLayout.createSequentialGroup()
                .addComponent(pathDataLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(dataPathTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 265, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(browsePathButton)
                .addContainerGap(27, Short.MAX_VALUE))
        );
        pathPanelLayout.setVerticalGroup(
            pathPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pathPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(pathPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(dataPathTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(pathDataLabel)
                    .addComponent(browsePathButton))
                .addContainerGap(14, Short.MAX_VALUE))
        );

        microscopedetailsPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Microscope details..."));

        nframesLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        nframesLabel.setText("nframes:");

        nconfigsLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        nconfigsLabel.setText("nconfigs:");

        sphericalAberrationLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        sphericalAberrationLabel.setText("Spherical Aberration:");

        amplitudeContrastLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        amplitudeContrastLabel.setText("Amp contrast:");

        sizePowerSpectrumLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        sizePowerSpectrumLabel.setText("Size power spectrum:");

        rotationAngleLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        rotationAngleLabel.setText("Rot Ang:");

        javax.swing.GroupLayout microscopedetailsPanelLayout = new javax.swing.GroupLayout(microscopedetailsPanel);
        microscopedetailsPanel.setLayout(microscopedetailsPanelLayout);
        microscopedetailsPanelLayout.setHorizontalGroup(
            microscopedetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(microscopedetailsPanelLayout.createSequentialGroup()
                .addGroup(microscopedetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(microscopedetailsPanelLayout.createSequentialGroup()
                        .addComponent(amplitudeContrastLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(amplitudeContrastTextField))
                    .addGroup(microscopedetailsPanelLayout.createSequentialGroup()
                        .addGroup(microscopedetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                            .addGroup(microscopedetailsPanelLayout.createSequentialGroup()
                                .addComponent(rotationAngleLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(rotationAngleTextField))
                            .addGroup(javax.swing.GroupLayout.Alignment.LEADING, microscopedetailsPanelLayout.createSequentialGroup()
                                .addComponent(nframesLabel)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(nframesTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 31, javax.swing.GroupLayout.PREFERRED_SIZE)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(nconfigsLabel)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(microscopedetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(microscopedetailsPanelLayout.createSequentialGroup()
                        .addComponent(nconfigsTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 47, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(sphericalAberrationLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(sphericalAberrationTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 48, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(microscopedetailsPanelLayout.createSequentialGroup()
                        .addComponent(sizePowerSpectrumLabel)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(sizePowerSpectrumTextField)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        microscopedetailsPanelLayout.setVerticalGroup(
            microscopedetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(microscopedetailsPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(microscopedetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(nframesLabel)
                    .addComponent(nframesTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(nconfigsLabel)
                    .addComponent(nconfigsTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(sphericalAberrationLabel)
                    .addComponent(sphericalAberrationTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(microscopedetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(amplitudeContrastLabel)
                    .addComponent(amplitudeContrastTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(sizePowerSpectrumLabel)
                    .addComponent(sizePowerSpectrumTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(microscopedetailsPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(rotationAngleLabel)
                    .addComponent(rotationAngleTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(14, Short.MAX_VALUE))
        );

        targetFilePanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Target file name..."));

        targetFileLabel.setFont(new java.awt.Font("Dialog", 1, 10));
        targetFileLabel.setText("Target file name:");

        viewTargetFileButton.setText("View");

        createDirPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Create Dir..."));

        noCreateDirRadioButton.setText("No");
        noCreateDirRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noCreateDirRadioButtonActionPerformed(evt);
            }
        });

        yesCreateDirRadioButton.setText("Yes");
        yesCreateDirRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesCreateDirRadioButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout createDirPanelLayout = new javax.swing.GroupLayout(createDirPanel);
        createDirPanel.setLayout(createDirPanelLayout);
        createDirPanelLayout.setHorizontalGroup(
            createDirPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(createDirPanelLayout.createSequentialGroup()
                .addComponent(noCreateDirRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(yesCreateDirRadioButton)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        createDirPanelLayout.setVerticalGroup(
            createDirPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(createDirPanelLayout.createSequentialGroup()
                .addGroup(createDirPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(noCreateDirRadioButton)
                    .addComponent(yesCreateDirRadioButton))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout targetFilePanelLayout = new javax.swing.GroupLayout(targetFilePanel);
        targetFilePanel.setLayout(targetFilePanelLayout);
        targetFilePanelLayout.setHorizontalGroup(
            targetFilePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(targetFilePanelLayout.createSequentialGroup()
                .addComponent(targetFileLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(targetFileTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 133, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(viewTargetFileButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(createDirPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
        );
        targetFilePanelLayout.setVerticalGroup(
            targetFilePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(targetFilePanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(targetFilePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(targetFileLabel)
                    .addComponent(targetFileTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(viewTargetFileButton)))
            .addComponent(createDirPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 50, javax.swing.GroupLayout.PREFERRED_SIZE)
        );

        preProcessingPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Pre-procesing directives..."));

        lookupDataRadioButton.setText("Lookup Data");
        lookupDataRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                lookupDataRadioButtonActionPerformed(evt);
            }
        });

        distributeRadioButton.setText("Distribute");
        distributeRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                distributeRadioButtonActionPerformed(evt);
            }
        });

        regroupRadioButton.setText("Regroup");
        regroupRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                regroupRadioButtonActionPerformed(evt);
            }
        });

        preProcessPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Pre-process data..."));

        preProcessRadioButton.setText("Pre-process");
        preProcessRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                preProcessRadioButtonActionPerformed(evt);
            }
        });

        organizeRadioButton.setText("Organize");
        organizeRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                organizeRadioButtonActionPerformed(evt);
            }
        });

        reOrganizeRadioButton.setText("Re-organize");
        reOrganizeRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                reOrganizeRadioButtonActionPerformed(evt);
            }
        });

        doNothingRadioButton.setText("Do Nothing");
        doNothingRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                doNothingRadioButtonActionPerformed(evt);
            }
        });

        unorganizeRadioButton.setText("Un-organize");
        unorganizeRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                unorganizeRadioButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout preProcessPanelLayout = new javax.swing.GroupLayout(preProcessPanel);
        preProcessPanel.setLayout(preProcessPanelLayout);
        preProcessPanelLayout.setHorizontalGroup(
            preProcessPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(preProcessRadioButton)
            .addGroup(preProcessPanelLayout.createSequentialGroup()
                .addComponent(organizeRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(reOrganizeRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(doNothingRadioButton))
            .addComponent(unorganizeRadioButton)
        );
        preProcessPanelLayout.setVerticalGroup(
            preProcessPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(preProcessPanelLayout.createSequentialGroup()
                .addComponent(preProcessRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(preProcessPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(organizeRadioButton)
                    .addComponent(reOrganizeRadioButton)
                    .addComponent(doNothingRadioButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(unorganizeRadioButton))
        );

        dataExtractionPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Data extraction..."));

        noDataExtractionRadioButton.setText("No");
        noDataExtractionRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noDataExtractionRadioButtonActionPerformed(evt);
            }
        });

        extractDataRadioButton.setText("Extract data");
        extractDataRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                extractDataRadioButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout dataExtractionPanelLayout = new javax.swing.GroupLayout(dataExtractionPanel);
        dataExtractionPanel.setLayout(dataExtractionPanelLayout);
        dataExtractionPanelLayout.setHorizontalGroup(
            dataExtractionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(dataExtractionPanelLayout.createSequentialGroup()
                .addComponent(extractDataRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(noDataExtractionRadioButton)
                .addContainerGap(99, Short.MAX_VALUE))
        );
        dataExtractionPanelLayout.setVerticalGroup(
            dataExtractionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(dataExtractionPanelLayout.createSequentialGroup()
                .addGroup(dataExtractionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(extractDataRadioButton)
                    .addComponent(noDataExtractionRadioButton))
                .addContainerGap(18, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout preProcessingPanelLayout = new javax.swing.GroupLayout(preProcessingPanel);
        preProcessingPanel.setLayout(preProcessingPanelLayout);
        preProcessingPanelLayout.setHorizontalGroup(
            preProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(preProcessingPanelLayout.createSequentialGroup()
                .addGroup(preProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(preProcessingPanelLayout.createSequentialGroup()
                        .addComponent(lookupDataRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(distributeRadioButton)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(regroupRadioButton))
                    .addComponent(dataExtractionPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(preProcessPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        preProcessingPanelLayout.setVerticalGroup(
            preProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(preProcessingPanelLayout.createSequentialGroup()
                .addGroup(preProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(preProcessPanel, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, preProcessingPanelLayout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(preProcessingPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(lookupDataRadioButton)
                            .addComponent(distributeRadioButton)
                            .addComponent(regroupRadioButton))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(dataExtractionPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .addContainerGap())
        );

        movieCorrectionPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Movie correction..."));

        unblurRadioButton.setText("Unblur-->Path:");
        unblurRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                unblurRadioButtonActionPerformed(evt);
            }
        });

        unblurPathSetButton.setText("Unblur Path");

        ctffindPathSetRadioButton.setText("CTFFind-->Path:");
        ctffindPathSetRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ctffindPathSetRadioButtonActionPerformed(evt);
            }
        });

        ctffindPathSetButton.setText("CTFFind4 Path");

        javax.swing.GroupLayout movieCorrectionPanelLayout = new javax.swing.GroupLayout(movieCorrectionPanel);
        movieCorrectionPanel.setLayout(movieCorrectionPanelLayout);
        movieCorrectionPanelLayout.setHorizontalGroup(
            movieCorrectionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(movieCorrectionPanelLayout.createSequentialGroup()
                .addGroup(movieCorrectionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(unblurRadioButton)
                    .addComponent(ctffindPathSetRadioButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(movieCorrectionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(ctffindPathTextField)
                    .addComponent(unblurPathTextField, javax.swing.GroupLayout.DEFAULT_SIZE, 273, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(movieCorrectionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(unblurPathSetButton)
                    .addComponent(ctffindPathSetButton))
                .addContainerGap(51, Short.MAX_VALUE))
        );
        movieCorrectionPanelLayout.setVerticalGroup(
            movieCorrectionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(movieCorrectionPanelLayout.createSequentialGroup()
                .addGroup(movieCorrectionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(unblurPathTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(unblurPathSetButton)
                    .addComponent(unblurRadioButton))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(movieCorrectionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(movieCorrectionPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(ctffindPathSetRadioButton)
                        .addComponent(ctffindPathTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(ctffindPathSetButton))
                .addContainerGap(13, Short.MAX_VALUE))
        );

        verbosityPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Verbosity..."));

        yesVerbosityRadioButton.setText("Yes");
        yesVerbosityRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesVerbosityRadioButtonActionPerformed(evt);
            }
        });

        noVerbosityRadioButton.setText("No");
        noVerbosityRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noVerbosityRadioButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout verbosityPanelLayout = new javax.swing.GroupLayout(verbosityPanel);
        verbosityPanel.setLayout(verbosityPanelLayout);
        verbosityPanelLayout.setHorizontalGroup(
            verbosityPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(verbosityPanelLayout.createSequentialGroup()
                .addComponent(yesVerbosityRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(noVerbosityRadioButton)
                .addContainerGap(75, Short.MAX_VALUE))
        );
        verbosityPanelLayout.setVerticalGroup(
            verbosityPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, verbosityPanelLayout.createSequentialGroup()
                .addGroup(verbosityPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(yesVerbosityRadioButton)
                    .addComponent(noVerbosityRadioButton))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        dataCleanUpPanel.setBorder(javax.swing.BorderFactory.createTitledBorder("Delete Directories..."));

        yesDeleteDirRadioButton.setText("Yes");
        yesDeleteDirRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                yesDeleteDirRadioButtonActionPerformed(evt);
            }
        });

        noDeleteDirRadioButton.setText("No");
        noDeleteDirRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                noDeleteDirRadioButtonActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout dataCleanUpPanelLayout = new javax.swing.GroupLayout(dataCleanUpPanel);
        dataCleanUpPanel.setLayout(dataCleanUpPanelLayout);
        dataCleanUpPanelLayout.setHorizontalGroup(
            dataCleanUpPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(dataCleanUpPanelLayout.createSequentialGroup()
                .addComponent(yesDeleteDirRadioButton)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(noDeleteDirRadioButton)
                .addContainerGap(107, Short.MAX_VALUE))
        );
        dataCleanUpPanelLayout.setVerticalGroup(
            dataCleanUpPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(dataCleanUpPanelLayout.createSequentialGroup()
                .addGroup(dataCleanUpPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(yesDeleteDirRadioButton)
                    .addComponent(noDeleteDirRadioButton))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        javax.swing.GroupLayout mainSimplePreProcPanelLayout = new javax.swing.GroupLayout(mainSimplePreProcPanel);
        mainSimplePreProcPanel.setLayout(mainSimplePreProcPanelLayout);
        mainSimplePreProcPanelLayout.setHorizontalGroup(
            mainSimplePreProcPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mainSimplePreProcPanelLayout.createSequentialGroup()
                .addGroup(mainSimplePreProcPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(mainSimplePreProcPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                        .addComponent(pathPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(microscopedetailsPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(targetFilePanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(preProcessingPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(movieCorrectionPanel, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(mainSimplePreProcPanelLayout.createSequentialGroup()
                        .addComponent(verbosityPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(dataCleanUpPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(51, Short.MAX_VALUE))
        );
        mainSimplePreProcPanelLayout.setVerticalGroup(
            mainSimplePreProcPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(mainSimplePreProcPanelLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(pathPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(microscopedetailsPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(targetFilePanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(preProcessingPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(movieCorrectionPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(mainSimplePreProcPanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(verbosityPanel, javax.swing.GroupLayout.PREFERRED_SIZE, 57, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(dataCleanUpPanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(353, Short.MAX_VALUE))
        );

        mainSimplePreProcScrollPanel.setViewportView(mainSimplePreProcPanel);

        add(mainSimplePreProcScrollPanel);
    }// </editor-fold>//GEN-END:initComponents

    private void lookupDataRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_lookupDataRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("lookupdata");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.FILE_PRE_HANDLER, value);
    }//GEN-LAST:event_lookupDataRadioButtonActionPerformed
    private void distributeRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_distributeRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("distribute");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.FILE_PRE_HANDLER, value);
    }//GEN-LAST:event_distributeRadioButtonActionPerformed
    private void regroupRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_regroupRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("regroup");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.FILE_PRE_HANDLER, value);
    }//GEN-LAST:event_regroupRadioButtonActionPerformed
    private void preProcessRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_preProcessRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("preprocess");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.FILE_PRE_HANDLER, value);
    }//GEN-LAST:event_preProcessRadioButtonActionPerformed
    private void organizeRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_organizeRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("or");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.FILE_PRE_HANDLER, value);
    }//GEN-LAST:event_organizeRadioButtonActionPerformed
    private void reOrganizeRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_reOrganizeRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("re");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.FILE_PRE_HANDLER, value);
    }//GEN-LAST:event_reOrganizeRadioButtonActionPerformed
    private void unorganizeRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_unorganizeRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("un");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.FILE_PRE_HANDLER, value);
    }//GEN-LAST:event_unorganizeRadioButtonActionPerformed
    private void doNothingRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_doNothingRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.FILE_PRE_HANDLER, value);
    }//GEN-LAST:event_doNothingRadioButtonActionPerformed
    private void extractDataRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_extractDataRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("extract_data");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.FILE_PST_HANDLER, value);
    }//GEN-LAST:event_extractDataRadioButtonActionPerformed
    private void noDataExtractionRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noDataExtractionRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("no");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.FILE_PST_HANDLER, value);
    }//GEN-LAST:event_noDataExtractionRadioButtonActionPerformed
    private void unblurRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_unblurRadioButtonActionPerformed
        // TODO: Done actionlistener taken care in initCustomComponenent method
    }//GEN-LAST:event_unblurRadioButtonActionPerformed
    private void ctffindPathSetRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ctffindPathSetRadioButtonActionPerformed
        // TODO: Done actionlistener taken care in initCustomComponenent method
    }//GEN-LAST:event_ctffindPathSetRadioButtonActionPerformed
    private void yesVerbosityRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesVerbosityRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("1");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.VERBOSE, value);
    }//GEN-LAST:event_yesVerbosityRadioButtonActionPerformed
    private void noVerbosityRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noVerbosityRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("0");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.VERBOSE, value);
    }//GEN-LAST:event_noVerbosityRadioButtonActionPerformed
    private void yesDeleteDirRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesDeleteDirRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("1");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.DELETE_DIR, value);
    }//GEN-LAST:event_yesDeleteDirRadioButtonActionPerformed
    private void noDeleteDirRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noDeleteDirRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("0");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.DELETE_DIR, value);
    }//GEN-LAST:event_noDeleteDirRadioButtonActionPerformed
    private void yesCreateDirRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_yesCreateDirRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("1");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.CREATE_DIR, value);
    }//GEN-LAST:event_yesCreateDirRadioButtonActionPerformed
    private void noCreateDirRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_noCreateDirRadioButtonActionPerformed
        EnvironmentVariableParameterValue value = new EnvironmentVariableParameterValue("0");
        parentPanel.getGridBeanModel().set(SimplePreProcGridBean.CREATE_DIR, value);
    }//GEN-LAST:event_noCreateDirRadioButtonActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel amplitudeContrastLabel;
    private javax.swing.JTextField amplitudeContrastTextField;
    private javax.swing.JButton browsePathButton;
    private javax.swing.JPanel createDirPanel;
    private javax.swing.ButtonGroup createDirbuttonGroup;
    private javax.swing.JButton ctffindPathSetButton;
    private javax.swing.JRadioButton ctffindPathSetRadioButton;
    private javax.swing.JTextField ctffindPathTextField;
    private javax.swing.JPanel dataCleanUpPanel;
    private javax.swing.JPanel dataExtractionPanel;
    private javax.swing.ButtonGroup dataExtractionbuttonGroup;
    private javax.swing.JTextField dataPathTextField;
    private javax.swing.ButtonGroup deleteDirYesNobuttonGroup;
    private javax.swing.JRadioButton distributeRadioButton;
    private javax.swing.JRadioButton doNothingRadioButton;
    private javax.swing.JRadioButton extractDataRadioButton;
    private javax.swing.ButtonGroup filePreProcDirectivesbuttonGroup;
    private javax.swing.JRadioButton lookupDataRadioButton;
    private javax.swing.JPanel mainSimplePreProcPanel;
    private javax.swing.JScrollPane mainSimplePreProcScrollPanel;
    private javax.swing.JPanel microscopedetailsPanel;
    private javax.swing.JPanel movieCorrectionPanel;
    private javax.swing.JLabel nconfigsLabel;
    private javax.swing.JTextField nconfigsTextField;
    private javax.swing.JLabel nframesLabel;
    private javax.swing.JTextField nframesTextField;
    private javax.swing.JRadioButton noCreateDirRadioButton;
    private javax.swing.JRadioButton noDataExtractionRadioButton;
    private javax.swing.JRadioButton noDeleteDirRadioButton;
    private javax.swing.JRadioButton noVerbosityRadioButton;
    private javax.swing.JRadioButton organizeRadioButton;
    private javax.swing.ButtonGroup organizeReorNothingbuttonGroup;
    private javax.swing.JLabel pathDataLabel;
    private javax.swing.JPanel pathPanel;
    private javax.swing.JPanel preProcessPanel;
    private javax.swing.JRadioButton preProcessRadioButton;
    private javax.swing.JPanel preProcessingPanel;
    private javax.swing.JRadioButton reOrganizeRadioButton;
    private javax.swing.JRadioButton regroupRadioButton;
    private javax.swing.JLabel rotationAngleLabel;
    private javax.swing.JTextField rotationAngleTextField;
    private javax.swing.JLabel sizePowerSpectrumLabel;
    private javax.swing.JTextField sizePowerSpectrumTextField;
    private javax.swing.JLabel sphericalAberrationLabel;
    private javax.swing.JTextField sphericalAberrationTextField;
    private javax.swing.JLabel targetFileLabel;
    private javax.swing.JPanel targetFilePanel;
    private javax.swing.JTextField targetFileTextField;
    private javax.swing.JButton unblurPathSetButton;
    private javax.swing.JTextField unblurPathTextField;
    private javax.swing.JRadioButton unblurRadioButton;
    private javax.swing.JRadioButton unorganizeRadioButton;
    private javax.swing.JPanel verbosityPanel;
    private javax.swing.ButtonGroup verbositybuttonGroup;
    private javax.swing.JButton viewTargetFileButton;
    private javax.swing.JRadioButton yesCreateDirRadioButton;
    private javax.swing.JRadioButton yesDeleteDirRadioButton;
    private javax.swing.JRadioButton yesVerbosityRadioButton;
    // End of variables declaration//GEN-END:variables

}
