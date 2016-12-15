package edu.kit.gridbeans.simpleeorecvol;

import java.net.URL;
import java.util.List;
import javax.xml.namespace.QName;
import com.intel.gpe.clients.api.Job;
import com.intel.gpe.clients.api.jsdl.gpe.GPEJob;
import com.intel.gpe.gridbeans.AbstractGridBean;
import com.intel.gpe.gridbeans.GPEConstants;
import com.intel.gpe.gridbeans.GridBeanException;
import com.intel.gpe.gridbeans.IGridBean;
import com.intel.gpe.gridbeans.parameters.GridBeanParameter;
import com.intel.gpe.gridbeans.parameters.GridBeanParameterType;
import com.intel.gpe.gridbeans.parameters.IGridBeanParameterValue;
import com.intel.gpe.gridbeans.parameters.InputFileParameterValue;
import com.intel.gpe.gridbeans.parameters.OutputFileParameterValue;
import com.intel.gpe.gridbeans.parameters.FileSetParameterValue;
import com.intel.gpe.gridbeans.parameters.processing.ParameterUtils;

public class SimpleEoRecVolGridBean extends AbstractGridBean implements IGridBean 
{
    /**
     * Add UID
     */
    private static final long serialVersionUID = 6426140630505599268L;

    private final static String QNAME_PREFIX = "SimpleEoRecVol";
    /** Application name */
    private static String APPLICATION_NAME = "SimpleEoRecVol";
    /** Application version */
    private static String APPLICATION_VERSION = "2.0";
    /** Project NAMESPACE */
    public final static String NAMESPACE = "http://edu.kit/gridbeans/SimpleEoRecVol";

    /** Declare input files */
    public final static QName INPUT_FILE = new QName(NAMESPACE, "INPUT_FILE", QNAME_PREFIX);
    public final static QName FILESET_IN   = new QName(NAMESPACE, "FILESET_IN",   QNAME_PREFIX);

    /** Declare output files */
    public final static QName OUTPUT_FILE = new QName(NAMESPACE, "OUTPUT_FILE", QNAME_PREFIX);
    public final static QName FILESET_OUT   = new QName(NAMESPACE, "FILESET_OUT",   QNAME_PREFIX);

    /** QName for eo_RecVol */
    public final static QName TITLE  = new QName(NAMESPACE, "TITLE" , QNAME_PREFIX);
    public final static QName STK    = new QName(NAMESPACE, "STK"   , QNAME_PREFIX);
    public final static QName MSK    = new QName(NAMESPACE, "MSK"   , QNAME_PREFIX);
    public final static QName SMPD   = new QName(NAMESPACE, "SMPD"  , QNAME_PREFIX);
    public final static QName ORITAB = new QName(NAMESPACE, "ORITAB", QNAME_PREFIX);
    public final static QName FRAC   = new QName(NAMESPACE, "FRAC"  , QNAME_PREFIX);
    public final static QName NTHR   = new QName(NAMESPACE, "NTHR"  , QNAME_PREFIX);
    public final static QName PGRP   = new QName(NAMESPACE, "PGRP"  , QNAME_PREFIX);
    public final static QName MUL    = new QName(NAMESPACE, "MUL"   , QNAME_PREFIX);
    public final static QName PART   = new QName(NAMESPACE, "PART"  , QNAME_PREFIX);
    public final static QName FROMP  = new QName(NAMESPACE, "FROMP" , QNAME_PREFIX);
    public final static QName TOP    = new QName(NAMESPACE, "TOP"   , QNAME_PREFIX);
    public final static QName CTF    = new QName(NAMESPACE, "CTF"   , QNAME_PREFIX);
    public final static QName KV     = new QName(NAMESPACE, "KV"    , QNAME_PREFIX);
    public final static QName FRACA  = new QName(NAMESPACE, "FRACA" , QNAME_PREFIX);
    public final static QName CS     = new QName(NAMESPACE, "CS"    , QNAME_PREFIX);
    public final static QName DEFTAB = new QName(NAMESPACE, "DEFTAB", QNAME_PREFIX);
    public final static QName STATE  = new QName(NAMESPACE, "STATE" , QNAME_PREFIX);
    public final static QName ALPHA  = new QName(NAMESPACE, "ALPHA" , QNAME_PREFIX);
    public final static QName MW     = new QName(NAMESPACE, "MW"    , QNAME_PREFIX);
    public final static QName AMSKLP = new QName(NAMESPACE, "AMSKLP", QNAME_PREFIX);
    public final static QName EDGE   = new QName(NAMESPACE, "EDGE"  , QNAME_PREFIX);
    public final static QName DENS   = new QName(NAMESPACE, "DENS"  , QNAME_PREFIX);
    public final static QName BOX    = new QName(NAMESPACE, "BOX"   , QNAME_PREFIX);
    public final static QName INNER  = new QName(NAMESPACE, "INNER" , QNAME_PREFIX);
    public final static QName WIDTH  = new QName(NAMESPACE, "WIDTH" , QNAME_PREFIX);

  /** Constructor */
  public SimpleEoRecVolGridBean() 
  {
    set(JOBNAME, APPLICATION_NAME);
    createInputEnvironmentVariables();
    createInputFiles();
    createOutputFiles();
  }

  
    /**
     * create input environment variables (appears in the Variables tab of the GridBean)
     * --workingDirectory=$WORKING_DIR \ #working directory local on server
     *--stk=$STK \                      # stk=<ptcls.ext>
     *--oritab=$ORITAB \                # oritab=<algndoc.txt>
     *--smpd=$SMPD \                    # smpd=<sampling distance(in A)>
     *--msk=$MSK \                      # msk=<mask radius(in pixels)>
     *--kv=$KV \                        # [kv=<acceleration voltage(in kV){300.}>]
     *--frac=$FRAC \                    # [frac=<fraction of ptcls to include{1.}>]
     *--pgrp=$PGRP \                    # [pgrp=<cn|dn|t|o|i{c1}>]
     *--mul=$MUL \                      # [mul=<shift multiplication factor{1}>]
     *--fromp=$FROMP \                  # [fromp=<start ptcl{1)>]
     *--top=$TOP \                      # [top=<stop ptcl{nptcls}>]
     *--ctf=$CTF \                      # [ctf=<yes|no|flip|mul{no}>]
     *--fraca=$FRACA \                  # [fraca=<frac amp contrast{0.07}>]
     *--cs=$CS \                        # [cs=<spherical aberration constant(in mm){2.7}>]
     *--deftab=$DEFTAB \                # [deftab=<text file defocus values>]
     *--state=$STATE \                  # [state=<state to reconstruct{all}>]
     *--alpha=$ALPHA \                  # [alpha=<oversampling ratio or padding factor{2.}>]
     *--mw=$MW \                        # [mw=<molecular weight(in kD)>]
     *--amsklp=$AMSKLP \                # [amsklp=<low-pass limit(in A){20}>]
     *--edge=$EDGE \                    # [edge=<edge size for softening molecular envelope(in pixels){3}>]
     *--dens=$DENS \                    # [dens=<density(e.g.9.368 Da/A3 4 gold clusters){0.}>]
     *--box=$BOX \                      # [box=<image size(in pixels)>]
     *--inner=$INNER \                  # [inner=<inner mask radius(in pixels)>]
     *--width=$WIDTH \                  # [width=<pixels falloff inner mask{10}>]
     *--nthr=$NTHR \                    # [nthr=<nr of openMP threads{1}>]
     *--part=$PART                      # [part=<partition number>]
     *
     */
    private void createInputEnvironmentVariables() {
        QName[] envVariables = { TITLE,
                                 STK,        //# stk=<ptcls.ext>
                                 MSK,        //# msk=<mask radius(in pixels)>
                                 SMPD,       //# smpd=<sampling distance(in A)>
                                 ORITAB,     //# oritab=<algndoc.txt>
                                 FRAC,       //# [frac=<fraction of ptcls to include{1.}>]
                                 NTHR,       //# [nthr=<nr of openMP threads{1}>]
                                 PGRP,       //# [pgrp=<cn|dn|t|o|i{c1}>]
                                 MUL,        //# [mul=<shift multiplication factor{1}>]
                                 PART,       //# [part=<partition number>]
                                 FROMP,      //# [fromp=<start ptcl{1)>]
                                 TOP,        //# [top=<stop ptcl{nptcls}>]
                                 CTF,        //# [ctf=<yes|no|flip|mul{no}>]
                                 KV,         //# [kv=<acceleration voltage(in kV){300.}>]
                                 FRACA,      //# [fraca=<frac amp contrast{0.07}>]
                                 CS,         //# [cs=<spherical aberration constant(in mm){2.7}>]
                                 DEFTAB,     //# [deftab=<text file defocus values>]
                                 STATE,      //# [state=<state to reconstruct{all}>]
                                 ALPHA,      //# [alpha=<oversampling ratio or padding factor{2.}>]
                                 MW,         //# [mw=<molecular weight(in kD)>]
                                 AMSKLP,     //# [amsklp=<low-pass limit(in A){20}>]
                                 EDGE,       //# [edge=<edge size for softening molecular envelope(in pixels){3}>]
                                 DENS,       //# [dens=<density(e.g.9.368 Da/A3 4 gold clusters){0.}>]
                                 BOX,        //# [box=<image size(in pixels)>]
                                 INNER,      //# [inner=<inner mask radius(in pixels)>]
                                 WIDTH};     //# [width=<pixels falloff inner mask{10}>]

        String[] initialValues = { "Title",
                                   "sumstack.mrc",
                                   "50",
                                   "1.101",
                                   "algndoc.txt",
                                   "0.7",
                                   "1",
                                   "c1",
                                   "1",
                                   "1",
                                   "1",
                                   "500",
                                   "no",
                                   "300.0",
                                   "0.07",
                                   "2.7",        //# [cs=<spherical aberration constant(in mm){2.7}>]
                                   "deftab.txt", //# [deftab=<text file defocus values>]
                                   "1",        //# [state=<state to reconstruct{all}>]
                                   "2.0",        //# [alpha=<oversampling ratio or padding factor{2.}>]
                                   "5.0",        //# [mw=<molecular weight(in kD)>]
                                   "20",         //# [amsklp=<low-pass limit(in A){20}>]
                                   "3",          //# [edge=<edge size for softening molecular envelope(in pixels){3}>]
                                   "0.0",        //# [dens=<density(e.g.9.368 Da/A3 4 gold clusters){0.}>]
                                   "256",        //# [box=<image size(in pixels)>]
                                   "128",        //# [inner=<inner mask radius(in pixels)>]
                                   "10"};        //# [width=<pixels falloff inner mask{10}>]

        getInputParameters().addAll(ParameterUtils.createEnvParameters(envVariables));

        List<IGridBeanParameterValue> values = ParameterUtils.createEnvParameterValues(envVariables, initialValues);

        // initialize input environment variables
        for (int i = 0; i < initialValues.length; i++) {
            set(envVariables[i], values.get(i));
        }

    }
  
    /** Note: in createInputFiles and createOutputFiles, the InputFileParameterValues and
        OutputFileParameterValues are the specific names of the actual files as they appear
        in the working directory.
    */

    private void createInputFiles() 
    {
//        set(INPUT_FILE, new InputFileParameterValue("Input_File"));
//
//        QName[] qnameArray = new QName[]{INPUT_FILE};
//        getInputParameters().addAll(ParameterUtils.createFileParameters(qnameArray, true));
//
//        set(FILESET_IN, new FileSetParameterValue());
//        getInputParameters().add(new GridBeanParameter(FILESET_IN, GridBeanParameterType.FILE_SET, true));
    }

    private void createOutputFiles() 
    {
//        set(OUTPUT_FILE, new OutputFileParameterValue("Output_File"));
//        QName[] qnameArray = new QName[]{OUTPUT_FILE};
//        getOutputParameters().addAll(ParameterUtils.createFileParameters(qnameArray, false));
//
//        set(FILESET_OUT, new FileSetParameterValue());
//        getOutputParameters().add(new GridBeanParameter(FILESET_OUT, GridBeanParameterType.FILE_SET, false));
    }


    /** Standard code for setting up job definition */
    public void setupJobDefinition(Job job) throws GridBeanException 
    {
        super.setupJobDefinition(job);
        if (job instanceof GPEJob) 
            {
                GPEJob gpeJob = (GPEJob) job;
                gpeJob.setApplicationName(APPLICATION_NAME);
                gpeJob.setApplicationVersion(APPLICATION_VERSION);
                gpeJob.setWorkingDirectory(GPEConstants.JobManagement.TEMPORARY_DIR_NAME);
            }
        else
            {
                throw new GridBeanException("Unsupported job class: " + job.getClass().getName());
            }
    }

    public String getName() {return APPLICATION_NAME;}
    public URL getIconURL() {return getIconURL(SimpleEoRecVolGridBean.class, "images/StarFish_helix_cropped.png");}

}
