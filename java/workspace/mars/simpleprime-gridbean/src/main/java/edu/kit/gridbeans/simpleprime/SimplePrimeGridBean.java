package edu.kit.gridbeans.simpleprime;

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

public class SimplePrimeGridBean extends AbstractGridBean implements IGridBean 
{

    private final static String QNAME_PREFIX = "SimplePrime";
    
    private static String APPLICATION_NAME = "SimplePrime";
  
    private static String APPLICATION_VERSION = "2.0";
  
    public final static String NAMESPACE = "http://edu.kit/gridbeans/SimplePrime";


    /** Declare input files */
    public final static QName INPUT_FILE = new QName(NAMESPACE, "INPUT_FILE", QNAME_PREFIX);
    public final static QName FILESET_IN   = new QName(NAMESPACE, "FILESET_IN",   QNAME_PREFIX);
    
    /** Declare output files */
    public final static QName OUTPUT_FILE = new QName(NAMESPACE, "OUTPUT_FILE", QNAME_PREFIX);
    public final static QName FILESET_OUT   = new QName(NAMESPACE, "FILESET_OUT",   QNAME_PREFIX);

    /**
     * QName
     */
    public final static QName TITLE          = new QName(NAMESPACE, "TITLE"         , QNAME_PREFIX);
    public final static QName STK            = new QName(NAMESPACE, "STK"           , QNAME_PREFIX);
    public final static QName SMPD           = new QName(NAMESPACE, "SMPD"          , QNAME_PREFIX);
    public final static QName MSK            = new QName(NAMESPACE, "MSK"           , QNAME_PREFIX);
    public final static QName VOL1           = new QName(NAMESPACE, "VOL1"          , QNAME_PREFIX);
    public final static QName VOL2           = new QName(NAMESPACE, "VOL2"          , QNAME_PREFIX);
    public final static QName ORITAB         = new QName(NAMESPACE, "ORITAB"        , QNAME_PREFIX);
    public final static QName TRS            = new QName(NAMESPACE, "TRS"           , QNAME_PREFIX);
    public final static QName LP             = new QName(NAMESPACE, "LP"            , QNAME_PREFIX);
    public final static QName DYNLP          = new QName(NAMESPACE, "DYNLP"         , QNAME_PREFIX);
    public final static QName NSTATES        = new QName(NAMESPACE, "NSTATES"       , QNAME_PREFIX);
    public final static QName FRAC           = new QName(NAMESPACE, "FRAC"          , QNAME_PREFIX);
    public final static QName MW             = new QName(NAMESPACE, "MW"            , QNAME_PREFIX);
    public final static QName NTHR           = new QName(NAMESPACE, "NTHR"          , QNAME_PREFIX);
    public final static QName STARTIT        = new QName(NAMESPACE, "STARTIT"       , QNAME_PREFIX);
    public final static QName REFINE         = new QName(NAMESPACE, "REFINE"        , QNAME_PREFIX);
    public final static QName LPSTOP         = new QName(NAMESPACE, "LPSTOP"        , QNAME_PREFIX);
    public final static QName DEFTAB         = new QName(NAMESPACE, "DEFTAB"        , QNAME_PREFIX);
    public final static QName NSPACE         = new QName(NAMESPACE, "NSPACE"        , QNAME_PREFIX);
    public final static QName EO             = new QName(NAMESPACE, "EO "           , QNAME_PREFIX);
    public final static QName AMSKLP         = new QName(NAMESPACE, "AMSKLP"        , QNAME_PREFIX);
    public final static QName PGRP           = new QName(NAMESPACE, "PGRP"          , QNAME_PREFIX);
    public final static QName CTF            = new QName(NAMESPACE, "CTF"           , QNAME_PREFIX);
    public final static QName KV             = new QName(NAMESPACE, "KV"            , QNAME_PREFIX);
    public final static QName CS             = new QName(NAMESPACE, "CS"            , QNAME_PREFIX);
    public final static QName FRACA          = new QName(NAMESPACE, "FRACA"         , QNAME_PREFIX);
    public final static QName HP             = new QName(NAMESPACE, "HP"            , QNAME_PREFIX);
    public final static QName NNN            = new QName(NAMESPACE, "NNN"           , QNAME_PREFIX);
    public final static QName DIVERSIFY      = new QName(NAMESPACE, "DIVERSIFY"     , QNAME_PREFIX);
    public final static QName MAXITS         = new QName(NAMESPACE, "MAXITS"        , QNAME_PREFIX);
    public final static QName FRACVALID      = new QName(NAMESPACE, "FRACVALID"     , QNAME_PREFIX);
    public final static QName LPVALID        = new QName(NAMESPACE, "LPVALID"       , QNAME_PREFIX);
    public final static QName EDGE           = new QName(NAMESPACE, "EDGE"          , QNAME_PREFIX);
    public final static QName FIND           = new QName(NAMESPACE, "FIND"          , QNAME_PREFIX);
    public final static QName FSTEP          = new QName(NAMESPACE, "FSTEP"         , QNAME_PREFIX);
    public final static QName NOISE          = new QName(NAMESPACE, "NOISE"         , QNAME_PREFIX);
    public final static QName TIME_PER_IMAGE = new QName(NAMESPACE, "TIME_PER_IMAGE", QNAME_PREFIX);
    public final static QName DENS           = new QName(NAMESPACE, "DENS"          , QNAME_PREFIX);
    public final static QName TRSSTEP        = new QName(NAMESPACE, "TRSSTEP"       , QNAME_PREFIX);
    public final static QName OUTFILE        = new QName(NAMESPACE, "OUTFILE"       , QNAME_PREFIX);
    public final static QName NVOX           = new QName(NAMESPACE, "NVOX"          , QNAME_PREFIX);
    public final static QName INNER          = new QName(NAMESPACE, "INNER"         , QNAME_PREFIX);
    public final static QName WIDTH          = new QName(NAMESPACE, "WIDTH"         , QNAME_PREFIX);
    public final static QName NOREC          = new QName(NAMESPACE, "NOREC"         , QNAME_PREFIX);

    /** Constructor */
    public SimplePrimeGridBean() 
    {
        set(JOBNAME, "SimplePrime");
        createInputEnvironmentVariables();
        createInputFiles();
        createOutputFiles();
    }

    /**
     * create input environment variables (appears in the Variables tab of the GridBean)
     * --workingDirectory=$WORKING_DIR \
     *--stk=$STK \                       #stk=<stack.ext>
     *--smpd=$SMPD \                     #smpd=<sampling distance(in A)>
     *--msk=$MSK \                       #msk=<mask radius(in pixels)>
     *--vol1=$VOL1 \                     #[vol1=<invol.ext>]
     *--vol2=$VOL2 \                     #[vol2=<refvol_2.ext> etc.]
     *--oritab=$ORITAB \                 #[oritab=<previous alignment doc>]
     *--deftab=$DEFTAB \                 #[deftab=<text file defocus values>]
     *--outfile=$OUTFILE \               #[outfile=<output alignment doc 4 parallell jobs>]
     *--kv=$KV \                         #[kv=<acceleration voltage(in kV){300.}>]
     *--cs=$CS \                         #[cs=<spherical aberration constant(in mm){2.7}>]
     *--fraca=$FRACA \                   #[fraca=<frac amp contrast{0.07}>]
     *--frac=$FRAC \                     #[frac=<fraction of ptcls to include{1}>]
     *--pgrp=$PGRP \                     #[pgrp=<cn|dn|t|o|i{c1}>]
     *--ctf=$CTF \                       #[ctf=<yes|no|flip|mul{no}>]
     *--hp=$HP \                         #[hp=<high-pass limit(in A)>]
     *--lp=$LP \                         #[lp=<low-pass limit{20}>]
     *--refine=$REFINE \                 #[refine=<no|exhaust|shift|neigh{no}>]
     *--dynlp=$DYNLP \                   #[dynlp=<yes|no{yes}>]
     *--noise=$NOISE \                   #[noise=<yes|no{no}>]
     *--diversify=$DIVERSIFY \           #[diversify=<yes|no>]
     *--eo=$EO \                         #[eo=<yes|no{no}>]
     *--norec=$NOREC                     #[norec=<yes|no{no}>]
     *--mw=$MW \                         #[mw=<molecular weight(in kD)>]
     *--nstates=$NSTATES \               #[nstates=<nstates to reconstruct>]
     *--startit=$STARTIT \               #[startit=<start iteration>]
     *--trs=$TRS \                       #[trs=<origin shift(in pixels){0}>]
     *--lpstop=$LPSTOP \                 #[lpstop=<stay at this low-pass limit (in A)>]
     *--nspace=$NSPACE \                 #[nspace=<nr reference sections{1000}>]
     *--amsklp=$AMSKLP \                 #[amsklp=<automask low-pass limit(in A)>] 
     *--nnn=$NNN \                       #[nnn=<nr nearest neighbors{300}>]
     *--maxits=$MAXITS \                 #[maxits=<max iterations{100}>]
     *--lpvalid=$LPVALID \               #[lpvalid=<low-pass limit validptcls{20}>]
     *--find=$FIND \                     #[find=<Fourier index>]
     *--fstep=$FSTEP \                   #[fstep=<Fourier step size>]
     *--fracvalid=$FRACVALID \           #[fracvalid=<fraction of particles 4 validation{0.}>]
     *--edge=$EDGE \                     #[edge=<edge size for softening molecular envelope(in pixels){3}>] 
     *--dens=$DENS \                     #[dens=<density(e.g. 9.368 Da/A3 4 gold clusters){0.}>]
     *--trsstep=$TRSSTEP \               #[trsstep=<origin shift stepsize{1}>] 
     *--nvox=$NVOX \                     #[nvox=<nr of voxels in mask{0}>] 
     *--inner=$INNER \                   #[inner=<inner mask radius(in pixels)>]
     *--time_per_image=$TIME_PER_IMAGE \ #[time_per_image=<{100}>] 
     *--width=$WIDTH \                   #[width=<pixels falloff inner mask{10}>]
     *--nthr=$NTHR                       #[nthr=<nr of OpenMP threads{1}>]
     *
     */
    private void createInputEnvironmentVariables() {
        QName[] envVariables = { TITLE,
                                 STK,            //#stk=<stack.ext>
                                 SMPD,           //#smpd=<sampling distance(in A)>
                                 MSK,            //#msk=<mask radius(in pixels)>
                                 VOL1,           //#[vol1=<invol.ext>]
                                 VOL2,           //#[vol2=<refvol_2.ext> etc.]
                                 ORITAB,         //#[oritab=<previous alignment doc>]
                                 TRS,            //#[trs=<origin shift(in pixels){0}>]
                                 LP,             //#[lp=<low-pass limit{20}>]
                                 DYNLP,          //#[dynlp=<yes|no{yes}>]
                                 NSTATES,        //#[nstates=<nstates to reconstruct>]
                                 FRAC,           //#[frac=<fraction of ptcls to include{1}>]
                                 MW,             //#[mw=<molecular weight(in kD)>]
                                 NTHR,           //#[nthr=<nr of OpenMP threads{1}>]
                                 STARTIT,        //#[startit=<start iteration>]
                                 REFINE,         //#[refine=<no|exhaust|shift|neigh{no}>]
                                 LPSTOP,         //#[lpstop=<stay at this low-pass limit (in A)>]
                                 DEFTAB,         //#[deftab=<text file defocus values>]
                                 NSPACE,         //#[nspace=<nr reference sections{1000}>]
                                 EO,             //#[eo=<yes|no{no}>]
                                 AMSKLP,         //#[amsklp=<automask low-pass limit(in A)>] 
                                 PGRP,           //#[pgrp=<cn|dn|t|o|i{c1}>]
                                 CTF,            //#[ctf=<yes|no|flip|mul{no}>]
                                 KV,             //#[kv=<acceleration voltage(in kV){300.}>]
                                 CS,             //#[cs=<spherical aberration constant(in mm){2.7}>]
                                 FRACA,          //#[fraca=<frac amp contrast{0.07}>]
                                 HP,             //#[hp=<high-pass limit(in A)>]
                                 NNN,            //#[nnn=<nr nearest neighbors{300}>]
                                 DIVERSIFY,      //#[diversify=<yes|no>]
                                 MAXITS,         //#[maxits=<max iterations{100}>]
                                 FRACVALID,      //#[fracvalid=<fraction of particles 4 validation{0.}>]
                                 LPVALID,        //#[lpvalid=<low-pass limit validptcls{20}>]
                                 EDGE,           //#[edge=<edge size for softening molecular envelope(in pixels){3}>] 
                                 FIND,           //#[find=<Fourier index>]
                                 FSTEP,          //#[fstep=<Fourier step size>]
                                 NOISE,          //#[noise=<yes|no{no}>]
                                 TIME_PER_IMAGE, //#[time_per_image=<{100}>] 
                                 DENS,           //#[dens=<density(e.g. 9.368 Da/A3 4 gold clusters){0.}>]
                                 TRSSTEP,        //#[trsstep=<origin shift stepsize{1}>] 
                                 OUTFILE,        //#[outfile=<output alignment doc 4 parallell jobs>]
                                 NVOX,           //#[nvox=<nr of voxels in mask{0}>] 
                                 INNER,          //#[inner=<inner mask radius(in pixels)>]
                                 WIDTH,          //#[width=<pixels falloff inner mask{10}>]
                                 NOREC };        //#[norec=<yes|no{no}>]
        String[] initialValues = { "Title",
                                   "sumstack.mrc",
                                   "1.101",
                                   "0.0",
                                   "invol.mrc",
                                   "refvol_2.mrc",
                                   "alignment_doc.txt",
                                   "0",
                                   "20",
                                   "yes",
                                   "1",
                                   "1",
                                   "0.0",
                                   "1",
                                   "1",
                                   "no",
                                   "7.0",
                                   "Deftabs_full_stack.asc",
                                   "1000",
                                   "no",
                                   "20.0",
                                   "c1",
                                   "no",
                                   "300.",
                                   "2.7",
                                   "0.07",
                                   "10.0",
                                   "300",
                                   "no",
                                   "100",
                                   "0.",
                                   "20",
                                   "3", 
                                   "1",
                                   "1",
                                   "no",
                                   "100", 
                                   "0.",
                                   "1", 
                                   "outfile.txt",
                                   "0", 
                                   "0.0",
                                   "10.0",
                                   "no"  };

        
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
//    set(INPUT_FILE, new InputFileParameterValue("Input_File"));
//    
//    QName[] qnameArray = new QName[]{INPUT_FILE};
//    getInputParameters().addAll(ParameterUtils.createFileParameters(qnameArray, true));
//
//    set(FILESET_IN, new FileSetParameterValue());
//	getInputParameters().add(new GridBeanParameter(FILESET_IN, GridBeanParameterType.FILE_SET, true));
  }

  private void createOutputFiles() 
  {
//    set(OUTPUT_FILE, new OutputFileParameterValue("Output_File"));
//    QName[] qnameArray = new QName[]{OUTPUT_FILE};
//    getOutputParameters().addAll(ParameterUtils.createFileParameters(qnameArray, false));
//
//    set(FILESET_OUT, new FileSetParameterValue());
//    getOutputParameters().add(new GridBeanParameter(FILESET_OUT, GridBeanParameterType.FILE_SET, false));
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


  public String getName() {return "SimplePrime";}
  public URL getIconURL() {return getIconURL(SimplePrimeGridBean.class, "images/Ribo_single_cropped.png");}


}
