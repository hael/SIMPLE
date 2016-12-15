package edu.kit.gridbeans.simplestackops;

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



public class SimpleStackopsGridBean extends AbstractGridBean implements IGridBean 
{

    /**
     * Add UID
     */
    private static final long serialVersionUID = -303699990439850682L;

    private final static String QNAME_PREFIX = "SimpleStackops";
    /** Application name */
    private static String APPLICATION_NAME = "SimpleStackops";
    /** Application version */
    private static String APPLICATION_VERSION = "2.0";
    /** Project NAMESPACE */
    public final static String NAMESPACE = "http://edu.kit/gridbeans/SimpleStackops";

    /** Declare input files */
    public final static QName INPUT_FILE = new QName(NAMESPACE, "INPUT_FILE", QNAME_PREFIX);
    public final static QName FILESET_IN = new QName(NAMESPACE, "FILESET_IN",   QNAME_PREFIX);

    /** Declare output files */
    public final static QName OUTPUT_FILE = new QName(NAMESPACE, "OUTPUT_FILE", QNAME_PREFIX);
    public final static QName FILESET_OUT = new QName(NAMESPACE, "FILESET_OUT",   QNAME_PREFIX);

    /**
     * QName
     */
    public final static QName TITLE      = new QName(NAMESPACE, "TITLE"     , QNAME_PREFIX);
    public final static QName STK        = new QName(NAMESPACE, "STK"       , QNAME_PREFIX);
    public final static QName STK2       = new QName(NAMESPACE, "STK2"      , QNAME_PREFIX);
    public final static QName ORITAB     = new QName(NAMESPACE, "ORITAB"    , QNAME_PREFIX);
    public final static QName OUTSTK     = new QName(NAMESPACE, "OUTSTK"    , QNAME_PREFIX);
    public final static QName DEFTAB     = new QName(NAMESPACE, "DEFTAB"    , QNAME_PREFIX);
    public final static QName OUTFILE    = new QName(NAMESPACE, "OUTFILE"   , QNAME_PREFIX);
    public final static QName SHALGN     = new QName(NAMESPACE, "SHALGN"    , QNAME_PREFIX);
    public final static QName MERGE      = new QName(NAMESPACE, "MERGE"     , QNAME_PREFIX);
    public final static QName ROALGN     = new QName(NAMESPACE, "ROALGN"    , QNAME_PREFIX);
    public final static QName MASSCEN    = new QName(NAMESPACE, "MASSCEN"   , QNAME_PREFIX);
    public final static QName VIS        = new QName(NAMESPACE, "VIS"       , QNAME_PREFIX);
    public final static QName BIN        = new QName(NAMESPACE, "BIN"       , QNAME_PREFIX);
    public final static QName ACF        = new QName(NAMESPACE, "ACF"       , QNAME_PREFIX);
    public final static QName PHRAND     = new QName(NAMESPACE, "PHRAND"    , QNAME_PREFIX);
    public final static QName HFUN       = new QName(NAMESPACE, "HFUN"      , QNAME_PREFIX);
    public final static QName NORM       = new QName(NAMESPACE, "NORM"      , QNAME_PREFIX);
    public final static QName NOISE_NORM = new QName(NAMESPACE, "NOISE_NORM", QNAME_PREFIX);
    public final static QName AVG        = new QName(NAMESPACE, "AVG"       , QNAME_PREFIX);
    public final static QName RANKIFY    = new QName(NAMESPACE, "RANKIFY"   , QNAME_PREFIX);
    public final static QName FILETAB    = new QName(NAMESPACE, "FILETAB"   , QNAME_PREFIX);
    public final static QName STATS      = new QName(NAMESPACE, "STATS"     , QNAME_PREFIX);
    public final static QName CTF        = new QName(NAMESPACE, "CTF"       , QNAME_PREFIX);
    public final static QName FT2IMG     = new QName(NAMESPACE, "FT2IMG"    , QNAME_PREFIX);
    public final static QName COMPARE    = new QName(NAMESPACE, "COMPARE"   , QNAME_PREFIX);
    public final static QName NEG        = new QName(NAMESPACE, "NEG"       , QNAME_PREFIX);
    public final static QName CTFSQ      = new QName(NAMESPACE, "CTFSQ"     , QNAME_PREFIX);
    public final static QName NPTCLS     = new QName(NAMESPACE, "NPTCLS"    , QNAME_PREFIX);
    public final static QName FROMP      = new QName(NAMESPACE, "FROMP"     , QNAME_PREFIX);
    public final static QName TRS        = new QName(NAMESPACE, "TRS"       , QNAME_PREFIX);
    public final static QName STATE      = new QName(NAMESPACE, "STATE"     , QNAME_PREFIX);
    public final static QName SYM_CLASS  = new QName(NAMESPACE, "SYM_CLASS" , QNAME_PREFIX);
    public final static QName TOP        = new QName(NAMESPACE, "TOP"       , QNAME_PREFIX);
    public final static QName NRAN       = new QName(NAMESPACE, "NRAN"      , QNAME_PREFIX);
    public final static QName NEWBOX     = new QName(NAMESPACE, "NEWBOX"    , QNAME_PREFIX);
    public final static QName FRAMEAVG   = new QName(NAMESPACE, "FRAMEAVG"  , QNAME_PREFIX);
    public final static QName CLIP       = new QName(NAMESPACE, "CLIP"      , QNAME_PREFIX);
    public final static QName BOX        = new QName(NAMESPACE, "BOX"       , QNAME_PREFIX);
    public final static QName SPLIT      = new QName(NAMESPACE, "SPLIT"     , QNAME_PREFIX);
    public final static QName NTHR       = new QName(NAMESPACE, "NTHR"      , QNAME_PREFIX);
    public final static QName SMPD       = new QName(NAMESPACE, "SMPD"      , QNAME_PREFIX);
    public final static QName HP         = new QName(NAMESPACE, "HP"        , QNAME_PREFIX);
    public final static QName MUL        = new QName(NAMESPACE, "MUL"       , QNAME_PREFIX);
    public final static QName LP         = new QName(NAMESPACE, "LP"        , QNAME_PREFIX);
    public final static QName FRAC       = new QName(NAMESPACE, "FRAC"      , QNAME_PREFIX);
    public final static QName SNR        = new QName(NAMESPACE, "SNR"       , QNAME_PREFIX);
    public final static QName MSK        = new QName(NAMESPACE, "MSK"       , QNAME_PREFIX);
    public final static QName KV         = new QName(NAMESPACE, "KV"        , QNAME_PREFIX);
    public final static QName SCALE      = new QName(NAMESPACE, "SCALE"     , QNAME_PREFIX);
    public final static QName FRACZERO   = new QName(NAMESPACE, "FRACZERO"  , QNAME_PREFIX);
    public final static QName FRACA      = new QName(NAMESPACE, "FRACA"     , QNAME_PREFIX);
    public final static QName CS         = new QName(NAMESPACE, "CS"        , QNAME_PREFIX);
    public final static QName TRES       = new QName(NAMESPACE, "TRES"      , QNAME_PREFIX);
    public final static QName INNER      = new QName(NAMESPACE, "INNER"     , QNAME_PREFIX);
    public final static QName WIDTH      = new QName(NAMESPACE, "WIDTH"     , QNAME_PREFIX);

    /** Constructor */
    public SimpleStackopsGridBean() 
    {
        set(JOBNAME, APPLICATION_NAME);
        createInputEnvironmentVariables();
        createInputFiles();
        createOutputFiles();
    }
  
    /**
     * create input environment variables (appears in the Variables tab of the GridBean)
     * --workingDirectory=$WORKING_DIR \
     *--stk=$STK \               #[stk=<stack.ext>]
     *--stk2=$STK2 \             #[stk2=<stack2.ext>]
     *--oritab=$ORITAB \         #[oritab=<SIMPLE alignment doc>]
     *--outstk=$OUTSTK \         #[outstk=<outstk.ext>]
     *--deftab=$DEFTAB \         #[deftab=<text file with defocus values>]
     *--filetab=$FILETAB \       #[filetab=<filenames.txt>]
     *--smpd=$SMPD \             #[smpd=<sampling distance(in A)>]
     *--kv=$KV \                 #[kv=<acceleration voltage(in kV){300.}>]
     *--cs=$CS \                 #[cs=<spherical aberration constant(in mm){2.7}>]
     *--fraca=$FRACA  \          #[fraca=<frac amp contrast{0.07}>] 
     *--nptcls=$NPTCLS \         #[nptcls=<nr of imgs>]
     *--frac=$FRAC \             #[frac=<fraction of ptcls to extract{1}>]
     *--sym_class=$SYM_CLASS \   #[class=<symmetry class>]
     *--ctf=$CTF \               #[ctf=<yes|no|flip|mul|abs{no}>]
     *--ctfsq=$CTFSQ \           #[ctfsq=<yes|no{no}>] 
     *--shalgn=$SHALGN \         #[shalgn=<yes|no{no}>]
     *--roalgn=$ROALGN  \        #[roalgn=<yes|no{no}>]
     *--masscen=$MASSCEN \       #[masscen=<yes|no{no}>]
     *--acf=$ACF \               #[acf=<yes|no{no}>]
     *--phrand=$PHRAND \         #[phrand=<yes|no{no}>]
     *--vis=$VIS \               #[vis=<yes|no>]
     *--hp=$HP \                 #[hp=<high-pass limit(in A)>]
     *--fromp=$FROMP \           #[fromp=<start ptcl>]
     *--mul=$MUL \               #[mul=<shift multiplication factor{1}>]
     *--trs=$TRS \               #[trs=<origin]
     *--lp=$LP \                 #[lp=<low-pass limitshift(in pixels){3}>]
     *--snr=$SNR \               #[snr=<signal2noise ratio>]
     *--state=$STATE \           #[state=<state to extract>]
     *--msk=$MSK  \              #[msk=<mask radius(in pixels){box/2}>]
     *--bin=$BIN \               #[bin=<binarize{no}>]
     *--top=$TOP \               #[top=<stop ptcl>]
     *--nran=$NRAN \             #[nran=<number of random images to select>]
     *--newbox=$NEWBOX \         #[newbox=<scaled box>]
     *--scale=$SCALE \           #[scale=<scale factor{1}>]
     *--hfun=$HFUN \             #[hfun=<sigm|tanh|lin{sigm}>]
     *--norm=$NORM \             #[norm=<yes|no{no}>]
     *--noise_norm=$NOISE_NORM \ #[noise_norm=<yes|no>]
     *--avg=$AVG \               #[avg=<yes|no>]
     *--rankify=$RANKIFY \       #[rankify=<yes|no>]
     *--stats=$STATS \           #[stats=<yes|no{yes}>]
     *--compare=$COMPARE \       #[compare=<yes|no{no}>]
     *--neg=$NEG \               #[neg=<yes|no{no}>]
     *--merge=$MERGE \           #[merge=<yes|no{no}>]
     *--ft2img=$FT2IMG \         #[ft2img=<yes|no{no}>]
     *--frameavg=$FRAMEAVG \     #[frameavg=<nr of frames to average{0}>]
     *--clip=$CLIP \             #[clip=<clipped box size{box}>] 
     *--box=$BOX \               #[box=<image size(in pixels)>]
     *--inner=$INNER \           #[inner=<inner mask radius(in pixels)>]
     *--width=$WIDTH \           #[width=<pixels falloff inner mask{10}>]
     *--outfile=$OUTFILE \       #[outfile=<output_params.txt>]
     *--fraczero=$FRACZERO \     #[fraczero=<fraction of zeroes{0.8}>]
     *--tres=$TRES \             #[tres=<threshold4bin[0,1]{0.6}>]
     *--split=$SPLIT \           #[split=<nr of partitions to split the stack into>]
     *--nthr=$NTHR               #[nthr=<nr of openMP threads{1}>]
     *
     */
    private void createInputEnvironmentVariables() {
        QName[] envVariables = { TITLE,
                                 STK,
                                 STK2,
                                 ORITAB,
                                 OUTSTK,
                                 DEFTAB,
                                 OUTFILE,
                                 SHALGN,
                                 MERGE,
                                 ROALGN,
                                 MASSCEN,
                                 VIS,
                                 BIN,
                                 ACF,
                                 PHRAND,
                                 HFUN,
                                 NORM,
                                 NOISE_NORM,
                                 AVG,
                                 RANKIFY,
                                 FILETAB,
                                 STATS,
                                 CTF,
                                 FT2IMG,
                                 COMPARE,
                                 NEG,
                                 CTFSQ,
                                 NPTCLS,
                                 FROMP,
                                 TRS,
                                 STATE,
                                 SYM_CLASS,
                                 TOP,
                                 NRAN,
                                 NEWBOX,
                                 FRAMEAVG,
                                 CLIP,
                                 BOX,
                                 SPLIT,
                                 NTHR,
                                 SMPD,
                                 HP,
                                 MUL,
                                 LP,
                                 FRAC,
                                 SNR,
                                 MSK,
                                 KV,
                                 SCALE,
                                 FRACZERO,
                                 FRACA ,
                                 CS,
                                 TRES,
                                 INNER,
                                 WIDTH };

        String[] initialValues = { "Title",
                                   "stack.ext",
                                   "stack2.ext",
                                   "SIMPLE alignment doc",
                                   "outstk.ext",
                                   "text file with defocus values",
                                   "output_params.txt",
                                   "no",
                                   "no",
                                   "no",
                                   "no",
                                   "no",
                                   "no",
                                   "no",
                                   "no",
                                   "sigm",
                                   "no",
                                   "no",
                                   "no",
                                   "no",
                                   "filenames.txt",
                                   "yes",
                                   "no",
                                   "no",
                                   "no",
                                   "no",
                                   "no", 
                                   "1",
                                   "0",
                                   "0",
                                   "1",
                                   "1",
                                   "0",
                                   "0",
                                   "0",
                                   "0",
                                   "128", 
                                   "256",
                                   "1",
                                   "1",
                                   "1.101",
                                   "100",
                                   "1",
                                   "20",
                                   "1",
                                   "1.0",
                                   "0",
                                   "300.",
                                   "1.",
                                   "0.8",
                                   "0.07", 
                                   "2.7",
                                   "0.6",
                                   "0.0",
                                   "10.0" };
		
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
    public URL getIconURL() {return getIconURL(SimpleStackopsGridBean.class, "images/Ribo_single_cropped.png");}

}
