/**
 * 
 */
package edu.kit.mmm.wrapper.simpleStackops;

import java.util.HashMap;

import org.apache.log4j.Logger;

/**
 * @author Frederic Bonnet
 * @Date 18th of November 2015
 *
 */
public class SimpleStackopsParameters {
    /** Logger */
    private static Logger log = Logger.getLogger(SimpleStackopsParameters.class.getName());
    
    private String stk;//=$STK \               #[stk=<stack.ext>]
    private String stk2;//=$STK2 \             #[stk2=<stack2.ext>]
    private String oritab;//=$ORITAB \         #[oritab=<SIMPLE alignment doc>]
    private String outstk;//=$OUTSTK \         #[outstk=<outstk.ext>]
    private String deftab;//=$DEFTAB \         #[deftab=<text file with defocus values>]
    private String outfile;//=$OUTFILE \       #[outfile=<output_params.txt>]
    private String shalgn;//=$SHALGN \         #[shalgn=<yes|no{no}>]
    private String merge;//=$MERGE \           #[merge=<yes|no{no}>]
    private String roalgn;//=$ROALGN  \        #[roalgn=<yes|no{no}>]
    private String masscen;//=$MASSCEN \       #[masscen=<yes|no{no}>]
    private String vis;//=$VIS \               #[vis=<yes|no>{no}]
    private String bin;//=$BIN \               #[bin=<binarize{no}>]
    private String acf;//=$ACF \               #[acf=<yes|no{no}>]
    private String phrand;//=$PHRAND \         #[phrand=<yes|no{no}>]
    private String hfun;//=$HFUN \             #[hfun=<sigm|tanh|lin{sigm}>]
    private String norm;//=$NORM \             #[norm=<yes|no{no}>]
    private String noise_norm;//=$NOISE_NORM \ #[noise_norm=<yes|no>{no}]
    private String avg;//=$AVG \               #[avg=<yes|no>{no}]
    private String rankify;//=$RANKIFY \       #[rankify=<yes|no>{no}]
    private String filetab;//=$FILETAB \       #[filetab=<filenames.txt>]
    private String stats;//=$STATS \           #[stats=<yes|no{yes}>]
    private String ctf;//=$CTF \               #[ctf=<yes|no|flip|mul|abs{no}>]
    private String ft2img;//=$FT2IMG \         #[ft2img=<yes|no{no}>]
    private String compare;//=$COMPARE \       #[compare=<yes|no{no}>]
    private String neg;//=$NEG \               #[neg=<yes|no{no}>]
    private String ctfsq;//=$CTFSQ \           #[ctfsq=<yes|no{no}>] 
    private int nptcls;//=$NPTCLS \            #[nptcls=<nr of imgs>{1}]
    private int fromp;//=$FROMP \              #[fromp=<start ptcl>{1}]
    private int trs;//=$TRS \                  #[trs=<origin]
    private int state;//=$STATE \              #[state=<state to extract>{1}]
    private int sym_class;//=$SYM_CLASS \      #[class=<symmetry class>{1}]
    private int top;//=$TOP \                  #[top=<stop ptcl>]
    private int nran;//=$NRAN \                #[nran=<number of random images to select>{0}]
    private int newbox;//=$NEWBOX \            #[newbox=<scaled box>{0}]
    private int frameavg;//=$FRAMEAVG \        #[frameavg=<nr of frames to average{0}>]
    private int clip;//=$CLIP \                #[clip=<clipped box size{box}>] 
    private int box;//=$BOX \                  #[box=<image size(in pixels)>]
    private int split;//=$SPLIT \              #[split=<nr of partitions to split the stack into>]
    private int nthr;//=$NTHR                  #[nthr=<nr of openMP threads{1}>]
    private float smpd;//=$SMPD \              #[smpd=<sampling distance(in A)>{1.101}]
    private float hp;//=$HP \                  #[hp=<high-pass limit(in A)>]
    private float mul;//=$MUL \                #[mul=<shift multiplication factor{1}>]
    private float lp;//=$LP \                  #[lp=<low-pass limitshift(in pixels){3}>]
    private float frac;//=$FRAC \              #[frac=<fraction of ptcls to extract{1}>]
    private float snr;//=$SNR \                #[snr=<signal2noise ratio>]
    private float msk;//=$MSK  \               #[msk=<mask radius(in pixels){box/2}>]
    private float kv;//=$KV \                  #[kv=<acceleration voltage(in kV){300.}>]
    private float scale;//=$SCALE \            #[scale=<scale factor{1.}>]
    private float fraczero;//=$FRACZERO \      #[fraczero=<fraction of zeroes{0.8}>]
    private float fraca;//=$FRACA  \           #[fraca=<frac amp contrast{0.07}>] 
    private float cs;//=$CS \                  #[cs=<spherical aberration constant(in mm){2.7}>]
    private float tres;//=$TRES \              #[tres=<threshold4bin[0,1]{0.6}>]
    private float inner;//=$INNER \            #[inner=<inner mask radius(in pixels)>]
    private float width;//=$WIDTH \            #[width=<pixels falloff inner mask{10}>]

    /**
     * @return the stk
     */
    public String getStk() {
        return stk;
    }
    /**
     * @param stk the stk to set
     */
    public void setStk(String stk) {
        this.stk = stk;
    }
    /**
     * @return the stk2
     */
    public String getStk2() {
        return stk2;
    }
    /**
     * @param stk2 the stk2 to set
     */
    public void setStk2(String stk2) {
        this.stk2 = stk2;
    }
    /**
     * @return the oritab
     */
    public String getOritab() {
        return oritab;
    }
    /**
     * @param oritab the oritab to set
     */
    public void setOritab(String oritab) {
        this.oritab = oritab;
    }
    /**
     * @return the outstk
     */
    public String getOutstk() {
        return outstk;
    }
    /**
     * @param outstk the outstk to set
     */
    public void setOutstk(String outstk) {
        this.outstk = outstk;
    }
    /**
     * @return the deftab
     */
    public String getDeftab() {
        return deftab;
    }
    /**
     * @param deftab the deftab to set
     */
    public void setDeftab(String deftab) {
        this.deftab = deftab;
    }
    /**
     * @return the outfile
     */
    public String getOutfile() {
        return outfile;
    }
    /**
     * @param outfile the outfile to set
     */
    public void setOutfile(String outfile) {
        this.outfile = outfile;
    }
    /**
     * @return the shalgn
     */
    public String getShalgn() {
        return shalgn;
    }
    /**
     * @param shalgn the shalgn to set
     */
    public void setShalgn(String shalgn) {
        this.shalgn = shalgn;
    }
    /**
     * @return the merge
     */
    public String getMerge() {
        return merge;
    }
    /**
     * @param merge the merge to set
     */
    public void setMerge(String merge) {
        this.merge = merge;
    }
    /**
     * @return the roalgn
     */
    public String getRoalgn() {
        return roalgn;
    }
    /**
     * @param roalgn the roalgn to set
     */
    public void setRoalgn(String roalgn) {
        this.roalgn = roalgn;
    }
    /**
     * @return the masscen
     */
    public String getMasscen() {
        return masscen;
    }
    /**
     * @param masscen the masscen to set
     */
    public void setMasscen(String masscen) {
        this.masscen = masscen;
    }
    /**
     * @return the vis
     */
    public String getVis() {
        return vis;
    }
    /**
     * @param vis the vis to set
     */
    public void setVis(String vis) {
        this.vis = vis;
    }
    /**
     * @return the bin
     */
    public String getBin() {
        return bin;
    }
    /**
     * @param bin the bin to set
     */
    public void setBin(String bin) {
        this.bin = bin;
    }
    /**
     * @return the acf
     */
    public String getAcf() {
        return acf;
    }
    /**
     * @param acf the acf to set
     */
    public void setAcf(String acf) {
        this.acf = acf;
    }
    /**
     * @return the phrand
     */
    public String getPhrand() {
        return phrand;
    }
    /**
     * @param phrand the phrand to set
     */
    public void setPhrand(String phrand) {
        this.phrand = phrand;
    }
    /**
     * @return the hfun
     */
    public String getHfun() {
        return hfun;
    }
    /**
     * @param hfun the hfun to set
     */
    public void setHfun(String hfun) {
        this.hfun = hfun;
    }
    /**
     * @return the norm
     */
    public String getNorm() {
        return norm;
    }
    /**
     * @param norm the norm to set
     */
    public void setNorm(String norm) {
        this.norm = norm;
    }
    /**
     * @return the noise_norm
     */
    public String getNoise_norm() {
        return noise_norm;
    }
    /**
     * @param noise_norm the noise_norm to set
     */
    public void setNoise_norm(String noise_norm) {
        this.noise_norm = noise_norm;
    }
    /**
     * @return the avg
     */
    public String getAvg() {
        return avg;
    }
    /**
     * @param avg the avg to set
     */
    public void setAvg(String avg) {
        this.avg = avg;
    }
    /**
     * @return the rankify
     */
    public String getRankify() {
        return rankify;
    }
    /**
     * @param rankify the rankify to set
     */
    public void setRankify(String rankify) {
        this.rankify = rankify;
    }
    /**
     * @return the filetab
     */
    public String getFiletab() {
        return filetab;
    }
    /**
     * @param filetab the filetab to set
     */
    public void setFiletab(String filetab) {
        this.filetab = filetab;
    }
    /**
     * @return the stats
     */
    public String getStats() {
        return stats;
    }
    /**
     * @param stats the stats to set
     */
    public void setStats(String stats) {
        this.stats = stats;
    }
    /**
     * @return the ctf
     */
    public String getCtf() {
        return ctf;
    }
    /**
     * @param ctf the ctf to set
     */
    public void setCtf(String ctf) {
        this.ctf = ctf;
    }
    /**
     * @return the ft2img
     */
    public String getFt2img() {
        return ft2img;
    }
    /**
     * @param ft2img the ft2img to set
     */
    public void setFt2img(String ft2img) {
        this.ft2img = ft2img;
    }
    /**
     * @return the compare
     */
    public String getCompare() {
        return compare;
    }
    /**
     * @param compare the compare to set
     */
    public void setCompare(String compare) {
        this.compare = compare;
    }
    /**
     * @return the neg
     */
    public String getNeg() {
        return neg;
    }
    /**
     * @param neg the neg to set
     */
    public void setNeg(String neg) {
        this.neg = neg;
    }
    /**
     * @return the ctfsq
     */
    public String getCtfsq() {
        return ctfsq;
    }
    /**
     * @param ctfsq the ctfsq to set
     */
    public void setCtfsq(String ctfsq) {
        this.ctfsq = ctfsq;
    }
    /**
     * @return the nptcls
     */
    public int getNptcls() {
        return nptcls;
    }
    /**
     * @param nptcls the nptcls to set
     */
    public void setNptcls(int nptcls) {
        this.nptcls = nptcls;
    }
    /**
     * @return the fromp
     */
    public int getFromp() {
        return fromp;
    }
    /**
     * @param fromp the fromp to set
     */
    public void setFromp(int fromp) {
        this.fromp = fromp;
    }
    /**
     * @return the trs
     */
    public int getTrs() {
        return trs;
    }
    /**
     * @param trs the trs to set
     */
    public void setTrs(int trs) {
        this.trs = trs;
    }
    /**
     * @return the state
     */
    public int getState() {
        return state;
    }
    /**
     * @param state the state to set
     */
    public void setState(int state) {
        this.state = state;
    }
    /**
     * @return the sym_class
     */
    public int getSym_class() {
        return sym_class;
    }
    /**
     * @param sym_class the sym_class to set
     */
    public void setSym_class(int sym_class) {
        this.sym_class = sym_class;
    }
    /**
     * @return the top
     */
    public int getTop() {
        return top;
    }
    /**
     * @param top the top to set
     */
    public void setTop(int top) {
        this.top = top;
    }
    /**
     * @return the nran
     */
    public int getNran() {
        return nran;
    }
    /**
     * @param nran the nran to set
     */
    public void setNran(int nran) {
        this.nran = nran;
    }
    /**
     * @return the newbox
     */
    public int getNewbox() {
        return newbox;
    }
    /**
     * @param newbox the newbox to set
     */
    public void setNewbox(int newbox) {
        this.newbox = newbox;
    }
    /**
     * @return the frameavg
     */
    public int getFrameavg() {
        return frameavg;
    }
    /**
     * @param frameavg the frameavg to set
     */
    public void setFrameavg(int frameavg) {
        this.frameavg = frameavg;
    }
    /**
     * @return the clip
     */
    public int getClip() {
        return clip;
    }
    /**
     * @param clip the clip to set
     */
    public void setClip(int clip) {
        this.clip = clip;
    }
    /**
     * @return the box
     */
    public int getBox() {
        return box;
    }
    /**
     * @param box the box to set
     */
    public void setBox(int box) {
        this.box = box;
    }
    /**
     * @return the split
     */
    public int getSplit() {
        return split;
    }
    /**
     * @param split the split to set
     */
    public void setSplit(int split) {
        this.split = split;
    }
    /**
     * @return the nthr
     */
    public int getNthr() {
        return nthr;
    }
    /**
     * @param nthr the nthr to set
     */
    public void setNthr(int nthr) {
        this.nthr = nthr;
    }
    /**
     * @return the smpd
     */
    public float getSmpd() {
        return smpd;
    }
    /**
     * @param smpd the smpd to set
     */
    public void setSmpd(float smpd) {
        this.smpd = smpd;
    }
    /**
     * @return the hp
     */
    public float getHp() {
        return hp;
    }
    /**
     * @param hp the hp to set
     */
    public void setHp(float hp) {
        this.hp = hp;
    }
    /**
     * @return the mul
     */
    public float getMul() {
        return mul;
    }
    /**
     * @param mul the mul to set
     */
    public void setMul(float mul) {
        this.mul = mul;
    }
    /**
     * @return the lp
     */
    public float getLp() {
        return lp;
    }
    /**
     * @param lp the lp to set
     */
    public void setLp(float lp) {
        this.lp = lp;
    }
    /**
     * @return the frac
     */
    public float getFrac() {
        return frac;
    }
    /**
     * @param frac the frac to set
     */
    public void setFrac(float frac) {
        this.frac = frac;
    }
    /**
     * @return the snr
     */
    public float getSnr() {
        return snr;
    }
    /**
     * @param snr the snr to set
     */
    public void setSnr(float snr) {
        this.snr = snr;
    }
    /**
     * @return the msk
     */
    public float getMsk() {
        return msk;
    }
    /**
     * @param msk the msk to set
     */
    public void setMsk(float msk) {
        this.msk = msk;
    }
    /**
     * @return the kv
     */
    public float getKv() {
        return kv;
    }
    /**
     * @param kv the kv to set
     */
    public void setKv(float kv) {
        this.kv = kv;
    }
    /**
     * @return the scale
     */
    public float getScale() {
        return scale;
    }
    /**
     * @param scale the scale to set
     */
    public void setScale(float scale) {
        this.scale = scale;
    }
    /**
     * @return the fraczero
     */
    public float getFraczero() {
        return fraczero;
    }
    /**
     * @param fraczero the fraczero to set
     */
    public void setFraczero(float fraczero) {
        this.fraczero = fraczero;
    }
    /**
     * @return the fraca
     */
    public float getFraca() {
        return fraca;
    }
    /**
     * @param fraca the fraca to set
     */
    public void setFraca(float fraca) {
        this.fraca = fraca;
    }
    /**
     * @return the cs
     */
    public float getCs() {
        return cs;
    }
    /**
     * @param cs the cs to set
     */
    public void setCs(float cs) {
        this.cs = cs;
    }
    /**
     * @return the tres
     */
    public float getTres() {
        return tres;
    }
    /**
     * @param tres the tres to set
     */
    public void setTres(float tres) {
        this.tres = tres;
    }
    /**
     * @return the inner
     */
    public float getInner() {
        return inner;
    }
    /**
     * @param inner the inner to set
     */
    public void setInner(float inner) {
        this.inner = inner;
    }
    /**
     * @return the width
     */
    public float getWidth() {
        return width;
    }
    /**
     * @param width the width to set
     */
    public void setWidth(float width) {
        this.width = width;
    }


    /**
     * Maps a {@link HashMap} which contains the Bigdft-Wrapper start parameters
     * to an instance of this class.
     * 
     * @param parameter
     *            the HashMap which contains the Bigdft-Wrapper start parameters
     */
    public void mapParameter(HashMap<String, String> parameter) {
        log.info("Map Wrapper parameter to SimpleStackopsParameters");

        //private String stk;//=$STK \               #[stk=<stack.ext>]
        if ( parameter.containsKey("--stk")) {
            stk = parameter.get("--stk");
        }
        //private String stk2;//=$STK2 \             #[stk2=<stack2.ext>]
        if ( parameter.containsKey("--stk2")) {
            stk2 = parameter.get("--stk2");
        }
        //private String oritab;//=$ORITAB \         #[oritab=<SIMPLE alignment doc>]
        if ( parameter.containsKey("--oritab")) {
            oritab = parameter.get("--oritab");
        }
        //private String outstk;//=$OUTSTK \         #[outstk=<outstk.ext>]
        if ( parameter.containsKey("--outstk")) {
            outstk = parameter.get("--outstk");
        }
        //private String deftab;//=$DEFTAB \         #[deftab=<text file with defocus values>]
        if ( parameter.containsKey("--stk")) {
            stk = parameter.get("--stk");
        }
        //private String outfile;//=$OUTFILE \       #[outfile=<output_params.txt>]
        if ( parameter.containsKey("--stk")) {
            stk = parameter.get("--stk");
        }
        //private String shalgn;//=$SHALGN \         #[shalgn=<yes|no{no}>]
        if ( parameter.containsKey("--shalgn")) {
            shalgn = parameter.get("--shalgn");
        }
        //private String merge;//=$MERGE \           #[merge=<yes|no{no}>]
        if ( parameter.containsKey("--merge")) {
            merge = parameter.get("--merge");
        }
        //private String roalgn;//=$ROALGN  \        #[roalgn=<yes|no{no}>]
        if ( parameter.containsKey("--roalgn")) {
            roalgn = parameter.get("--roalgn");
        }
        //private String masscen;//=$MASSCEN \       #[masscen=<yes|no{no}>]
        if ( parameter.containsKey("--masscen")) {
            masscen = parameter.get("--smasscen");
        }
        //private String vis;//=$VIS \               #[vis=<yes|no>{no}]
        if ( parameter.containsKey("--vis")) {
            vis = parameter.get("--vis");
        }
        //private String bin;//=$BIN \               #[bin=<binarize{no}>]
        if ( parameter.containsKey("--bin")) {
            bin = parameter.get("--bin");
        }
        //private String acf;//=$ACF \               #[acf=<yes|no{no}>]
        if ( parameter.containsKey("--acf")) {
            acf = parameter.get("--acf");
        }
        //private String phrand;//=$PHRAND \         #[phrand=<yes|no{no}>]
        if ( parameter.containsKey("--phrand")) {
            phrand = parameter.get("--phrand");
        }
        //private String hfun;//=$HFUN \             #[hfun=<sigm|tanh|lin{sigm}>]
        if ( parameter.containsKey("--hfun")) {
            hfun = parameter.get("--hfun");
        }
        //private String norm;//=$NORM \             #[norm=<yes|no{no}>]
        if ( parameter.containsKey("--norm")) {
            norm = parameter.get("--norm");
        }
        //private String noise_norm;//=$NOISE_NORM \ #[noise_norm=<yes|no>{no}]
        if ( parameter.containsKey("--noise_norm")) {
            noise_norm = parameter.get("--noise_norm");
        }
        //private String avg;//=$AVG \               #[avg=<yes|no>{no}]
        if ( parameter.containsKey("--avg")) {
            avg = parameter.get("--avg");
        }
        //private String rankify;//=$RANKIFY \       #[rankify=<yes|no>{no}]
        if ( parameter.containsKey("--rankify")) {
            rankify = parameter.get("--rankify");
        }
        //private String filetab;//=$FILETAB \       #[filetab=<filenames.txt>]
        if ( parameter.containsKey("--filetab")) {
            filetab = parameter.get("--filetab");
        }
        //private String stats;//=$STATS \           #[stats=<yes|no{yes}>]
        if ( parameter.containsKey("--stats")) {
            stats = parameter.get("--stats");
        }
        //private String ctf;//=$CTF \               #[ctf=<yes|no|flip|mul|abs{no}>]
        if ( parameter.containsKey("--ctf")) {
            ctf = parameter.get("--ctf");
        }
        //private String ft2img;//=$FT2IMG \         #[ft2img=<yes|no{no}>]
        if ( parameter.containsKey("--ft2img")) {
            ft2img = parameter.get("--ft2img");
        }
        //private String compare;//=$COMPARE \       #[compare=<yes|no{no}>]
        if ( parameter.containsKey("--compare")) {
            compare = parameter.get("--compare");
        }
        //private String neg;//=$NEG \               #[neg=<yes|no{no}>]
        if ( parameter.containsKey("--neg")) {
            neg = parameter.get("--neg");
        }
        //private String ctfsq;//=$CTFSQ \           #[ctfsq=<yes|no{no}>] 
        if ( parameter.containsKey("--ctfsq")) {
            ctfsq = parameter.get("--ctfsq");
        }
        //private int nptcls;//=$NPTCLS \            #[nptcls=<nr of imgs>{1}]
        if ( parameter.containsKey("--nptcls")) {
            nptcls = Integer.parseInt(parameter.get("--nptcls"));
        }
        //private int fromp;//=$FROMP \              #[fromp=<start ptcl>{1}]
        if ( parameter.containsKey("--fromp")) {
            fromp = Integer.parseInt(parameter.get("--fromp"));
        }
        //private int trs;//=$TRS \                  #[trs=<origin]
        if ( parameter.containsKey("--trs")) {
            trs = Integer.parseInt(parameter.get("--trs"));
        }
        //private int state;//=$STATE \              #[state=<state to extract>{1}]
        if ( parameter.containsKey("--state")) {
            state = Integer.parseInt(parameter.get("--state"));
        }
        //private int sym_class;//=$SYM_CLASS \      #[class=<symmetry class>{1}]
        if ( parameter.containsKey("--sym_class")) {
            sym_class = Integer.parseInt(parameter.get("--sym_class"));
        }
        //private int top;//=$TOP \                  #[top=<stop ptcl>]
        if ( parameter.containsKey("--top")) {
            top = Integer.parseInt(parameter.get("--top"));
        }
        //private int nran;//=$NRAN \                #[nran=<number of random images to select>{0}]
        if ( parameter.containsKey("--nran")) {
            nran = Integer.parseInt(parameter.get("--nran"));
        }
        //private int newbox;//=$NEWBOX \            #[newbox=<scaled box>{0}]
        if ( parameter.containsKey("--newbox")) {
            newbox = Integer.parseInt(parameter.get("--newbox"));
        }
        //private int frameavg;//=$FRAMEAVG \        #[frameavg=<nr of frames to average{0}>]
        if ( parameter.containsKey("--frameavg")) {
           frameavg = Integer.parseInt(parameter.get("--frameavg"));
        }
        //private int clip;//=$CLIP \                #[clip=<clipped box size{box}>] 
        if ( parameter.containsKey("--clip")) {
            clip = Integer.parseInt(parameter.get("--clip"));
        }
        //private int box;//=$BOX \                  #[box=<image size(in pixels)>]
        if ( parameter.containsKey("--box")) {
            box = Integer.parseInt(parameter.get("--box"));
        }
        //private int split;//=$SPLIT \              #[split=<nr of partitions to split the stack into>]
        if ( parameter.containsKey("--split")) {
            split = Integer.parseInt(parameter.get("--split"));
        }
        //private int nthr;//=$NTHR                  #[nthr=<nr of openMP threads{1}>]
        if ( parameter.containsKey("--nthr")) {
            nthr = Integer.parseInt(parameter.get("--nthr"));
        }
        //private float smpd;//=$SMPD \              #[smpd=<sampling distance(in A)>{1.101}]
        if (parameter.containsKey("--smpd")) {
            smpd = Float.parseFloat(parameter.get("--smpd"));
        }
        //private float hp;//=$HP \                  #[hp=<high-pass limit(in A)>]
        if (parameter.containsKey("--hp")) {
            hp = Float.parseFloat(parameter.get("--hp"));
        }
        //private float mul;//=$MUL \                #[mul=<shift multiplication factor{1}>]
        if (parameter.containsKey("--mul")) {
            mul = Float.parseFloat(parameter.get("--mul"));
        }
        //private float lp;//=$LP \                  #[lp=<low-pass limitshift(in pixels){3}>]
        if (parameter.containsKey("--lp")) {
            lp = Float.parseFloat(parameter.get("--lp"));
        }
        //private float frac;//=$FRAC \              #[frac=<fraction of ptcls to extract{1}>]
        if (parameter.containsKey("--frac")) {
            frac = Float.parseFloat(parameter.get("--frac"));
        }
        //private float snr;//=$SNR \                #[snr=<signal2noise ratio>]
        if (parameter.containsKey("--snr")) {
            snr = Float.parseFloat(parameter.get("--snr"));
        }
        //private float msk;//=$MSK  \               #[msk=<mask radius(in pixels){box/2}>]
        if (parameter.containsKey("--msk")) {
            msk = Float.parseFloat(parameter.get("--msk"));
        }
        //private float kv;//=$KV \                  #[kv=<acceleration voltage(in kV){300.}>]
        if (parameter.containsKey("--kv")) {
            kv = Float.parseFloat(parameter.get("--kv"));
        }
        //private float scale;//=$SCALE \            #[scale=<scale factor{1.}>]
        if (parameter.containsKey("--scale")) {
            scale = Float.parseFloat(parameter.get("--scale"));
        }
        //private float fraczero;//=$FRACZERO \      #[fraczero=<fraction of zeroes{0.8}>]
        if (parameter.containsKey("--fraczero")) {
            fraczero = Float.parseFloat(parameter.get("--fraczero"));
        }
        //private float fraca;//=$FRACA  \           #[fraca=<frac amp contrast{0.07}>] 
        if (parameter.containsKey("--fraca")) {
            fraca = Float.parseFloat(parameter.get("--fraca"));
        }
        //private float cs;//=$CS \                  #[cs=<spherical aberration constant(in mm){2.7}>]
        if (parameter.containsKey("--cs")) {
            cs = Float.parseFloat(parameter.get("--cs"));
        }
        //private float tres;//=$TRES \              #[tres=<threshold4bin[0,1]{0.6}>]
        if (parameter.containsKey("--tres")) {
            tres = Float.parseFloat(parameter.get("--tres"));
        }
        //private float inner;//=$INNER \            #[inner=<inner mask radius(in pixels)>]
        if (parameter.containsKey("--inner")) {
            inner = Float.parseFloat(parameter.get("--inner"));
        }
        //private float width;//=$WIDTH \            #[width=<pixels falloff inner mask{10}>]
        if (parameter.containsKey("--width")) {
            width = Float.parseFloat(parameter.get("--width"));
        }

    }
}
