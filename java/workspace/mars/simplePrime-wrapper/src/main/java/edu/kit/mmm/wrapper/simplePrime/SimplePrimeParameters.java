/**
 * 
 */
package edu.kit.mmm.wrapper.simplePrime;

import java.util.HashMap;

import org.apache.log4j.Logger;

/**
 * @author Frederic Bonnet
 * @Date 18th of November 2015
 *
 */
public class SimplePrimeParameters {
    /** Logger */
    private static Logger log = Logger.getLogger(SimplePrimeParameters.class.getName());
	
    private String stk;//=$STK \                      #stk=<stack.ext{sumstack.mrc}>
    private String vol1;//=$VOL1 \                    #[vol1=<invol.ext>]
    private String vol2;//=$VOL2 \                    #[vol2=<refvol_2.ext> etc.]
    private String oritab;//=$ORITAB \                #[oritab=<previous alignment doc>]
    private String dynlp;//=$DYNLP \                  #[dynlp=<yes|no{yes}>]
    private String refine;//=$REFINE \                #[refine=<no|exhaust|shift|neigh{no}>]
    private String deftab;//=$DEFTAB \                #[deftab=<text file defocus values>]
    private String eo;//=$EO \                        #[eo=<yes|no{no}>]
    private String pgrp;//=$PGRP \                    #[pgrp=<cn|dn|t|o|i{c1}>]
    private String ctf;//=$CTF \                      #[ctf=<yes|no|flip|mul{no}>]
    private String diversify;//=$DIVERSIFY \          #[diversify=<yes|no>]
    private String noise;//=$NOISE \                  #[noise=<yes|no{no}>]
    private String outfile;//=$OUTFILE \              #[outfile=<output alignment doc 4 parallell jobs>]
    private String norec;//=$NOREC                    #[norec=<yes|no{no}>]
    private int	nstates;//=$NSTATES \                 #[nstates=<nstates to reconstruct>]
    private int	nthr;//=$NTHR \                       #[nthr=<nr of OpenMP threads{1}>]
    private int	startit;//=$STARTIT \                 #[startit=<start iteration>]
    private int	nspace;//=$NSPACE \                   #[nspace=<nr reference sections{1000}>]
    private int	nnn;//=$NNN \                         #[nnn=<nr nearest neighbors{300}>]
    private int	maxits;//=$MAXITS \                   #[maxits=<max iterations{100}>]
    private int	edge;//=$EDGE \                       #[edge=<edge size for softening molecular envelope(in pixels){3}>] 
    private int	find;//=$FIND \                       #[find=<Fourier index>]
    private int	fstep;//=$FSTEP \                     #[fstep=<Fourier step size>]
    private int	trsstep;//=$TRSSTEP \                 #[trsstep=<origin shift stepsize{1}>] 
    private int	nvox;//=$NVOX \                       #[nvox=<nr of voxels in mask{0}>] 
    private float smpd;//=$SMPD \                     #smpd=<sampling distance(in A){1.101}>
    private float msk;//=$MSK \                       #msk=<mask radius(in pixels)>
    private float trs;//=$TRS \                       #[trs=<origin shift(in pixels){0}>]
    private float lp;//=$LP \                         #[lp=<low-pass limit{20}>]
    private float frac;//=$FRAC \                     #[frac=<fraction of ptcls to include{1}>]
    private float mw;//=$MW \                         #[mw=<molecular weight(in kD)>]
    private float lpstop;//=$LPSTOP \                 #[lpstop=<stay at this low-pass limit (in A)>]
    private float amsklp;//=$AMSKLP \                 #[amsklp=<automask low-pass limit(in A)>] 
    private float kv;//=$KV \                         #[kv=<acceleration voltage(in kV){300.}>]
    private float cs;//=$CS \                         #[cs=<spherical aberration constant(in mm){2.7}>]
    private float fraca;//=$FRACA \                   #[fraca=<frac amp contrast{0.07}>]
    private float hp;//=$HP \                         #[hp=<high-pass limit(in A)>]
    private float fracvalid;//=$FRACVALID \           #[fracvalid=<fraction of particles 4 validation{0.}>]
    private float lpvalid;//=$LPVALID \               #[lpvalid=<low-pass limit validptcls{20}>]
    private float time_per_image;//=$TIME_PER_IMAGE \ #[time_per_image=<{100}>] 
    private float dens;//=$DENS \                     #[dens=<density(e.g. 9.368 Da/A3 4 gold clusters){0.}>]
    private float inner;//=$INNER \                   #[inner=<inner mask radius(in pixels)>]
    private float width;//=$WIDTH \                   #[width=<pixels falloff inner mask{10}>]
	
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
     * @return the vol1
     */
    public String getVol1() {
        return vol1;
    }

    /**
     * @param vol1 the vol1 to set
     */
    public void setVol1(String vol1) {
        this.vol1 = vol1;
    }

    /**
     * @return the vol2
     */
    public String getVol2() {
        return vol2;
    }

    /**
     * @param vol2 the vol2 to set
     */
    public void setVol2(String vol2) {
        this.vol2 = vol2;
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
     * @return the dynlp
     */
    public String getDynlp() {
        return dynlp;
    }

    /**
     * @param dynlp the dynlp to set
     */
    public void setDynlp(String dynlp) {
        this.dynlp = dynlp;
    }

    /**
     * @return the refine
     */
    public String getRefine() {
        return refine;
    }

    /**
     * @param refine the refine to set
     */
    public void setRefine(String refine) {
        this.refine = refine;
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
     * @return the eo
     */
    public String getEo() {
        return eo;
    }

    /**
     * @param eo the eo to set
     */
    public void setEo(String eo) {
        this.eo = eo;
    }

    /**
     * @return the pgrp
     */
    public String getPgrp() {
        return pgrp;
    }

    /**
     * @param pgrp the pgrp to set
     */
    public void setPgrp(String pgrp) {
        this.pgrp = pgrp;
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
     * @return the diversify
     */
    public String getDiversify() {
        return diversify;
    }

    /**
     * @param diversify the diversify to set
     */
    public void setDiversify(String diversify) {
        this.diversify = diversify;
    }

    /**
     * @return the noise
     */
    public String getNoise() {
        return noise;
    }

    /**
     * @param noise the noise to set
     */
    public void setNoise(String noise) {
        this.noise = noise;
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
     * @return the norec
     */
    public String getNorec() {
        return norec;
    }

    /**
     * @param norec the norec to set
     */
    public void setNorec(String norec) {
        this.norec = norec;
    }

    /**
     * @return the nstates
     */
    public int getNstates() {
        return nstates;
    }

    /**
     * @param nstates the nstates to set
     */
    public void setNstates(int nstates) {
        this.nstates = nstates;
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
     * @return the startit
     */
    public int getStartit() {
        return startit;
    }

    /**
     * @param startit the startit to set
     */
    public void setStartit(int startit) {
        this.startit = startit;
    }

    /**
     * @return the nspace
     */
    public int getNspace() {
        return nspace;
    }

    /**
     * @param nspace the nspace to set
     */
    public void setNspace(int nspace) {
        this.nspace = nspace;
    }

    /**
     * @return the nnn
     */
    public int getNnn() {
        return nnn;
    }

    /**
     * @param nnn the nnn to set
     */
    public void setNnn(int nnn) {
        this.nnn = nnn;
    }

    /**
     * @return the maxits
     */
    public int getMaxits() {
        return maxits;
    }

    /**
     * @param maxits the maxits to set
     */
    public void setMaxits(int maxits) {
        this.maxits = maxits;
    }

    /**
     * @return the edge
     */
    public int getEdge() {
        return edge;
    }

    /**
     * @param edge the edge to set
     */
    public void setEdge(int edge) {
        this.edge = edge;
    }

    /**
     * @return the find
     */
    public int getFind() {
        return find;
    }

    /**
     * @param find the find to set
     */
    public void setFind(int find) {
        this.find = find;
    }

    /**
     * @return the fstep
     */
    public int getFstep() {
        return fstep;
    }

    /**
     * @param fstep the fstep to set
     */
    public void setFstep(int fstep) {
        this.fstep = fstep;
    }

    /**
     * @return the trsstep
     */
    public int getTrsstep() {
        return trsstep;
    }

    /**
     * @param trsstep the trsstep to set
     */
    public void setTrsstep(int trsstep) {
        this.trsstep = trsstep;
    }

    /**
     * @return the nvox
     */
    public int getNvox() {
        return nvox;
    }

    /**
     * @param nvox the nvox to set
     */
    public void setNvox(int nvox) {
        this.nvox = nvox;
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
     * @return the trs
     */
    public float getTrs() {
        return trs;
    }

    /**
     * @param trs the trs to set
     */
    public void setTrs(float trs) {
        this.trs = trs;
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
     * @return the mw
     */
    public float getMw() {
        return mw;
    }

    /**
     * @param mw the mw to set
     */
    public void setMw(float mw) {
        this.mw = mw;
    }

    /**
     * @return the lpstop
     */
    public float getLpstop() {
        return lpstop;
    }

    /**
     * @param lpstop the lpstop to set
     */
    public void setLpstop(float lpstop) {
        this.lpstop = lpstop;
    }

    /**
     * @return the amsklp
     */
    public float getAmsklp() {
        return amsklp;
    }

    /**
     * @param amsklp the amsklp to set
     */
    public void setAmsklp(float amsklp) {
        this.amsklp = amsklp;
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
     * @return the fracvalid
     */
    public float getFracvalid() {
        return fracvalid;
    }

    /**
     * @param fracvalid the fracvalid to set
     */
    public void setFracvalid(float fracvalid) {
        this.fracvalid = fracvalid;
    }

    /**
     * @return the lpvalid
     */
    public float getLpvalid() {
        return lpvalid;
    }

    /**
     * @param lpvalid the lpvalid to set
     */
    public void setLpvalid(float lpvalid) {
        this.lpvalid = lpvalid;
    }

    /**
     * @return the time_per_image
     */
    public float getTime_per_image() {
        return time_per_image;
    }

    /**
     * @param time_per_image the time_per_image to set
     */
    public void setTime_per_image(float time_per_image) {
        this.time_per_image = time_per_image;
    }

    /**
     * @return the dens
     */
    public float getDens() {
        return dens;
    }

    /**
     * @param dens the dens to set
     */
    public void setDens(float dens) {
        this.dens = dens;
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
        log.info("Map Wrapper parameter to SimplePrimeParameters");
        
        //private String stk;//=$STK \                      #stk=<stack.ext{sumstack.mrc}>
        if ( parameter.containsKey("--stk")) {
            stk = parameter.get("--stk");
        }
        //private String vol1;//=$VOL1 \                    #[vol1=<invol.ext>]
        if ( parameter.containsKey("--vol1")) {
            vol1 = parameter.get("--vol1");
        }
        //private String vol2;//=$VOL2 \                    #[vol2=<refvol_2.ext> etc.]
        if ( parameter.containsKey("--vol2")) {
            vol2 = parameter.get("--vol2");
        }
        //private String oritab;//=$ORITAB \                #[oritab=<previous alignment doc>]
        if ( parameter.containsKey("--oritab")) {
            oritab = parameter.get("--oritab");
        }
        //private String dynlp;//=$DYNLP \                  #[dynlp=<yes|no{yes}>]
        if ( parameter.containsKey("--dynlp")) {
            dynlp = parameter.get("--dynlp");
        }
        //private String refine;//=$REFINE \                #[refine=<no|exhaust|shift|neigh{no}>]
        if ( parameter.containsKey("--refine")) {
            refine = parameter.get("--refine");
        }
        //private String deftab;//=$DEFTAB \                #[deftab=<text file defocus values>]
        if ( parameter.containsKey("--deftab")) {
            deftab = parameter.get("--deftab");
        }
        //private String eo;//=$EO \                        #[eo=<yes|no{no}>]
        if ( parameter.containsKey("--eo")) {
            eo = parameter.get("--eo");
        }
        //private String pgrp;//=$PGRP \                    #[pgrp=<cn|dn|t|o|i{c1}>]
        if ( parameter.containsKey("--pgrp")) {
            pgrp = parameter.get("--pgrp");
        }
        //private String ctf;//=$CTF \                      #[ctf=<yes|no|flip|mul{no}>]
        if ( parameter.containsKey("--stk")) {
            stk = parameter.get("--stk");
        }
        //private String diversify;//=$DIVERSIFY \          #[diversify=<yes|no>]
        if ( parameter.containsKey("--diversify")) {
            diversify = parameter.get("--diversify");
        }
        //private String noise;//=$NOISE \                  #[noise=<yes|no{no}>]
        if ( parameter.containsKey("--noise")) {
            noise = parameter.get("--noise");
        }
        //private String outfile;//=$OUTFILE \              #[outfile=<output alignment doc 4 parallell jobs>]
        if ( parameter.containsKey("--outfile")) {
            outfile = parameter.get("--outfile");
        }
        //private String norec;//=$NOREC                    #[norec=<yes|no{no}>]
        if ( parameter.containsKey("--norec")) {
            norec = parameter.get("--norec");
        }
        //private int nstates;//=$NSTATES \                 #[nstates=<nstates to reconstruct>]
        if ( parameter.containsKey("--nstates")) {
            nstates = Integer.parseInt(parameter.get("--nstates"));
        }
        //private int nthr;//=$NTHR \                       #[nthr=<nr of OpenMP threads{1}>]
        if ( parameter.containsKey("--nthr")) {
            nthr = Integer.parseInt(parameter.get("--nthr"));
        }
        //private int startit;//=$STARTIT \                 #[startit=<start iteration>]
        if ( parameter.containsKey("--startit")) {
            startit = Integer.parseInt(parameter.get("--startit"));
        }
        //private int nspace;//=$NSPACE \                   #[nspace=<nr reference sections{1000}>]
        if ( parameter.containsKey("--nspace")) {
            nspace = Integer.parseInt(parameter.get("--nspace"));
        }
        //private int nnn;//=$NNN \                         #[nnn=<nr nearest neighbors{300}>]
        if ( parameter.containsKey("--nnnn")) {
            nnn = Integer.parseInt(parameter.get("--nnn"));
        }
        //private int maxits;//=$MAXITS \                   #[maxits=<max iterations{100}>]
        if ( parameter.containsKey("--maxits")) {
            maxits = Integer.parseInt(parameter.get("--maxits"));
        }
        //private int edge;//=$EDGE \                       #[edge=<edge size for softening molecular envelope(in pixels){3}>] 
        if ( parameter.containsKey("--edge")) {
            edge = Integer.parseInt(parameter.get("--edge"));
        }
        //private int find;//=$FIND \                       #[find=<Fourier index>]
        if ( parameter.containsKey("--find")) {
            find = Integer.parseInt(parameter.get("--find"));
        }
        //private int fstep;//=$FSTEP \                     #[fstep=<Fourier step size>]
        if ( parameter.containsKey("--fstep")) {
            fstep = Integer.parseInt(parameter.get("--fstep"));
        }
        //private int trsstep;//=$TRSSTEP \                 #[trsstep=<origin shift stepsize{1}>] 
        if ( parameter.containsKey("--trsstep")) {
           trsstep = Integer.parseInt(parameter.get("--trsstep"));
        }
        //private int nvox;//=$NVOX \                       #[nvox=<nr of voxels in mask{0}>] 
        if ( parameter.containsKey("--nvox")) {
            nvox = Integer.parseInt(parameter.get("--nvox"));
        }
        //private float smpd;//=$SMPD \                     #smpd=<sampling distance(in A){1.101}>
        if (parameter.containsKey("--smpd")) {
            smpd = Float.parseFloat(parameter.get("--smpd"));
        }
        //private float msk;//=$MSK \                       #msk=<mask radius(in pixels)>
        if (parameter.containsKey("--msk")) {
            msk = Float.parseFloat(parameter.get("--msk"));
        }
        //private float trs;//=$TRS \                       #[trs=<origin shift(in pixels){0}>]
        if (parameter.containsKey("--trs")) {
            trs = Float.parseFloat(parameter.get("--trs"));
        }
        //private float lp;//=$LP \                         #[lp=<low-pass limit{20}>]
        if (parameter.containsKey("--lp")) {
            lp = Float.parseFloat(parameter.get("--lp"));
        }
        //private float frac;//=$FRAC \                     #[frac=<fraction of ptcls to include{1}>]
        if (parameter.containsKey("--frac")) {
            frac = Float.parseFloat(parameter.get("--frac"));
        }
        //private float mw;//=$MW \                         #[mw=<molecular weight(in kD)>]
        if (parameter.containsKey("--mw")) {
            mw = Float.parseFloat(parameter.get("--mw"));
        }
        //private float lpstop;//=$LPSTOP \                 #[lpstop=<stay at this low-pass limit (in A)>]
        if (parameter.containsKey("--lpstop")) {
            lpstop = Float.parseFloat(parameter.get("--lpstop"));
        }
        //private float amsklp;//=$AMSKLP \                 #[amsklp=<automask low-pass limit(in A)>] 
        if (parameter.containsKey("--amsklp")) {
            amsklp = Float.parseFloat(parameter.get("--amsklp"));
        }
        //private float kv;//=$KV \                         #[kv=<acceleration voltage(in kV){300.}>]
        if (parameter.containsKey("--kv")) {
            kv = Float.parseFloat(parameter.get("--kv"));
        }
        //private float cs;//=$CS \                         #[cs=<spherical aberration constant(in mm){2.7}>]
        if (parameter.containsKey("--cs")) {
            cs = Float.parseFloat(parameter.get("--cs"));
        }
        //private float fraca;//=$FRACA \                   #[fraca=<frac amp contrast{0.07}>]
        if (parameter.containsKey("--fraca")) {
            fraca = Float.parseFloat(parameter.get("--fraca"));
        }
        //private float hp;//=$HP \                         #[hp=<high-pass limit(in A)>]
        if (parameter.containsKey("--hp")) {
            hp = Float.parseFloat(parameter.get("--hp"));
        }
        //private float fracvalid;//=$FRACVALID \           #[fracvalid=<fraction of particles 4 validation{0.}>]
        if (parameter.containsKey("--fracvalid")) {
            fracvalid = Float.parseFloat(parameter.get("--fracvalid"));
        }
        //private float lpvalid;//=$LPVALID \               #[lpvalid=<low-pass limit validptcls{20}>]
        if (parameter.containsKey("--lpvalid")) {
            lpvalid = Float.parseFloat(parameter.get("--lpvalid"));
        }
        //private float time_per_image;//=$TIME_PER_IMAGE \ #[time_per_image=<{100}>] 
        if (parameter.containsKey("--time_per_image")) {
            time_per_image = Float.parseFloat(parameter.get("--time_per_image"));
        }
        //private float dens;//=$DENS \                     #[dens=<density(e.g. 9.368 Da/A3 4 gold clusters){0.}>]
        if (parameter.containsKey("--dens")) {
            dens = Float.parseFloat(parameter.get("--dens"));
        }
        //private float inner;//=$INNER \                   #[inner=<inner mask radius(in pixels)>]
        if (parameter.containsKey("--inner")) {
            inner = Float.parseFloat(parameter.get("--inner"));
        }
        //private float width;//=$WIDTH \                   #[width=<pixels falloff inner mask{10}>]
        if (parameter.containsKey("--width")) {
            width = Float.parseFloat(parameter.get("--width"));
        }

    }

}
