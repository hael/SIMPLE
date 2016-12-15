/**
 * 
 */
package edu.kit.mmm.wrapper.simpleEoRecVol;

import java.util.HashMap;
import java.util.logging.Logger;

/**
 * @author Frederic Bonnet
 * @Date 18th of November 2015
 *
 */
public class SimpleEoRecVolParameters {
/** Logger */
	private static Logger log = Logger.getLogger(SimpleEoRecVolParameters.class.getName());

    private String stk;//=$STK \       # stk=<ptcls.ext>
    private float msk; //=$MSK \       # msk=<mask radius(in pixels)>
    private float smpd;//=$SMPD \      # smpd=<sampling distance(in A)>
    private String oritab;//=$ORITAB \ # oritab=<algndoc.txt>
    private float frac;//=$FRAC \      # [frac=<fraction of ptcls to include{1.}>]
    private int nthr;//=$NTHR \        # [nthr=<nr of openMP threads{1}>]
    private String pgrp;//=$PGRP \     # [pgrp=<cn|dn|t|o|i{c1}>]
    private float mul;//=$MUL \        # [mul=<shift multiplication factor{1}>]
    private int part;//=$PART \        # [part=<partition number>]
    private int	fromp;//=$FROMP \      # [fromp=<start ptcl{1)>]
    private int	top;//=$TOP \          # [top=<stop ptcl{nptcls}>]
    private String ctf;//=$CTF \       # [ctf=<yes|no|flip|mul{no}>]
    private float kv;//=$KV \          # [kv=<acceleration voltage(in kV){300.}>]
    private float fraca;//=$FRACA  \   # [fraca=<frac amp contrast{0.07}>]
    private float cs;//=$CS \          # [cs=<spherical aberration constant(in mm){2.7}>]
    private String deftab;//=$DEFTAB \ # [deftab=<text file defocus values>]
    private int	state;//=$STATE \      # [state=<state to reconstruct{all}>]
    private float alpha;//=$ALPHA \    # [alpha=<oversampling ratio or padding factor{2.}>]
    private float mw;//=$MW \          # [mw=<molecular weight(in kD)>]
    private float amsklp;//=$AMSKLP \  # [amsklp=<low-pass limit(in A){20}>]
    private int	edge;//=$EDGE  \       # [edge=<edge size for softening molecular envelope(in pixels){3}>]
    private float dens;//=$DENS \      # [dens=<density(e.g.9.368 Da/A3 4 gold clusters){0.}>]
    private int	box;//=$BOX \          # [box=<image size(in pixels)>]
    private float inner;//=$INNER \    # [inner=<inner mask radius(in pixels)>]
    private float width;//=$WIDTH      # [width=<pixels falloff inner mask{10}>]
	
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
     * @return the part
     */
    public int getPart() {
        return part;
    }

    /**
     * @param part the part to set
     */
    public void setPart(int part) {
        this.part = part;
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
     * @return the alpha
     */
    public float getAlpha() {
        return alpha;
    }

    /**
     * @param alpha the alpha to set
     */
    public void setAlpha(float alpha) {
        this.alpha = alpha;
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
        log.info("Map Wrapper parameter to SimpleEoRecVolParameters");

        // private String stk;//=$STK \       # stk=<ptcls.ext>
        if ( parameter.containsKey("--stk")) {
            stk = parameter.get("--stk");
        }
        // private float msk; //=$MSK \       # msk=<mask radius(in pixels)>
        if (parameter.containsKey("--msk")) {
            msk = Float.parseFloat(parameter.get("--msk"));
        }
        // private float smpd;//=$SMPD \      # smpd=<sampling distance(in A)>
        if (parameter.containsKey("--smpd")) {
            smpd = Float.parseFloat(parameter.get("--smpd"));
        }
        // private String oritab;//=$ORITAB \ # oritab=<algndoc.txt>
        if ( parameter.containsKey("--oritab")) {
            oritab = parameter.get("--oritab");
        }
        // private float frac;//=$FRAC \      # [frac=<fraction of ptcls to include{1.}>]
        if (parameter.containsKey("--frac")) {
            frac = Float.parseFloat(parameter.get("--frac"));
        }
        // private int nthr;//=$NTHR \        # [nthr=<nr of openMP threads{1}>]
        if ( parameter.containsKey("--nthr")) {
            nthr = Integer.parseInt(parameter.get("--nthr"));
        }
        // private String pgrp;//=$PGRP \     # [pgrp=<cn|dn|t|o|i{c1}>]
        if ( parameter.containsKey("--pgrp")) {
            pgrp = parameter.get("--pgrp");
        }
        // private float mul;//=$MUL \        # [mul=<shift multiplication factor{1}>]
        if (parameter.containsKey("--mul")) {
            mul = Float.parseFloat(parameter.get("--mul"));
        }
        // private int part;//=$PART \        # [part=<partition number>]
        if ( parameter.containsKey("--part")) {
            part = Integer.parseInt(parameter.get("--part"));
        }
        // private int fromp;//=$FROMP \      # [fromp=<start ptcl{1)>]
        if ( parameter.containsKey("--fromp")) {
            fromp = Integer.parseInt(parameter.get("--fromp"));
        }
        // private int top;//=$TOP \          # [top=<stop ptcl{nptcls}>]
        if ( parameter.containsKey("--top")) {
            top = Integer.parseInt(parameter.get("--top"));
        }
        // private String ctf;//=$CTF \       # [ctf=<yes|no|flip|mul{no}>]
        if ( parameter.containsKey("--ctf")) {
            ctf = parameter.get("--ctf");
        }
        // private float kv;//=$KV \          # [kv=<acceleration voltage(in kV){300.}>]
        if (parameter.containsKey("--kv")) {
            kv = Float.parseFloat(parameter.get("--kv"));
        }
        // private float fraca;//=$FRACA  \   # [fraca=<frac amp contrast{0.07}>]
        if (parameter.containsKey("--fraca")) {
            fraca = Float.parseFloat(parameter.get("--fraca"));
        }
        // private float cs;//=$CS \          # [cs=<spherical aberration constant(in mm){2.7}>]
        if (parameter.containsKey("--cs")) {
            cs = Float.parseFloat(parameter.get("--cs"));
        }
        // private String deftab;//=$DEFTAB \ # [deftab=<text file defocus values>]
        if ( parameter.containsKey("--deftab")) {
            deftab = parameter.get("--deftab");
        }
        // private int state;//=$STATE \      # [state=<state to reconstruct{all}>]
        if ( parameter.containsKey("--state")) {
            state = Integer.parseInt(parameter.get("--state"));
        }
        // private float alpha;//=$ALPHA \    # [alpha=<oversampling ratio or padding factor{2.}>]
        if (parameter.containsKey("--alpha")) {
            alpha = Float.parseFloat(parameter.get("--alpha"));
        }
        // private float mw;//=$MW \          # [mw=<molecular weight(in kD)>]
        if (parameter.containsKey("--mw")) {
            mw = Float.parseFloat(parameter.get("--mw"));
        }
        // private float amsklp;//=$AMSKLP \  # [amsklp=<low-pass limit(in A){20}>]
        if (parameter.containsKey("--amsklp")) {
            amsklp = Float.parseFloat(parameter.get("--amsklp"));
        }
        // private int edge;//=$EDGE  \       # [edge=<edge size for softening molecular envelope(in pixels){3}>]
        if ( parameter.containsKey("--edge")) {
            edge = Integer.parseInt(parameter.get("--edge"));
        }
        // private float dens;//=$DENS \      # [dens=<density(e.g.9.368 Da/A3 4 gold clusters){0.}>]
        if (parameter.containsKey("--dens")) {
            dens = Float.parseFloat(parameter.get("--dens"));
        }
        // private int box;//=$BOX \          # [box=<image size(in pixels)>]
        if ( parameter.containsKey("--box")) {
            box = Integer.parseInt(parameter.get("--box"));
        }
        // private float inner;//=$INNER \    # [inner=<inner mask radius(in pixels)>]
        if (parameter.containsKey("--inner")) {
            inner = Float.parseFloat(parameter.get("--inner"));
        }
        // private float width;//=$WIDTH      # [width=<pixels falloff inner mask{10}>]
        if (parameter.containsKey("--width")) {
            width = Float.parseFloat(parameter.get("--width"));
        }

    }
}
