/**
 * 
 */
package edu.kit.mmm.wrapper.simpleEoVolAssemble;

import java.util.HashMap;
import java.util.logging.Logger;

/**
 * @author Frederic Bonnet
 * @Date 18th of November 2015
 *
 */
public class SimpleEoVolAssembleParameters {
    /** Logger */
    private static Logger log = Logger.getLogger(SimpleEoVolAssembleParameters.class.getName());

    private String stk;//=$STK \      #stk=<ptcls.ext>
    private int npart;//t=$NPART \     #npart=<number of partitions to assemble{1}>
    private int nstates;//=$NSTATES \ #nstates=<nr of states{1}> 
    private int nthr;//=$NTHR \       #[nthr=<nr of openMP threads{1}>]
    private float msk;//=$MSK \       #msk=<mask radius(in pixels)> smpd=<sampling distance(in A){0}>
    private float mw;//=$MW \         #[mw=<molecular weight(in kD){0.0}>]
    private float lpstop;//=$LPSTOP \ #[lpstop=<stay at this low-pass limit(in A){7.0}>]
    private float inner;//=$INNER \   #[inner=<inner mask radius(in pixels){0.0}>]
    private float width;//=$WIDTH     #[width=<pixels falloff inner mask{10.0}>]'
    
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
     * @return the npart
     */
    public int getNpart() {
        return npart;
    }

    /**
     * @param npart the npart to set
     */
    public void setNpart(int npart) {
        this.npart = npart;
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
        log.info("Map Wrapper parameter to SimpleEoVolAssembleParameters");
        
        //private String stk;//=$STK \      #stk=<ptcls.ext>
        if ( parameter.containsKey("--stk")) {
            stk = parameter.get("--stk");
        }
        //private int npart;//t=$NPART \     #npart=<number of partitions to assemble{1}>
        if ( parameter.containsKey("--npart")) {
            npart = Integer.parseInt(parameter.get("--npart"));
        }
        //private int nstates;//=$NSTATES \ #nstates=<nr of states{1}> 
        if ( parameter.containsKey("--nstates")) {
            nstates = Integer.parseInt(parameter.get("--nstates"));
        }
        //private int nthr;//=$NTHR \       #[nthr=<nr of openMP threads{1}>]
        if ( parameter.containsKey("--nthr")) {
            nthr = Integer.parseInt(parameter.get("--nthr"));
        }
        //private float msk;//=$MSK \       #msk=<mask radius(in pixels)> smpd=<sampling distance(in A){0}>
        if (parameter.containsKey("--msk")) {
            msk = Float.parseFloat(parameter.get("--msk"));
        }
        //private float mw;//=$MW \         #[mw=<molecular weight(in kD){0.0}>]
        if (parameter.containsKey("--mw")) {
            mw = Float.parseFloat(parameter.get("--mw"));
        }
        //private float lpstop;//=$LPSTOP \ #[lpstop=<stay at this low-pass limit(in A){7.0}>]
        if (parameter.containsKey("--lpstop")) {
            lpstop = Float.parseFloat(parameter.get("--lpstop"));
        }
        //private float inner;//=$INNER \   #[inner=<inner mask radius(in pixels){0.0}>]
        if (parameter.containsKey("--inner")) {
            inner = Float.parseFloat(parameter.get("--inner"));
        }
        //private float width;//=$WIDTH     #[width=<pixels falloff inner mask{10.0}>]'
        if (parameter.containsKey("--width")) {
            width = Float.parseFloat(parameter.get("--width"));
        }

    }
}
