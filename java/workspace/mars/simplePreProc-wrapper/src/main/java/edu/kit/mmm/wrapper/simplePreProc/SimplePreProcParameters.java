/**
 * 
 */
package edu.kit.mmm.wrapper.simplePreProc;

import java.util.HashMap;
import java.util.logging.Logger;

/**
 * @author Frederic Bonnet
 * @Date 18th of November 2015
 *
 */
public class SimplePreProcParameters {
    /** Logger */
    private static Logger log = Logger.getLogger(SimplePreProcParameters.class.getName());
    private String workingDirectory; //--workingDirectory=`pwd` \
    private String simple_data_path; //=$SIMPLE_DATA_PATH \  #='path to data';
    private int nframes;             //=$NFRAMES \           #= 7;
    private int nconfig;             //=$NCONFIG \           #= 452; #not used for checking only, real nconfig extracted
    private double sph_abe;          //=$SPH_ABE \           #= 2.7; #sphererical aberration
    private double amp_con;          //=$AMP_CON \           #= 0.07;#Amplitude constrast
    private int sze_pwr_spc;         //=$SZE_PWR_SPC \       #= 512; #Size of prower spectrum 
    private int data_rot;            //=$DATA_ROT \          #= -90; #Rotation angle for the Size of prower spectrum 
    private String target_file;      //=$TARGET_FILE \       #= 'target.txt';#Name of output file list
    private int create_dir;          //=$CREATE_DIR \        #= 0;  #set to 1 if they not already created otherwise set to 0 
    private String file_pre_handler; //=$FILE_PRE_HANDLER \  # = "preprocess";
    private String file_organise_dir;//=$FILE_ORGANISE_DIR \ #= "or";
    private int unblur;              //=$UNBLUR \            #= 1;
    private int ctffind;             //=$CTFFIND \           #= 1;
    private String unblur_dir;       //=$UNBLUR_DIR \        #= "/opt/Devel_tools/unblur_1.0.2/bin";
    private String ctffind_dir;      //=$CTFFIND_DIR \       #= "/opt/Devel_tools/ctffind-4.0.16/build";
    private String file_pst_handler; //=$FILE_PST_HANDLER \  #= "extract_data";
    private int have_boxfile;        //=$HAVE_BOXFILE \      #= "0";
    private int verbose;             //=$VERBOSE \           #= "1";
    private int delete_dir;          //=$DELETE_DIR          #= 0;
    
    /**
     * @return the workingDirectory
     */
    public String getWorkingDirectory() {
        return workingDirectory;
    }

    /**
     * @param workingDirectory the workingDirectory to set
     */
    public void setWorkingDirectory(String workingDirectory) {
        this.workingDirectory = workingDirectory;
    }

    /**
     * @return the simple_data_path
     */
    public String getSimple_data_path() {
        return simple_data_path;
    }

    /**
     * @param simple_data_path the simple_data_path to set
     */
    public void setSimple_data_path(String simple_data_path) {
        this.simple_data_path = simple_data_path;
    }

    /**
     * @return the nframes
     */
    public int getNframes() {
        return nframes;
    }

    /**
     * @param nframes the nframes to set
     */
    public void setNframes(int nframes) {
        this.nframes = nframes;
    }

    /**
     * @return the nconfig
     */
    public int getNconfig() {
        return nconfig;
    }

    /**
     * @param nconfig the nconfig to set
     */
    public void setNconfig(int nconfig) {
        this.nconfig = nconfig;
    }

    /**
     * @return the sph_abe
     */
    public double getSph_abe() {
        return sph_abe;
    }

    /**
     * @param sph_abe the sph_abe to set
     */
    public void setSph_abe(double sph_abe) {
        this.sph_abe = sph_abe;
    }

    /**
     * @return the amp_con
     */
    public double getAmp_con() {
        return amp_con;
    }

    /**
     * @param amp_con the amp_con to set
     */
    public void setAmp_con(double amp_con) {
        this.amp_con = amp_con;
    }

    /**
     * @return the sze_pwr_spc
     */
    public int getSze_pwr_spc() {
        return sze_pwr_spc;
    }

    /**
     * @param sze_pwr_spc the sze_pwr_spc to set
     */
    public void setSze_pwr_spc(int sze_pwr_spc) {
        this.sze_pwr_spc = sze_pwr_spc;
    }

    /**
     * @return the data_rot
     */
    public int getData_rot() {
        return data_rot;
    }

    /**
     * @param data_rot the data_rot to set
     */
    public void setData_rot(int data_rot) {
        this.data_rot = data_rot;
    }

    /**
     * @return the target_file
     */
    public String getTarget_file() {
        return target_file;
    }

    /**
     * @param target_file the target_file to set
     */
    public void setTarget_file(String target_file) {
        this.target_file = target_file;
    }

    /**
     * @return the create_dir
     */
    public int getCreate_dir() {
        return create_dir;
    }

    /**
     * @param create_dir the create_dir to set
     */
    public void setCreate_dir(int create_dir) {
        this.create_dir = create_dir;
    }

    /**
     * @return the file_pre_handler
     */
    public String getFile_pre_handler() {
        return file_pre_handler;
    }

    /**
     * @param file_pre_handler the file_pre_handler to set
     */
    public void setFile_pre_handler(String file_pre_handler) {
        this.file_pre_handler = file_pre_handler;
    }

    /**
     * @return the file_organise_dir
     */
    public String getFile_organise_dir() {
        return file_organise_dir;
    }

    /**
     * @param file_organise_dir the file_organise_dir to set
     */
    public void setFile_organise_dir(String file_organise_dir) {
        this.file_organise_dir = file_organise_dir;
    }

    /**
     * @return the unblur
     */
    public int getUnblur() {
        return unblur;
    }

    /**
     * @param unblur the unblur to set
     */
    public void setUnblur(int unblur) {
        this.unblur = unblur;
    }

    /**
     * @return the ctffind
     */
    public int getCtffind() {
        return ctffind;
    }

    /**
     * @param ctffind the ctffind to set
     */
    public void setCtffind(int ctffind) {
        this.ctffind = ctffind;
    }

    /**
     * @return the unblur_dir
     */
    public String getUnblur_dir() {
        return unblur_dir;
    }

    /**
     * @param unblur_dir the unblur_dir to set
     */
    public void setUnblur_dir(String unblur_dir) {
        this.unblur_dir = unblur_dir;
    }

    /**
     * @return the ctffind_dir
     */
    public String getCtffind_dir() {
        return ctffind_dir;
    }

    /**
     * @param ctffind_dir the ctffind_dir to set
     */
    public void setCtffind_dir(String ctffind_dir) {
        this.ctffind_dir = ctffind_dir;
    }

    /**
     * @return the file_pst_handler
     */
    public String getFile_pst_handler() {
        return file_pst_handler;
    }

    /**
     * @param file_pst_handler the file_pst_handler to set
     */
    public void setFile_pst_handler(String file_pst_handler) {
        this.file_pst_handler = file_pst_handler;
    }

    /**
     * @return the have_boxfile
     */
    public int getHave_boxfile() {
        return have_boxfile;
    }

    /**
     * @param have_boxfile the have_boxfile to set
     */
    public void setHave_boxfile(int have_boxfile) {
        this.have_boxfile = have_boxfile;
    }

    /**
     * @return the verbose
     */
    public int getVerbose() {
        return verbose;
    }

    /**
     * @param verbose the verbose to set
     */
    public void setVerbose(int verbose) {
        this.verbose = verbose;
    }

    /**
     * @return the delete_dir
     */
    public int getDelete_dir() {
        return delete_dir;
    }

    /**
     * @param delete_dir the delete_dir to set
     */
    public void setDelete_dir(int delete_dir) {
        this.delete_dir = delete_dir;
    }

    /**
     * Maps a {@link HashMap} which contains the SimplePreProc-Wrapper start parameters
     * to an instance of this class.
     * 
     * @param parameter
     *            the HashMap which contains the SimplePreProc-Wrapper start parameters
     */
    public void mapParameter(HashMap<String, String> parameter) {
        log.info("Map Wrapper parameter to SimplePreProcParameters");

        //private String simple_data_path;
        if ( parameter.containsKey("--simple_data_path")) {
            simple_data_path = parameter.get("--simple_data_path");
        }
        //private int nframes;
        if ( parameter.containsKey("--nframes")) {
            nframes = Integer.parseInt(parameter.get("--nframes"));
        }
        //private int nconfig;
        if ( parameter.containsKey("--nconfig")) {
            nconfig = Integer.parseInt(parameter.get("--nconfig"));
        }
        //private double sph_abe;
        if ( parameter.containsKey("--sph_abe")) {
            sph_abe = Double.parseDouble(parameter.get("--sph_abe"));
        }
        //private double amp_con;
        if ( parameter.containsKey("--amp_con")) {
            amp_con = Double.parseDouble(parameter.get("--amp_con"));
        }
        //private int sze_pwr_spc;
        if ( parameter.containsKey("--sze_pwr_spc")) {
            sze_pwr_spc = Integer.parseInt(parameter.get("--sze_pwr_spc"));
        }
        //private int data_rot;
        if ( parameter.containsKey("--data_rot")) {
            data_rot = Integer.parseInt(parameter.get("--data_rot"));
        }
        //private String target_file;
        if ( parameter.containsKey("--target_file")) {
            target_file = parameter.get("--target_file");
        }
        //private int create_dir;
        if ( parameter.containsKey("--create_dir")) {
            create_dir = Integer.parseInt(parameter.get("--create_dir"));
        }
        //private String file_pre_handler;
        if ( parameter.containsKey("--file_pre_handler")) {
            file_pre_handler = parameter.get("--file_pre_handler");
        }
        //private String file_organise_dir;
        if ( parameter.containsKey("--file_organise_dir")) {
            file_organise_dir = parameter.get("--file_organise_dir");
        }
        //private int unblur;
        if ( parameter.containsKey("--unblur")) {
            unblur = Integer.parseInt(parameter.get("--unblur"));
        }
        //private int ctffind;
        if ( parameter.containsKey("--ctffind")) {
            ctffind = Integer.parseInt(parameter.get("--ctffind"));
        }
        //private String unblur_dir;
        if ( parameter.containsKey("--unblur_dir")) {
            unblur_dir = parameter.get("--unblur_dir");
        }
        //private String ctffind_dir;
        if ( parameter.containsKey("--ctffind_dir")) {
            ctffind_dir = parameter.get("--ctffind_dir");
        }
        //private String file_pst_handler;
        if ( parameter.containsKey("--file_pst_handler")) {
            file_pst_handler = parameter.get("--file_pst_handler");
        }
        //private int have_boxfile;
        if ( parameter.containsKey("--have_boxfile")) {
            have_boxfile = Integer.parseInt(parameter.get("--have_boxfile"));
        }
        //private int verbose;
        if ( parameter.containsKey("--verbose")) {
            verbose = Integer.parseInt(parameter.get("--verbose"));
        }
        //private int delete_dir;
        if ( parameter.containsKey("--delete_dir")) {
            delete_dir = Integer.parseInt(parameter.get("--delete_dir"));
        }
        
        // working environment parameters value mapping
        //--workingDirectory=`pwd` \

        if (parameter.containsKey("--workingDirectory")) {
            workingDirectory = parameter.get("--workingDirectory");
        }
 
        
    }

}
