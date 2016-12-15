package edu.kit.mmm.wrapper.simpleEoVolAssemble;

import org.apache.log4j.Logger;

public class simpleEoVolAssembleLauncher 
{
  private static Logger log = Logger.getLogger(simpleEoVolAssembleLauncher.class.getName());
	
  public static void main(String[] args) 
  {
		
	    log.info("                                                     ");
	    log.info("          Instance of simpleCalculation              ");
	    log.info("                                                     ");
	    log.info("   ####      #    #    #  #####   #       ######     ");
	    log.info("  #          #    ##  ##  #    #  #       #          ");
	    log.info("   ####      #    # ## #  #    #  #       #####      ");
	    log.info("       #     #    #    #  #####   #       #          ");
	    log.info("  #    #     #    #    #  #       #       #          ");
	    log.info("   ####      #    #    #  #       ######  ######     ");
	    log.info("                                                     ");
	    log.info("By Frederic Bonnet 4th of March 2015                 ");
	    log.info("                                                     ");

    SimpleEoVolAssembleCalculation eoVolAssembleCalc = new SimpleEoVolAssembleCalculation(args);
    eoVolAssembleCalc.startWrapper();
		
  }
}
