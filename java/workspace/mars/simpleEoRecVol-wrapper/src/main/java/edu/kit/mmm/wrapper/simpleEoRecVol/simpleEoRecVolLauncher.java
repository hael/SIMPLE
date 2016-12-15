package edu.kit.mmm.wrapper.simpleEoRecVol;

import org.apache.log4j.Logger;

public class simpleEoRecVolLauncher 
{
  private static Logger log = Logger.getLogger(simpleEoRecVolLauncher.class.getName());
	
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

    SimpleEoRecVolCalculation eoRecVolCalc = new SimpleEoRecVolCalculation(args);
    eoRecVolCalc.startWrapper();
		
  }
}
