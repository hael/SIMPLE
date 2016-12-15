package edu.kit.mmm.wrapper.simple;

import org.apache.log4j.Logger;

public class simpleLauncher 
{
  private static Logger log = Logger.getLogger(simpleLauncher.class.getName());
	
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

    simpleCalculation simpleCalc = new simpleCalculation(args);
    simpleCalc.startWrapper();
		
  }
}
