package edu.kit.mmm.wrapper.simplePrime;

import org.apache.log4j.Logger;

public class simplePrimeLauncher 
{
  private static Logger log = Logger.getLogger(simplePrimeLauncher.class.getName());
	
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

    SimplePrimeCalculation simplePrimeCalc = new SimplePrimeCalculation(args);
    simplePrimeCalc.startWrapper();
		
  }
}
