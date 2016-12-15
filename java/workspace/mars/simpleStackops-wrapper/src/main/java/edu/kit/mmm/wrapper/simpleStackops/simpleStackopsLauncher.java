package edu.kit.mmm.wrapper.simpleStackops;

import org.apache.log4j.Logger;

public class simpleStackopsLauncher 
{
  private static Logger log = Logger.getLogger(simpleStackopsLauncher.class.getName());
	
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

    SimpleStackopsCalculation stackopsCalc = new SimpleStackopsCalculation(args);
    stackopsCalc.startWrapper();
		
  }
}
