package edu.kit.mmm.wrapper.simplePreProc;

import org.apache.log4j.Logger;

public class simplePreProcLauncher 
{
  private static Logger log = Logger.getLogger(simplePreProcLauncher.class.getName());
	
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

	    simplePreProcCalculation simplePreProcCalc = new simplePreProcCalculation(args);
	    simplePreProcCalc.startWrapper();
		
  }
}
