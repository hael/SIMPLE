package org.openmolgrid.wrapper;

import java.util.Date;
import java.util.HashMap;

import org.apache.log4j.Logger;

/**
 * Base class for all application wrapper. It provides some common functionality like printing application parameters, provide a shutdownHook
 * and a logging mechanism for runtime and memory consumption.
 * 
 * @author Stefan Bozic
 */
public abstract class ApplicationWrapper {

	/** Logger */
	private static Logger log = Logger.getLogger(ApplicationWrapper.class.getName());

	/** List with the application parameter */
	protected HashMap<String, String> javaParameter = new HashMap<String, String>();

	/**
	 * Constructor
	 */
	public ApplicationWrapper() {
		
		//Add a shutdownhook that can react if the user pressed ctrl-c or the system send a TERM signal
		Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
		    public void run() {
		        log.warn("ShutdownHook triggered! The ApplicationWrapper receives a term signal! Do some cleanup to avoid data lost!");
		    }
		}));
	}
	
	/**
	 * Starts the ApplicationWrapper which consists of 3 step lifecycle phase
	 * Phase 1: preprocessing
	 * Phase 2: applciation execution
	 * Phase 3: postprocessing
	 * 
	 * The Wrapper starts also an instance of {@link MemoryConsumptionLogger}
	 * to monitor the memory consumption.
	 */
	public void startWrapper() {
		Thread memoryConsumptionThread = null;
		MemoryConsumptionLogger memoryLogger = new MemoryConsumptionLogger(30000);
		
		try {
			//Starts the MemoryConsumptionLogger in a thread
			memoryConsumptionThread = new Thread(memoryLogger);
			memoryConsumptionThread.start();
			log.info("Start memory Logging");
			
			Date startTime = new Date();
			log.info("Application starts at:" + startTime.toString());

			// Execute the abstract methods that are implemented by every subclass of ApplicationWrapper
			preprocess();
			execute();
			postprocess();
			
			Date stopTime = new Date();
			log.info("Application finish at: " + stopTime.toString());

			long elapsedTime = stopTime.getTime() - startTime.getTime();
			log.info("The application wrapper runs " + getElapsedTimeHoursMinutesSecondsString(elapsedTime));
		} catch (Throwable t) {
			log.error("System.exit(1)! An error occurs during application execution. ", t);
			System.exit(1);
		} finally {
			log.info("Finally cleanup application.");
			if (memoryLogger != null && memoryConsumptionThread != null && memoryConsumptionThread.isAlive()) {				
				memoryLogger.setStop(true);
				
				synchronized(memoryLogger) {
					memoryLogger.notify();
				}
			}
		}
	}

	/**
	 * Starts the proprocess. Typical usecases in this lifecycle phase are:
	 * <ul> 
	 * <li>validating application parameters
	 * <li>obtaining application dta from external data sources
	 * <li>generation of application specific input data
	 * </ul>
	 * 
	 * @throws Throwable An Exception that might occurs. 
	 */
	protected abstract void preprocess();
	
	/**
	 * Starts the the application run. Typical usecases  in this lifecycle phase are:
	 * <ul>
	 * <li>Create new processes/threads for the execution of an application
	 * <li>Monitors the stdout and stderr and react on events in the streams
	 * </ul>
	 * 
	 * @throws Throwable An Exception that might occurs. 
	 */
	protected abstract void execute() throws Throwable; 
	
	/**
	 * Starts the postprocessing()
	 * <ul> 
	 * <li>Error handling
	 * <li>Store application data in external data source
	 * <li>Generates an application neutral data format of molecular data
	 * <li>Providing an exit code
	 * </ul>
	 * 
	 * @throws Throwable An Exception that might occurs. 
	 */
	protected abstract void postprocess();
	
	/**
	 * Returns a {@link String} representation of elapsed time in hours/minutes/seconds
	 * @param elapsedTime the elapsedTime in Milliseconds
	 * 
	 * @return String
	 */
	public String getElapsedTimeHoursMinutesSecondsString(long elapsedTime) {
		String format = String.format("%%0%dd", 2);
		elapsedTime = elapsedTime / 1000;
		String seconds = String.format(format, elapsedTime % 60);
		String minutes = String.format(format, (elapsedTime % 3600) / 60);
		String hours = String.format(format, elapsedTime / 3600);
		String time = hours + ":" + minutes + ":" + seconds;
		return time;
	}

	/**
	 * Prints the array of arguments
	 * 
	 * @param args the application arguments
	 */
	protected void printArgs(String[] args) {

		// the args can be null
		if (null != args) {
			log.info("List of program parameters: ");
			for (String arg : args) {
				log.info(arg);
				if (arg.contains("=")) {
					String[] argValue = arg.split("=");

					if (argValue.length == 2) {
						javaParameter.put(argValue[0], argValue[1]);
					}
				}
			}
		}
	}

}
