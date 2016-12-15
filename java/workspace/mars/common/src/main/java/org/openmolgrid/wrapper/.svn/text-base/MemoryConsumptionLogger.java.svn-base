package org.openmolgrid.wrapper;

import org.apache.log4j.Logger;

/**
 * Class which logs the MemoryUsage in a specified interval.
 * 
 * @author Stefan Bozic
 */
public class MemoryConsumptionLogger implements Runnable {

	/** Logger */
	private static Logger log = Logger.getLogger(MemoryConsumptionLogger.class.getName());

	/** Default interval is one minute */
	private long loggingInterval = 60000;

	/**
	 * Constructor
	 */
	public MemoryConsumptionLogger() {
	}

	/**
	 * Constructor
	 * 
	 * @param loggingInterval specifies the logging interval
	 */
	public MemoryConsumptionLogger(long loggingInterval) {
		this.loggingInterval = loggingInterval;
	}

	/**
	 * @return the loggingInterval
	 */
	public long getLoggingInterval() {
		return loggingInterval;
	}

	/**
	 * @param loggingInterval the loggingInterval to set
	 */
	public void setLoggingInterval(long loggingInterval) {
		this.loggingInterval = loggingInterval;
	}

	/** This flag will stop the logging. */
	private boolean stop = false;

	/**
	 * @return the stop
	 */
	public boolean isStop() {
		return stop;
	}

	/**
	 * @param stop the stop to set
	 */
	public void setStop(boolean stop) {
		this.stop = stop;
	}

	/**
	 * @see java.lang.Runnable#run()
	 */
	public void run() {
		int mb = 1024 * 1024;
		
		synchronized (this) {
			while (!stop) {
				Runtime runtime = Runtime.getRuntime();

				log.info("Memory information [MB]:");
				log.info("------------------------------------");
				log.info("Used Memory:" + (runtime.totalMemory() - runtime.freeMemory()) / mb);
				log.info("Free Memory:" + runtime.freeMemory() / mb);
				log.info("Max Memory:" + runtime.maxMemory() / mb);
				log.info("Total Memory:" + runtime.totalMemory() / mb);
				log.info("------------------------------------");

				try {
					wait(loggingInterval);
				} catch (InterruptedException e) {
					log.error("Sleep interrupted:" + e);
				}
			}

			log.info("Stop MemoryLogger.");
		}
	}

}
