/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.process;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import org.apache.log4j.Logger;
import org.openmolgrid.util.CUtil;

/**
 * ProcessWrapper is a utility for running external programs. It starts process and waits for the termination while
 * collecting output from standard output and error.
 * 
 * This code was inspired by the following article: http://www.javaworld.com/javaworld/jw-12-2000/jw-1229-traps_p.html
 * 
 * @author Sulev Sild
 */
public class ProcessWrapper {

	/** Logger */
	private static Logger log = Logger.getLogger(ProcessWrapper.class.getName());

	/** Whether the std out and err should be printes by System.out.println */
	private boolean quiet;

	/** The working directory of the process */
	private File workingDirectory;

	/** The process instance created by this class */
	private Process process;

	/** StdOut buffer */
	private StringBuffer procStdout;

	/** StdErr buffer */
	private StringBuffer procStderr;

	/** the execution time of the process */
	private long executionTime;

	/**
	 * Constructor for ProcessWrapper.
	 * 
	 * @param quiet - if false, then stdout/stderr will be printed to console
	 */
	public ProcessWrapper(boolean quiet) {
		procStdout = new StringBuffer();
		procStderr = new StringBuffer();
		executionTime = 0;
		this.quiet = quiet;
	}

	/**
	 * Default constructor without printing the stdout and stderr to the logger
	 */
	public ProcessWrapper() {
		this(true);
	}

	/**
	 * Constructor for ProcessWrapper
	 * 
	 * @param workingDir the workingDirectory for the process
	 */
	public ProcessWrapper(File workingDir) {
		this(true);
		this.workingDirectory = workingDir;
	}

	/**
	 * Executes a given command string. This method saves output from the standard output and error. When command is
	 * finished, the return error from the command is returned. The output from the command is accessible with getStdout
	 * and getStderr methods. The command will be executed in the current working directory where the JVM was started,
	 * unless an alternative directory has been created with createWorkingDirectory() method.
	 * 
	 * @param command String with command to be executed
	 * @return error code from the command
	 * @throws IOException if execution fails (e.g. program not found)
	 */
	public int execute(String command) throws IOException {
		return execute(command, workingDirectory);
	}

	/**
	 * Executes a given command string. This method saves output from the standard output and error. When command is
	 * finished, the return error from the command is returned. The output from the command is accessible with getStdout
	 * and getStderr methods.
	 * 
	 * @param command String with command to be executed
	 * @param workDir path to the working directory
	 * @return error code from the command
	 * @throws IOException if execution fails (e.g. program not found)
	 */
	public int execute(String command, File workDir) throws IOException {
		procStdout.delete(0, procStdout.length());
		procStderr.delete(0, procStderr.length());
		executionTime = System.currentTimeMillis();

		try {
			process = Runtime.getRuntime().exec(command, null, workDir);
			InputStream output = process.getInputStream();
			Thread out = new Thread(new StreamReaderRunnable("stdout", output, procStdout, quiet));
			out.start();
			InputStream error = process.getErrorStream();
			Thread err = new Thread(new StreamReaderRunnable("stderr", error, procStderr, quiet));
			err.start();

			// wait for threads to finish
			err.join();
			out.join();
			int retCode = process.waitFor();
			executionTime = System.currentTimeMillis() - executionTime;
			return retCode;
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new IOException("Thread interrupted!");
		}
	}

	/**
	 * Executes a given command string. With this its possible to create pipes in Linux! This method saves output from
	 * the standard output and error. When command is finished, the return error from the command is returned. The
	 * output from the command is accessible with getStdout and getStderr methods. The command will be executed in the
	 * current working directory where the JVM was started, unless an alternative directory has been created with
	 * createWorkingDirectory() method.
	 * 
	 * @param command String with command to be executed
	 * @return error code from the command
	 * @throws IOException if execution fails (e.g. program not found)
	 */
	public int execute(String[] command) throws IOException {
		return execute(command, workingDirectory);
	}

	/**
	 * Executes a given command array. With this its possible to create pipes in Linux! This method saves output from
	 * the standard output and error. When command is finished, the return error from the command is returned. The
	 * output from the command is accessible with getStdout and getStderr methods.
	 * 
	 * @param command String with command to be executed
	 * @param workDir path to the working directory
	 * @return error code from the command
	 * @throws IOException if execution fails (e.g. program not found)
	 */
	public int execute(String[] command, File workDir) throws IOException {
		procStdout.delete(0, procStdout.length());
		procStderr.delete(0, procStderr.length());

		executionTime = System.currentTimeMillis();

		try {
			log.debug("Inside execute");
			log.debug("command: " + command);
			log.debug("workDir: " + workDir);

			// ProcessBuilder procBuilder = new ProcessBuilder(command);
			// process = procBuilder.start();

			process = Runtime.getRuntime().exec(command, null, workDir);

			log.debug("Proc created");
			log.debug(process.toString());

			InputStream output = process.getInputStream();
			StreamReaderRunnable stdStreamRunnable = new StreamReaderRunnable("stdout", output, procStdout, quiet);
			Thread out = new Thread(stdStreamRunnable);
			out.start();
			log.debug("Thread out started");

			InputStream error = process.getErrorStream();
			Thread err = new Thread(new StreamReaderRunnable("stderr", error, procStderr, quiet));
			err.start();
			log.debug("Thread err started");

			// wait for threads to finish
			err.join();
			log.debug("err joined");

			out.join();
			log.debug("out joined");

			log.debug("wait for");
			int retCode = process.waitFor();
			log.debug("wait ends");

			executionTime = System.currentTimeMillis() - executionTime;
			return retCode;
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new IOException("Thread interrupted!");
		}
	}

	/**
	 * Executes a given command array in a new {@link Process}. This method will not wait for the Process to end.
	 * 
	 * @param command String with command to be executed
	 * @throws IOException if execution fails (e.g. program not found)
	 */
	public Process executeAsync(String[] command) throws IOException {
		return executeAsync(command, workingDirectory);
	}

	/**
	 * Executes a given command array in a new {@link Process}. This method will not wait for the Process to end.
	 * 
	 * @param command String with command to be executed
	 * @param workDir path to the working directory
	 * 
	 * @throws IOException if execution fails (e.g. program not found)
	 */
	public Process executeAsync(String[] command, File workDir) throws IOException {
		process = Runtime.getRuntime().exec(command, null, workDir);

		return process;
	}

	/**
	 * Gets data that was printed to the standard output by the process.
	 * 
	 * @return StringBuffer with the contents of the standard output
	 */
	public StringBuffer getStdout() {
		return procStdout;
	}

	/**
	 * Gets data that was printed to the standard error by the process.
	 * 
	 * @return StringBuffer with the contents of the standard error
	 */
	public StringBuffer getStderr() {
		return procStderr;
	}

	/**
	 * Returns time that was spent executing the command.
	 * 
	 * @return execution time in milliseconds
	 */
	public long getExecutionTime() {
		return executionTime;
	}

	/**
	 * @return the process
	 */
	public Process getProcess() {
		return process;
	}

	/**
	 * @param process the process to set
	 */
	public void setProcess(Process process) {
		this.process = process;
	}

	/**
	 * Creates a new working directory. This method causes the subsequent commands to be executed in that directory. The
	 * directory where the new directory will be created is determined by the java.io.tmpdir system property.
	 * 
	 * @return a File object corresponding to this directory.
	 * @throws IOException If directory can not be created
	 */
	public File createWorkingDirectory() throws IOException {
		File tmp = java.io.File.createTempFile("omgtmp", "");
		if (!tmp.delete()) {
			throw new IOException("Unable to delete directory: " + tmp.getAbsolutePath());
		}
		if (!tmp.mkdir()) {
			throw new IOException("Unable to create directory: " + tmp.getAbsolutePath());
		}

		workingDirectory = tmp;
		return workingDirectory;
	}

	/**
	 * Deletes the working directory. After this method commands will be executed in the current working directory where
	 * the JVM was started.
	 * 
	 * @return true if directory was deleted, false otherwise
	 */
	public boolean deleteWorkingDirectory() {
		if (null != workingDirectory) {
			File tmp = workingDirectory;
			workingDirectory = null;
			return CUtil.deleteAll(tmp);
		}

		return false;
	}
}
