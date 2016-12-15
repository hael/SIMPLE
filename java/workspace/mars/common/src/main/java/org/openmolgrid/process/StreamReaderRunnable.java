package org.openmolgrid.process;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 * Helper class to save output that is produced by the running process. This class can also inform handlers whether a
 * specific String has occured in the {@link InputStream}.
 */
public class StreamReaderRunnable implements Runnable {

	/** Logger */
	private static Logger log = Logger.getLogger(StreamReaderRunnable.class.getName());

	/** The name of the stream */
	private String streamName;

	/** Whether the stream should be traced to the logger */
	private boolean quiet;

	/** The input stream */
	private InputStream is;

	/** The buffer for the reader */
	private StringBuffer saveBuf;

	/** If this String is read the */
	private String notifyString;

	/** List with handler that will be informed if a specific {@link String} has been read. */
	private List<StreamReaderHandler> streamReaderHandlerList;

	/**
	 * Constructor
	 * 
	 * @param name the name of the stream
	 * @param is an {@link InputStream} of the process
	 * @param saveBuf the buffer to copy data from the {@link InputStream}
	 * @param quiet Whether the stream should be traced to the logger
	 */
	public StreamReaderRunnable(String name, InputStream is, StringBuffer saveBuf, boolean quiet) {
		this.streamName = name;
		this.streamReaderHandlerList = new ArrayList<StreamReaderHandler>();
		this.is = is;
		this.saveBuf = saveBuf;
		this.quiet = quiet;
	}

	/**
	 * Constructor
	 * 
	 * @param name the name of the stream
	 * @param is an {@link InputStream} of the process
	 * @param saveBuf the buffer to copy data from the {@link InputStream}
	 * @param quiet Whether the stream should be traced to the logger
	 */
	public StreamReaderRunnable(String name, InputStream is, StringBuffer saveBuf, boolean quiet, String notifyString) {
		this(name, is, saveBuf, quiet);
		this.notifyString = notifyString;
		this.streamReaderHandlerList = new ArrayList<StreamReaderHandler>();
	}

	/**
	 * Add a {@link StreamReaderHandler} to the handler list.
	 * 
	 * @param handler the handler to addd
	 */
	public void addStreamReaderhandler(StreamReaderHandler handler) {
		streamReaderHandlerList.add(handler);
	}

	/**
	 * Removes a {@link StreamReaderHandler} from the handler list.
	 * 
	 * @param handler the handler to remove
	 * @return <code>true</code> if this list contained the specified element
	 */
	public boolean removeStreamReaderhandler(StreamReaderHandler handler) {
		return streamReaderHandlerList.remove(handler);
	}

	/**
	 * Notify all handlers if a specific {@link String} has been read by the Streamreader.
	 */
	private void notifyHandlers() {
		for (StreamReaderHandler handler : streamReaderHandlerList) {
			handler.stringRead();
		}
	}

	/**
	 * Copies the {@link InputStream} from the process' out and err to the instance buffer
	 */
	public void run() {
		BufferedReader in = null;
		try {
			in = new BufferedReader(new InputStreamReader(is));
			String temp = null;

			while ((temp = in.readLine()) != null) {
				saveBuf.append(temp + '\n');
				if (!quiet) {
					log.info(streamName + ": " + temp);
				}

				if (notifyString != null) {
					if (temp.contains(notifyString)) {
						log.info("Found notifyString in InputStream. Notify all registered Handlers.");
						notifyHandlers();
					}
				}
			}
		} catch (IOException e) {
			log.error(e);
		} finally {
			try {
				if (in != null) {
					in.close();
				}

				if (is != null) {
					is.close();
				}
			} catch (IOException e) {
				log.error("An error occurs while closing Reader", e);
			}
		}
	}
}
