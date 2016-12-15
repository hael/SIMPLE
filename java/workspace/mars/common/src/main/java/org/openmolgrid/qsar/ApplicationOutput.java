/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.qsar;

import java.io.*;

import org.openmolgrid.util.CUtil;

public class ApplicationOutput implements Serializable {
	private String filePath = null;
	private String content = null;

	public ApplicationOutput(String filePath, String content) {
		this.filePath = filePath;
		this.content = content;
	}

	public ApplicationOutput(ApplicationOutput m) {
		filePath = m.filePath;
		content = m.content;
	}

	public String getContent() {
		if (content == null) {
			return CUtil.file2String(new File(filePath));
		} else
			return content;
	}

	public void writeCDR(PrintWriter pw) {
		BufferedReader reader;
		if (content == null) {
			try {
				reader = new BufferedReader(new FileReader(filePath));
			} catch (FileNotFoundException e) {
				throw new RuntimeException("Can't open: " + filePath, e);
			}
		} else {
			reader = new BufferedReader(new StringReader(content));
		}

		try {
			pw.write("<![CDATA[");
			char[] buf = new char[8 * 1024];
			int n;
			while ((n = reader.read(buf)) != -1) {
				pw.write(buf, 0, n);
			}
			pw.write("]]>");
		} catch (IOException e) {
			throw new RuntimeException("Can't read from: " + filePath, e);
		} finally {
			try {
				reader.close();
			} catch (IOException e) {
				throw new RuntimeException("Can't close: " + filePath, e);
			}
		}
	}
}
