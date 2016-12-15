/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common.table;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.lang.ref.SoftReference;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.openmolgrid.format.slf.StructureReader;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.model.ChemLibException;
import org.openmolgrid.qsar.StructureConsumer;

/**
 * Disk based cache for large SLF files.
 *  
 * @author Sulev Sild
 */
public class StructureListCache {
	
	private List<Entry> index =  Collections.synchronizedList(new ArrayList<Entry>());
	private FileOutputStream fos;
	private long offset = 0L;
	private File file = null;
	
	private static class Entry {
		final SoftReference<CStructure> struct;
		final long offset;
		final int size;
		
		public Entry(CStructure struct, long offset, int size) {
			this.struct = new SoftReference<CStructure>(struct);
			this.offset = offset;
			this.size = size;
		}
	}

	public int cacheSize() {
		return index.size();
	}

	public synchronized void cacheAdd(int i, CStructure s) throws IOException {
		if (fos == null) {
			if (file == null) {
				String tmpDirectory = System.getProperty("openmolgrid.tmpdir");
				if (tmpDirectory == null)
					file = File.createTempFile("_slfcache", ".tmp");
				else
					file = File.createTempFile("_slfcache", ".tmp", new File(tmpDirectory));
				file.deleteOnExit();
			}
			fos = new FileOutputStream(file);
		}
		
		Entry entry = createEntry(s, i<50);
		index.add(i, entry);
	}

	public synchronized CStructure cacheSet(int i, CStructure s) throws IOException {
		CStructure old = cacheGet(i);
		Entry entry = createEntry(s, true);
		index.set(i, entry);
		return old;
	}
	
	public CStructure cacheRemove(int i) throws IOException {
		CStructure old = cacheGet(i);
		index.remove(i);
		return old;
	}
	
	private int writeEntry(CStructure s) throws IOException {
		ByteArrayOutputStream buf = new ByteArrayOutputStream(8*1024);
		PrintWriter pw = new PrintWriter(new OutputStreamWriter(buf));
		s.writeCDR(pw);
		pw.flush();
		
		buf.writeTo(fos);
		fos.flush();
		return buf.size();
	}
	
	private Entry createEntry(CStructure s, boolean memCached) throws IOException {
		int size = writeEntry(s);
		Entry e = new Entry(memCached ? s : null, offset, size);
		offset += size;
		return e;
	}
	
	public void close() throws IOException {
		fos.close();
	}
	
	public synchronized CStructure cacheGet(int i) throws IOException {
		Entry e = index.get(i);
		
		CStructure s = e.struct.get(); 
		if (s != null) return s;

		RandomAccessFile f = new RandomAccessFile(file, "r");
		byte[] b = new byte[e.size];
		try {
			f.seek(e.offset);
			int off = 0;
			while (off < e.size) {
				int n = f.read(b, off, e.size - off);
				if (n == -1) throw new IOException("Unexpected end of cached file");
				off += n;
			}
		} finally {
			f.close();
		}

		final CStructure[] r = new CStructure[1];
		ByteArrayInputStream is = new ByteArrayInputStream(b);
		StructureReader sr = new StructureReader(is, new StructureConsumer() {
			public void startHandling() {}
			public void endHandling() {}
			public void consumeStructure(CStructure s) {
				r[0] = s;
			}
		});
		try {
			sr.provideStructures();
		} catch (ChemLibException exception) {
			throw new IOException(exception);
		}
		
		index.set(i, new Entry(r[0], e.offset, e.size));
		
		return r[0];
		
		// TODO: replace the above with: return StructureReader.asList(is).get(0);
	}

	public synchronized void cacheClear() throws IOException {
		index.clear();
		offset = 0L;
		if (fos != null) {
			fos = new FileOutputStream(file);
		}
	}
	
}
