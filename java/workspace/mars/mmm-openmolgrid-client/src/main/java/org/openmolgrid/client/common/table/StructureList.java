/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common.table;

import java.io.IOException;
import java.util.AbstractList;
import java.util.HashMap;
import java.util.concurrent.CopyOnWriteArraySet;

import org.openmolgrid.model.CStructure;
import org.openmolgrid.qsar.CDescriptorMetaMap;
import org.openmolgrid.qsar.CPropertyMetaMap;



/**
 * This class is for working with a structure list.  
 *
 * @author Sulev Sild
 */
public class StructureList extends AbstractList<CStructure> {
	// descriptor/property meta data
	private CDescriptorMetaMap descMeta = new CDescriptorMetaMap();
	private CPropertyMetaMap propMeta = new CPropertyMetaMap();
	
	// index structures by structure ids
	private HashMap<String, Integer> idMap = new HashMap<String, Integer>();  

	
	public CDescriptorMetaMap getDescriptorMetaData() {
		return descMeta;
	}
	
	public CPropertyMetaMap getPropertyMetaData() {
		return propMeta;
	}

	public void setDescriptorMetaData(CDescriptorMetaMap descMetaMap) {
		descMeta = descMetaMap;
	}

	public void setPropertyMetaData(CPropertyMetaMap propertyMetaMap) {
		propMeta = propertyMetaMap;
	}
	
	public int indexById(String id) {
		Integer i = idMap.get(id);
		if (i != null) {
			if (!id.equals(get(i).getId()))
				throw new IllegalStateException();
			return i;
		}
		return -1;
	}

	// Change listener support ------------------------------------------------
	
	public interface Listener {
		void handleStructureAdded(int index, CStructure s);
		void handleStructureRemoved(int index, CStructure s);
		void handleStructuresCleared();
		void handleStructureChanged(int index, CStructure s);
	}
	
	private CopyOnWriteArraySet<Listener> listeners = 
		new CopyOnWriteArraySet<Listener>();
	
	public void registerListener(Listener l) {
		listeners.add(l);
	}
	
	private void fireAddStructure(int i, CStructure s) {
		for (Listener l: listeners) {
			l.handleStructureAdded(i, s);
		}
	}
	
	private void fireRemoveStructure(int i, CStructure s) {
		for (Listener l: listeners) {
			l.handleStructureRemoved(i, s);
		}
	}
	
	private void fireStructuresCleared() {
		for (Listener l: listeners) {
			l.handleStructuresCleared();
		}
	}

	private void fireStructureChanged(int i, CStructure s) {
		for (Listener l: listeners) {
			l.handleStructureChanged(i, s);
		}
	}
	
	// AbstractList implementation --------------------------------------------
	
	private StructureListCache data = new StructureListCache();
	
	@Override
	public void add(int i, CStructure s) {
		try {
			data.cacheAdd(i, s);
			idMap.put(s.getId(), i);
			fireAddStructure(i, s);
		} catch (IOException e) {
			throw new RuntimeException("I/O error in StructureList.add(): "+e.getMessage());
		}
	}
	
	@Override
	public CStructure remove(int i) {
		try {
			CStructure s = data.cacheRemove(i);
			idMap.remove(get(i).getId());
			fireRemoveStructure(i, s);
			return s;
		} catch (IOException e) {
			throw new RuntimeException("I/O error in StructureList.remove(): "+e.getMessage());
		}
	}
	
	@Override
	public CStructure set(int i, CStructure s) {
		try {
			CStructure old = data.cacheSet(i, s);
			idMap.remove(old.getId());
			idMap.put(s.getId(), i);
			fireStructureChanged(i, s);
			return old;
		} catch (IOException e) {
			throw new RuntimeException("I/O error in StructureList.set(): "+e.getMessage());
		}
	}
	
	@Override
	public void clear() {
		try {
			data.cacheClear();
			idMap.clear();
		} catch (IOException e) {
			throw new RuntimeException("I/O error in StructureList.clear(): "+e.getMessage());
		}
		setDescriptorMetaData(new CDescriptorMetaMap());
		setPropertyMetaData(new CPropertyMetaMap());
		fireStructuresCleared();
	}
	
	@Override
	public CStructure get(int i) {
		try {
			return data.cacheGet(i);
		} catch (IOException e) {
			throw new RuntimeException("I/O error in StructureList.get(): "+e.getMessage());
		}
	}

	@Override
	public int size() {
		return data.cacheSize();
	}

	// -----------------------------------------------------------------------

}
