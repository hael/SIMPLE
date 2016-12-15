/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common.table;

import org.openmolgrid.qsar.CProperty;
import org.openmolgrid.qsar.CPropertyList;

class ColumnProperty extends Column {
	private final String id;

	public ColumnProperty(StructureList sl, String id, String name) {
		super(sl, (name != null) ? name : id);
		this.id = id;
	}

	public Double get(int row) {
		CPropertyList pl = sl.get(row).getPropertyList();
		if (pl != null) {
			CProperty p = pl.getPropertyById(id);
			if (p != null)
				return p.getValue();
		}
		return null;
	}

}
