/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common.table;

import org.openmolgrid.qsar.CDescriptor;
import org.openmolgrid.qsar.CDescriptorList;


class ColumnDescriptor extends Column {
	private final String id;

	public ColumnDescriptor(StructureList sl, String id, String name) {
		super(sl, (name != null) ? name : id);
		this.id = id;
	}

	public Double get(int row) {
		CDescriptorList pl = sl.get(row).getDescriptorList();
		if (pl != null) {
			CDescriptor p = pl.getDescriptorById(id);
			if (p != null)
				return p.getValue();
		}
		return Double.NaN;
	}
}
