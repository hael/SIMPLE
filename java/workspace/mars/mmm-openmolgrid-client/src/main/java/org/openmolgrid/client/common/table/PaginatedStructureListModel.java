/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.client.common.table;

import org.openmolgrid.model.CStructure;


public class PaginatedStructureListModel extends StructureListModel {

	private static final int SIZE = 30;
	private int currentPage;
	
	public PaginatedStructureListModel(StructureList strList) {
		super(strList);
	}
	
	@Override
	public int getRowCount() {
		int N = getTotalSize();
		int n = N - currentPage*SIZE;
		return n < SIZE ? n : SIZE;
	}
	
	@Override
	public Object getValueAt(int row, int col) {
		row += currentPage*SIZE;
		return super.getValueAt(row, col);
	}
	
	@Override
	public CStructure getStructure(int row) {
		return super.getStructure(currentPage*SIZE+row);
	}

	public int getPageIndex() {
		return currentPage*SIZE + 1;
	}

	public int getTotalSize() {
		return getStructureList().size();
	}

	public int getNumberOfPages() {
		int n = getTotalSize()/SIZE;
		return n;
	}

	public void gotoFirstPage() {
		currentPage = 0;
		fireTableDataChanged();
	}
	
	public boolean hasNextPage() {
		return currentPage < getNumberOfPages();
	}
	
	public boolean hasPreviousPage() {
		return currentPage > 0;
	}
	
	public void gotoPreviousPage() {
		if (hasPreviousPage()) {
			currentPage--;
			fireTableDataChanged();
		} else {
			throw new IndexOutOfBoundsException();
		}

	}
	
	public void gotoNextPage() {
		if (hasNextPage()) {
			currentPage++;
			fireTableDataChanged();
		} else {
			throw new IndexOutOfBoundsException();
		}
	}

	public void gotoLastPage() {
		currentPage = getNumberOfPages();
		fireTableDataChanged();
	}
	
	public void gotoRecord(int n) {
		// TODO: highlight n-th table row
		currentPage = n/SIZE;
		fireTableDataChanged();
	}
	
	@Override
	public void handleStructureAdded(final int row, final CStructure s) {
		if (row > currentPage*SIZE+SIZE) return;
		super.handleStructureAdded(row, s);
	}

}
