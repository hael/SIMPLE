/*
 * Copyright 2003-2010 University of Tartu
 */
package org.openmolgrid.model;

/**
 * Exception that is thrown by Chem classes when some problem has occurred.
 */
public class ChemLibException extends Exception {
	public ChemLibException(String message) {
		super(message);
	}

	public ChemLibException(String message, Throwable cause) {
		super(message, cause);
	}
}
