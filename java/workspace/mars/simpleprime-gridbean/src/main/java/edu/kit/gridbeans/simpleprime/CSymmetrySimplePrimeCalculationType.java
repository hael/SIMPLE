/**
 * 
 */
package edu.kit.gridbeans.simpleprime;

import java.util.Vector;

/**
 * @author Frederic Bonnet
 * 30th of November 2015
 *
 */
public enum CSymmetrySimplePrimeCalculationType {
    c1,c2, c3, c4, c5, c6, c7, c8, c9;
    public static Vector<String> getAllTypes() {
        CSymmetrySimplePrimeCalculationType[] values = CSymmetrySimplePrimeCalculationType.values();
        Vector<String> stringValues = new Vector<String>();

        for (int i = 0; i < values.length; i++) {
            stringValues.add(values[i].name());
        }
        
        return stringValues;
    }

}
