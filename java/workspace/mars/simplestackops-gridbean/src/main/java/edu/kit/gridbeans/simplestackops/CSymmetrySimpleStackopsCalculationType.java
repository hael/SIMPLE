/**
 * 
 */
package edu.kit.gridbeans.simplestackops;

import java.util.Vector;

/**
 * @author Frederic Bonnet
 * @Date 30th of November 2015
 *
 */
public enum CSymmetrySimpleStackopsCalculationType {
    c1,c2, c3, c4, c5, c6, c7, c8, c9;
    public static Vector<String> getAllTypes() {
        CSymmetrySimpleStackopsCalculationType[] values = CSymmetrySimpleStackopsCalculationType.values();
        Vector<String> stringValues = new Vector<String>();

        for (int i = 0; i < values.length; i++) {
            stringValues.add(values[i].name());
        }
        
        return stringValues;
    }

}
