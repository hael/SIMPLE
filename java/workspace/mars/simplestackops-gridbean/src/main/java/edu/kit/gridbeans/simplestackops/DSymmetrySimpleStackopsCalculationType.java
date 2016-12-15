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
public enum DSymmetrySimpleStackopsCalculationType {
    d1,d2, d3, d4, d5, d6, d7, d8, d9;
    public static Vector<String> getAllTypes() {
        DSymmetrySimpleStackopsCalculationType[] values = DSymmetrySimpleStackopsCalculationType.values();
        Vector<String> stringValues = new Vector<String>();

        for (int i = 0; i < values.length; i++) {
            stringValues.add(values[i].name());
        }

        return stringValues;
    }

}
