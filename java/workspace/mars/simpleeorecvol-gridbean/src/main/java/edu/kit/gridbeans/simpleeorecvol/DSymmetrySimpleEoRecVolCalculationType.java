/**
 * 
 */
package edu.kit.gridbeans.simpleeorecvol;

import java.util.Vector;

/**
 * @author Frederic Bonnet
 * @Date 1st of December 2015
 *
 */
public enum DSymmetrySimpleEoRecVolCalculationType {
    d1,d2, d3, d4, d5, d6, d7, d8, d9;
    public static Vector<String> getAllTypes() {
        DSymmetrySimpleEoRecVolCalculationType[] values = DSymmetrySimpleEoRecVolCalculationType.values();
        Vector<String> stringValues = new Vector<String>();

        for (int i = 0; i < values.length; i++) {
            stringValues.add(values[i].name());
        }

        return stringValues;
    }

}
