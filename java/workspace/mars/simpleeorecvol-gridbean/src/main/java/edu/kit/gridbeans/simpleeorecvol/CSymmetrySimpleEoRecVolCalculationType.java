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
public enum CSymmetrySimpleEoRecVolCalculationType {
    c1,c2, c3, c4, c5, c6, c7, c8, c9;
    public static Vector<String> getAllTypes() {
        CSymmetrySimpleEoRecVolCalculationType[] values = CSymmetrySimpleEoRecVolCalculationType.values();
        Vector<String> stringValues = new Vector<String>();

        for (int i = 0; i < values.length; i++) {
            stringValues.add(values[i].name());
        }
        
        return stringValues;
    }

}
