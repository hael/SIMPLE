/**
 * 
 */
package edu.kit.gridbeans.simpleeovolassemble;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Frederic Bonnet
 * @Date 3rd of December 2015
 *
 */
public enum NthrSimpleEoVolAssembleCalculationType {
//    [1,2, 3, 4, 5, 6, 7, 8, 9,10];
    singlet(0), doublet(1), triplet(2), quartet(3), pentet(4), sextet(5), septet(6), octet(7), nontet(8);
    private int numberOfThreads;
    
    NthrSimpleEoVolAssembleCalculationType(int numberofthreads){
        this.numberOfThreads = numberofthreads;
    }
    public static List<Integer> getIntegerList() {
        ArrayList<Integer> retVal = new ArrayList<Integer>();
        for (NthrSimpleEoVolAssembleCalculationType value : NthrSimpleEoVolAssembleCalculationType.values()) {
            retVal.add(value.getNthr());
        }
        return retVal;
    }
    public static List<String> getStringList(){
        ArrayList<String> retVal = new ArrayList<String>();

        for (NthrSimpleEoVolAssembleCalculationType value : NthrSimpleEoVolAssembleCalculationType.values()) {
            retVal.add(value.name());
        }
        
        return retVal;
    }
    
    public static List<String> getEvenStringList(){
        ArrayList<String> retVal = new ArrayList<String>();

        for (NthrSimpleEoVolAssembleCalculationType value : NthrSimpleEoVolAssembleCalculationType.values()) {
            if(value.numberOfThreads%2 == 0) retVal.add(value.name());
        }
        
        return retVal;
    }
    
    public static List<String> getOddStringList(){
        ArrayList<String> retVal = new ArrayList<String>();

        for (NthrSimpleEoVolAssembleCalculationType value : NthrSimpleEoVolAssembleCalculationType.values()) {
            if(value.numberOfThreads%2 == 1) retVal.add(value.name());
        }
        
        return retVal;
    }
    public int getNthr() {
        return numberOfThreads;
    }
//   public static Vector<String> getAllTypes() {
//        NthrSimpleEoVolAssembleCalculationType[] values = NthrSimpleEoVolAssembleCalculationType.values();
//        Vector<String> stringValues = new Vector<String>();
//
//        for (int i = 0; i < values.length; i++) {
//            stringValues.add(values[i].name());
//        }
//        
//        return stringValues;
//    }

}
