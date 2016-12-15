package org.openmolgrid.format.xyz;



import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;

import org.openmolgrid.model.CAtom;
import org.openmolgrid.model.CStructure;
import org.openmolgrid.util.AtomDatabase;

public class XyzReader {
	
	
	public static ArrayList<AtomDatabase> readXyzToAtoms (String input) throws NumberFormatException, Exception {
		ArrayList<AtomDatabase> structure = new ArrayList<AtomDatabase>();
		String[] parsedText = input.split("\n");
		int atomsNumber = 0;
		String[] atomsNumberString = parsedText[0].split("[ ]+");
		try {
			atomsNumber = Integer.parseInt(atomsNumberString[0]);
		} catch (NumberFormatException e){ 
			NumberFormatException e1 = new NumberFormatException("XyzReader.readXyzToAtoms: Can not parse \"" + parsedText[0] + "\" as Integer");
			throw e1;
		}
		int elementEntryIndex = 0;
		for (int i = 2; i < atomsNumber+2; i++) {
			String[] fields = parsedText[i].split("[ ]+");
			if (fields[0].length()== 0) elementEntryIndex = 1; else elementEntryIndex =	0;
			String element = fields[elementEntryIndex];
			/**make first letter to uppercase*/
			char[] tempCharArray = element.toCharArray();
			tempCharArray[0] = Character.toUpperCase(tempCharArray[0]);
			element = new String(tempCharArray);
			
			try{
				AtomDatabase at = AtomDatabase.valueOf(element); 
				structure.add(at);
			} catch(Exception e){ 
				Exception e1 = new Exception("Element Symbol "+ element +" does not exist");
				throw e1;
			}
		} 
		return structure;
	}
	
	public static CStructure readXyz (String input) throws NumberFormatException{
		CStructure structure = new CStructure();
		String[] parsedText = input.split("\n");
		int atomsNumber = 0;
		String[] atomsNumberString = parsedText[0].split("[ ]+");
		try {
			atomsNumber = Integer.parseInt(atomsNumberString[0]);
		} catch (NumberFormatException e){ 
			NumberFormatException e1 = new NumberFormatException("XyzReader.readXyz: Can not parse \"" + parsedText[0] + "\" as Integer");
			throw e1;
		}
		int elementEntryIndex = 0;
		for (int i = 2; i < atomsNumber+2; i++) {
			String[] fields = parsedText[i].split("[ ]+");
			if (fields[0].length()== 0) elementEntryIndex = 1; else elementEntryIndex =	0;
			String element = fields[elementEntryIndex];
			double x=0, y=0, z=0;
			try {
				x = Double.parseDouble(fields[elementEntryIndex+1]);
				y = Double.parseDouble(fields[elementEntryIndex+2]);
				z = Double.parseDouble(fields[elementEntryIndex+3]);
			} catch (NumberFormatException e){ 
				NumberFormatException e1 = new NumberFormatException("XyzReader.readXyz: Can not parse \"" + parsedText[i] + "\" as Coordinates");
				throw e1;
			}
			
			/**make first letter to uppercase*/
			char[] tempCharArray = element.toCharArray();
			tempCharArray[0] = Character.toUpperCase(tempCharArray[0]);
			element = new String(tempCharArray);
			
			CAtom  cAt = new CAtom(element, x, y, z, i-1);
			if (fields.length > elementEntryIndex + 4){
				try{
				double charge = Integer.parseInt(fields[elementEntryIndex+4]);
				cAt.setPartialCharge(charge);
				}  catch (NumberFormatException e){ 
					/** Ignore the failure, because this feature is not officiall xyz feature*/
				}
				
			}
			structure.addAtom(cAt);
		} 
		return structure;
	}
	
	public static CStructure readXyzFile(File file) throws Exception{
		try{
			FileReader ifi = new FileReader(file);
			BufferedReader reader = new BufferedReader(ifi);

			String sContent = new String();

			for (String line = reader.readLine(); line != null; line= reader.readLine()){
				sContent = sContent.concat(line);
				sContent = sContent.concat("\n");
			}
			reader.close();

			return XyzReader.readXyz(sContent);

		} catch (Exception e){
			Exception e1 = new Exception("Could not read the file" + file.getAbsolutePath());
			e1.setStackTrace(e.getStackTrace());
			throw e1;
		} 
	}
}
