/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * Version.java
 *
 * Created on Mar 5, 2012, 11:18:07 PM
 */
package edu.kit.gridbeans.simpleeorecvol.plugin;

/**
 * Class that returns the version numbers of the projects
 * @author Frederic Bonnet
 * @Date 20th of November 2015
 *
 */
public class VersionEoRecVol {

	/**
	 * version for the Eo_RecVol GridBean
	 * @return v1.0:String
	 */
	public static String getVersionEoRecVolGridBean() {
		return "v1.0";
	}
	/**
	 * version for the Eo_RecVol application
	 * @return November 20th of November 2015 - version 1.0:String
	 */
	public static String getVersionEoRecVol() {
		return "May 16th 2015 - version 2.1";
	}
	/**
	 * method to return the version number of 
	 * @return March 7. 2012 - version v1.0a0:String
	 */
	public static String getVersionParsingClasses() {
		return "August 24. 2013 - version v1.0";
	}
	/**
	 * method to return the header message of the parsing class SnakeYaml
	 * @return header message:String
	 */
	public static String getHeaderSnakeYaml() {
		return 	" * Copyright (c) 2008-2012, http://www.snakeyaml.org"+
				" *"+
				" * Licensed under the Apache License, Version 2.0 (the "+"'License'"+");"+
				" * you may not use this file except in compliance with the License."+
				" * You may obtain a copy of the License at"+
				" *"+
				" *     http://www.apache.org/licenses/LICENSE-2.0"+
				" *"+
				" * Unless required by applicable law or agreed to in writing, software"+
				" * distributed under the License is distributed on an 'AS IS' BASIS,"+
				" * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied."+
				" * See the License for the specific language governing permissions and"+
				" * limitations under the License.";
	}
	
	
	
}
