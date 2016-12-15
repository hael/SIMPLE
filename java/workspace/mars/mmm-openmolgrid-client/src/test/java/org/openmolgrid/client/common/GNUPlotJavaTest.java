/**
 * 
 */
package org.openmolgrid.client.common;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.Map;

import junit.framework.Assert;

import org.junit.Test;
import org.yaml.snakeyaml.Yaml;

import com.panayotis.gnuplot.JavaPlot;
import com.panayotis.gnuplot.plot.Graph3D;

/**
 * @author Frederic Bonnet
 *
 */
public class GNUPlotJavaTest {
	private Yaml yaml = new Yaml();
	//input test files 
	private static final String INPUT_TIME_YAML_FILE = "time.yaml";
	private static final String INPUT_LOG_AGBULK_YAML_FILE = "log_AgBulk.yaml";
	private static final String INPUT_LOG_MULTIPT_YAML_FILE = "log_multiPt.yaml";
	private static final String INPUT_LOG_SINGLEPT_YAML_FILE = "log_singlePt.yaml";
	//output test files
	private static final String CML_OUTPUT_GEOMETRY_MULTIPT_YAML_FILE = "geometry_multiPt.cml";

	/**
	 * test the Yaml file
	 */
	@Test
	public void testYamlfile() {
		Reader reader = null;

		try {
			String 
			path = "src" + File.separator
					+ "test" + File.separator
					+ "resources" + File.separator
					+ "yaml" + File.separator
					+ INPUT_TIME_YAML_FILE;
			
			File testFile = new File(path);
	
			reader = new FileReader(testFile);
			Assert.assertTrue(testFile.exists());
			Iterable<Object> allYaml = yaml.loadAll(new FileInputStream(testFile));

			Map<?, ?> yamlMap = null;
			Object lastYaml = null;
			for (Object currentYaml : allYaml) {
				//System.out.println(currentYaml);
				lastYaml = currentYaml;
			}

			yamlMap = (Map<?, ?>) lastYaml;
			//System.out.println(yaml.dump(yamlMap));

			Assert.assertNotNull(yamlMap);
			
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

	}
	/**
	 * test if the Yaml file exist
	 */
	@Test
	public void testYamlExist() {
		Reader reader = null;
	
		try {
			File testFile = new File("src" + File.separator + "test" + File.separator + "resources" + File.separator
					+ "yaml" + File.separator + INPUT_TIME_YAML_FILE);
	
			reader = new FileReader(testFile);
			Assert.assertTrue(testFile.exists());
				
		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	 }
	/**
	 * test simple sin(x) function plot
	 */
	//@Test
	public void testSinxJavaPlot() {

		JavaPlot p = new JavaPlot();
        p.addPlot("sin(x)");
        p.plot();

	}
	/**
	 * test if the class JavaPLot can be used in the tester class
	 */
	@Test
	public void testGnuPlot() {
		Reader reader = null;

		try {
			String path = "src" + File.separator
					+ "test" + File.separator
					+ "resources" + File.separator
					+ "yaml" + File.separator
					+ INPUT_TIME_YAML_FILE;

		} catch (Exception e) {
			e.printStackTrace();
			Assert.fail();
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	Graph3D graph3di = new Graph3D();
	
			
}
