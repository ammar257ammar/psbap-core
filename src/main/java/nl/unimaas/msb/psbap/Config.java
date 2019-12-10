/**
* binding Pocket's SNPS effect on Binding Affinity Project (PSBAP) 
* 
*Copyright (C) 2019  Ammar Ammar <ammar257ammar@gmail.com>
*
*This program is free software: you can redistribute it and/or modify
*it under the terms of the GNU Affero General Public License as published by
*the Free Software Foundation, either version 3 of the License, or
*(at your option) any later version.
*
*This program is distributed in the hope that it will be useful,
*but WITHOUT ANY WARRANTY; without even the implied warranty of
*MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*GNU Affero General Public License for more details.
*
*You should have received a copy of the GNU Affero General Public License
*along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*/

package nl.unimaas.msb.psbap;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

/**
 * Base class to read the configuration from a config file 
 * and load it into a Properties object
 * 
 * The class assumes a file named "config.properties" next to 
 * the application JAR and tries to load it
 * 
 * @author Ammar Ammar
 *
 */
public class Config {
	
	private static Properties properties = new Properties();
	private static InputStream input = null;

	static {

		try {

			input = new FileInputStream("config.properties");

			properties.load(input);

		} catch (IOException ex) {

			ex.printStackTrace();

		} finally {

			if (input != null) {
				try {
					input.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * Get the property value of the requested key
	 * 
	 * @param key is a string which you want to retrieve the value for
	 * @return a string holding the value of the requested key
	 */
	public static String getProperty(String key) {
		return properties.getProperty(key);
	}

}
