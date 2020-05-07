/**
* binding Pocket's SNPs effect on Binding Affinity Project (PSBAP) 
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

package nl.unimaas.msb.psbap.utils;

import java.util.ArrayList;
import java.util.List;

import com.univocity.parsers.csv.CsvParser;
import com.univocity.parsers.csv.CsvParserSettings;

/**
 * A utility class to use 48 amino acid properties in creating features for AA residues and their surrounding
 * It can calculate the AA properties for all residues surrounding a specific residue with 8A radius
 * and return the results as a list of double values
 * 
 * @author Ammar Ammar
 * 
 */
public class AAprops {
		
	
	/**
	 * create a header from the AA properties names in the resource CSV file
	 * @return the header as a list of strings
	 */
	public static List<String> getAApropsHeader() {

		List<String> header = new ArrayList<String>();

		CsvParserSettings settings = new CsvParserSettings();

		settings.getFormat().setLineSeparator("\n");
		settings.getFormat().setDelimiter(',');

		settings.setNumberOfRowsToSkip(1);

		CsvParser parser = new CsvParser(settings);

		List<String[]> rows = parser.parseAll(AAprops.class.getClass().getClassLoader().getResourceAsStream("/AAprops.csv"));
		

		for (String[] row : rows) {

			header.add(row[0]);
		}

		return header;
	}

}
