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

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.List;

import com.univocity.parsers.tsv.TsvWriter;
import com.univocity.parsers.tsv.TsvWriterSettings;

/**
 * 
 * @author Ammar Ammar
 * 
 * A utility class to perform operations on datasets of type List<String[]>
 * It can print first rows of the dataset, return number of rows and columns, 
 * and write the dataset to TSV file
 */
public class DataHandler {
	
	

	/**
	 * write a dataset to TSV file
	 * @param dataset to be written to the filesystem
	 * @param path of the file to write the dataset to it
	 * @param header row for the dataset
	 */
	public static void writeDatasetToTSV(List<String[]> dataset, String path, String[] header){
		
		Writer outputWriter;
		
		try {
			
			outputWriter = new FileWriter(new File(path));

			TsvWriter writer = new TsvWriter(outputWriter, new TsvWriterSettings());
		
			if(header != null){
				writer.writeHeaders(header);
			}
			
			writer.writeStringRowsAndClose(dataset);
		
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	/**
	 * Overload writeDatasetToTSV(dataset, path, header) to write datasets without header
	 * @param dataset to be written to the filesystem
	 * @param path of the file to write the dataset to it
	 */
	public static void writeDatasetToTSV(List<String[]> dataset, String path){
		writeDatasetToTSV(dataset, path, null);
	}
	
	

}
