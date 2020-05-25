
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

package nl.unimaas.msb.psbap;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.LineIterator;

import nl.unimaas.msb.psbap.Config;
import nl.unimaas.msb.psbap.model.PdbBindDataset.PdbbindAttribute;


/**
 * This class map between the Uniprot variants dataset with PdbBind dataset and filter the 
 * missence variants that occur in PdbBind proteins only
 * 
 * @author Ammar Ammar
 * 
 */
public class UniProtVariantsMapper {
	
	/**
	 * A method to filter missense variants that occur in PdbBind proteins
	 * @param pdbbindDataset to filter UniProt variants against
	 * @return a dataset with all missence variants in PdbBind proteins
	 */
	public static List<String[]> mapMissenseVariantsToPdbbindDataset(List<String[]> pdbbindDataset) {
		
		List<String[]> pdbbindVariantsRows = new ArrayList<String[]>();
		
		try {
			LineIterator it = FileUtils.lineIterator(new File(Config.getProperty("UNIPROT_VARIANTS_PATH")), "UTF-8");
		
			try {
				
			    while (it.hasNext()) {
			    				    	
			    	String line = it.nextLine();

			    	String[] lineArr = line.split("\t");

			        if(lineArr.length > 5 && lineArr[4].toLowerCase().contains("missense")) {

			        	for(String[] row: pdbbindDataset) {

			        		if(lineArr[1].trim().toLowerCase().equals(row[PdbbindAttribute.UNIPROT.ordinal()].trim().toLowerCase())) {
			        			
			        			pdbbindVariantsRows.add(new String[] {lineArr[1], 
					        			lineArr[2], 
					        			lineArr[3],
					        			lineArr[4],
					        			row[PdbbindAttribute.PDB.ordinal()],
					        			row[PdbbindAttribute.RESOLUTION.ordinal()],
					        			row[PdbbindAttribute.LIGAND.ordinal()]
					        			});
			        		}
			        	}
			        }	    
			    }
				
			} finally {
				it.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return pdbbindVariantsRows;
	}
	
}
