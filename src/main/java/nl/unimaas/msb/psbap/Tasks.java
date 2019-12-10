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

import nl.unimaas.msb.psbap.model.PdbBindDataset;
import nl.unimaas.msb.psbap.model.PdbBindDataset.PdbbindAttribute;
import nl.unimaas.msb.psbap.utils.DataHandler;

/**
 * A class with static methods to perform the steps in the PSBAP project via calling them from command line
 * 
 * @author Ammar Ammar
 *
 */
public class Tasks {
	
	/**
	 * This is usually the first step in the workflow. 
	 * It generates PdbBind dataset filtered, sorted and grouped by UniProt ID.
	 * It uses a threshold of 2.51 Angstrom and remove PDBs obtained with NMR and PDBs withou UniProt ID
	 * 
	 * @return PdbBindDataset which is the main block for the project
	 */
	public static PdbBindDataset generatePdbBindDataset(){
		
		PdbBindDataset pdbbindData = PdbBindDataset.create().
											loadData().
											filterStringNotEqual(PdbbindAttribute.UNIPROT, "------").
											filterStringNotEqual(PdbbindAttribute.RESOLUTION, "NMR").
											sortBy(PdbbindAttribute.RESOLUTION).
											keepAsFolderMatch().
											filterDoubleCutoff(PdbbindAttribute.RESOLUTION, 2.51).
											groupByUniProtAndKeepMinResolution(true);
		
		return pdbbindData;		
	}
	
	/**
	 * A method to write PdbBindDataset list to a TSV file
	 * @param pdbbindData which is the generated PdbBindDataset object
	 * @param path which is the path of TSV file
	 */
	public static void writePdbBindDatasetToFile(PdbBindDataset pdbbindData, String path){
		
    	DataHandler.writeDatasetToTSV(pdbbindData.getData(), path + "/pdbbind_entries_data.tsv");
	}

}
