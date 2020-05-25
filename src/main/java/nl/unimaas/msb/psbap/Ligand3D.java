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
import java.util.List;

import com.google.common.io.Files;

import nl.unimaas.msb.psbap.model.PdbBindDataset;
import nl.unimaas.msb.psbap.model.PdbBindDataset.PdbbindAttribute;


/**
 * A class to prepare ligands folder structure and prepare a dataset of ligands ChEMBL IDs and Tanimoto similarities
 * 
 * @author Ammar Ammar
 *
 */
public class Ligand3D {
	
	
	/**
	 * A method to generate folder structure and ligands files for OpenBabel
	 * @param foldxPath the path for FoldX folder where selected PDBbind entries reside
	 * @param pdbEntriesPath the PdbBind entries folder path to get the original ligands files
	 * @param output the ligands output folder to be used by OpenBabel
	 * @throws IOException in case of error in IO operations
	 */
	public static void prepareLigandsFolder(String foldxPath, String pdbEntriesPath, String output)
			throws IOException {

		PdbBindDataset pdbbindData = PdbBindDataset.create().loadData()
				.filterStringNotEqual(PdbbindAttribute.UNIPROT, "------")
				.filterStringNotEqual(PdbbindAttribute.RESOLUTION, "NMR").sortBy(PdbbindAttribute.RESOLUTION)
				.filterDoubleCutoff(PdbbindAttribute.RESOLUTION, 2.51).keepAsFolderMatch()
				.groupByUniProtAndKeepMinResolution(true);

		File casf = new File(foldxPath);
		File[] mols = casf.listFiles();

		List<String[]> pdbbindDataset = pdbbindData.getData();

		for (File molFolder : mols) {
			if (molFolder.isDirectory()) {

				for (String[] row : pdbbindDataset) {
					if (row[0].equals(molFolder.getName())) {

						String[] ligands = row[2].split(";");

						for (String ligand : ligands) {

							String[] ligandArr = ligand.split(":");

							new File(output + "/" + row[0]).mkdir();

							File ligandSrc = new File(
									pdbEntriesPath + "/" + ligandArr[0] + "/" + ligandArr[0] + "_ligand.sdf");
							File ligandDest = new File(output + "/" + row[0] + "/" + ligandArr[0] + "_ligand.sdf");

							Files.copy(ligandSrc, ligandDest);

							System.out.println("File copied: " + ligandArr[0]);

						}
					}
				}
			}
		}
	}

}
