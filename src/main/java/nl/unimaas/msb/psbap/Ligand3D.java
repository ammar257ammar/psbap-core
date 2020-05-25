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
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;

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
	
	/**
	 * A method to extract ChEMBL IDs from the similar ligands SDF files selected with OpenBabel
	 * @param ligandsPath the OpenBabel-selected ligands folder path
	 * @throws IOException in case of error in IO operations
	 * @throws FileNotFoundException in case file not found
	 * @return a list of string arrays holding the ligands file names and IDs
	 */
	public static List<String[]> getLigandIDsFromFiles(String ligandsPath) throws FileNotFoundException, IOException {

		List<String[]> ligandsDataset = new ArrayList<String[]>();

		File casf = new File(ligandsPath);
		File[] mols = casf.listFiles();

		for (File molFolder : mols) {
			if (molFolder.isDirectory()) {

				File splitted = new File(molFolder.getAbsolutePath() + "/results");

				if (splitted.exists()) {

					File[] similarLigands = splitted.listFiles();

					for (File similarLigand : similarLigands) {

						int index = 1;

						try (IteratingSDFReader reader = new IteratingSDFReader(new FileInputStream(similarLigand),
								DefaultChemObjectBuilder.getInstance())) {

							while (reader.hasNext()) {

								IAtomContainer ac = (IAtomContainer) reader.next();

								for (Map.Entry<Object, Object> entry : ac.getProperties().entrySet()) {
									if (entry.getKey().equals("chembl_id")) {
										ligandsDataset.add(new String[] { molFolder.getName(),
												similarLigand.getName().substring(0,
														similarLigand.getName().lastIndexOf(".")) + "_" + index++,
												entry.getValue().toString() });
									}
								}

							}
						}
					}
				}
			}
		}
		return ligandsDataset;
	}
	
	/**
	 * A method to remove duplicated ligands from a selected ChEMBL IDs list
	 * @param similarLigands the OpenBabel-selected ligands list
	 * @return a list of string arrays holding the filtered ligands names and IDs 
	 */
	public static List<String[]> getLigandsIDsFiltered(List<String[]> similarLigands) {

		Map<Object, List<String[]>> map = similarLigands.stream().collect(Collectors.groupingBy(line -> line[0]));

		similarLigands.clear();

		Iterator<?> it = map.entrySet().iterator();

		while (it.hasNext()) {

			@SuppressWarnings("unchecked")
			Map.Entry<String, List<String[]>> pair = (Map.Entry<String, List<String[]>>) it.next();

			List<String[]> ligandAndRes = (List<String[]>) pair.getValue();

			Map<Object, List<String[]>> map2 = ligandAndRes.stream().collect(Collectors.groupingBy(line -> line[2]));

			Iterator<?> it2 = map2.entrySet().iterator();

			while (it2.hasNext()) {

				@SuppressWarnings("unchecked")
				Map.Entry<String, List<String[]>> pair2 = (Map.Entry<String, List<String[]>>) it2.next();

				List<String[]> ligandAndRes2 = (List<String[]>) pair2.getValue();

				String[] ligandAndResMin = ligandAndRes2.get(0);
				ligandAndRes2.remove(0);

				for (String[] row : ligandAndRes2) {
					new File(Config.getProperty("LIGANDS_PATH") + "/" + row[0] + "/splitted/" + row[1]).delete();
				}

				similarLigands.add(ligandAndResMin);

				it2.remove(); // avoids a ConcurrentModificationException
			}

			it.remove();
		}

		return similarLigands;
	}
	
	/**
	 * A wrapper method to extract ChEMBL IDs from the similar ligands and remove duplicates
	 * @param ligandsPath the OpenBabel-selected ligands folder path
	 * @throws IOException in case of error in IO operations
	 * @throws FileNotFoundException in case file not found
	 * @return a list of string arrays holding the filtered ligands names and IDs 
	 */
	public static List<String[]> getLigandsIDsFiltered(String ligandsPath) throws FileNotFoundException, IOException {

		List<String[]> similarLigands = Ligand3D.getLigandIDsFromFiles(ligandsPath);

		similarLigands = Ligand3D.getLigandsIDsFiltered(similarLigands);

		return similarLigands;
	}

}
