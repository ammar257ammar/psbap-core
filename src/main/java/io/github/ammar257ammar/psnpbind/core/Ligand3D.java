/**
* Binding Pocket SNPs' effect on Binding Affinity Database Project (PSnpBind)
* 
*Copyright (C) 2019-2021  Ammar Ammar <ammar257ammar@gmail.com> ORCID:0000-0002-8399-8990
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

package io.github.ammar257ammar.psnpbind.core;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.jooq.lambda.Seq;
import org.jooq.lambda.tuple.Tuple2;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;

import com.google.common.io.Files;
import com.univocity.parsers.tsv.TsvParser;
import com.univocity.parsers.tsv.TsvParserSettings;

import io.github.ammar257ammar.psnpbind.core.model.PdbBindDataset;
import io.github.ammar257ammar.psnpbind.core.model.PdbBindDataset.PdbbindAttribute;


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
										ligandsDataset.add(
												new String[] { 
													molFolder.getName(),
													similarLigand.getName().substring(0,similarLigand.getName().lastIndexOf(".")) + "_" + index++,
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
	public static List<String[]> getLigandsIDsFiltered(List<String[]> similarLigands, boolean keepall) {

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
					new File(Config.getProperty("LIGANDS_PATH") + "/" + row[0] + "/splitted/" + row[1] + ".mol2").delete();
					new File(Config.getProperty("LIGANDS_PATH") + "/" + row[0] + "/splitted-smi/" + row[1] + ".smi").delete();
					
					if(keepall) {
						similarLigands.add(row);						
					}
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
	public static List<String[]> getLigandsIDsFiltered(String ligandsPath, boolean keepall) throws FileNotFoundException, IOException {

		List<String[]> similarLigands = Ligand3D.getLigandIDsFromFiles(ligandsPath);

		similarLigands = Ligand3D.getLigandsIDsFiltered(similarLigands, keepall);

		return similarLigands;
	}

	
	/**
	 * A method to combine ligands information with Tanimoto similarity
	 * @param ligandsPath the OpenBabel-selected ligands folder path
	 * @param idsFile ligands IDs file
	 * @throws IOException in case of error in IO operations
	 * @throws FileNotFoundException in case file not found
	 * @return a list of string arrays holding the filtered ligands names, IDs and Tanimoto similarity 
	 */
	public static List<String[]> combineIDsAndTanimotoOfLigands(String ligandsPath, String idsFile, boolean keepall)
			throws FileNotFoundException, IOException {

		List<String[]> ligandsDataset = new ArrayList<String[]>();

		List<String[]> logsDataset = new ArrayList<String[]>();

		File casf = new File(ligandsPath);
		File[] mols = casf.listFiles();

		for (File molFolder : mols) {
			if (molFolder.isDirectory()) {

				File logsFolder = new File(molFolder.getAbsolutePath() + "/logs");

				if (logsFolder.exists()) {

					File[] logs = logsFolder.listFiles();

					for (File log : logs) {

						try (BufferedReader reader = new BufferedReader(new FileReader(log))) {
							String line;
							boolean moleculeDefined = false;
							String molecule = "";

							while ((line = reader.readLine()) != null) {

								if (line.startsWith(">")) {

									if (!moleculeDefined) {
										molecule = line.substring(1, 5);
										moleculeDefined = true;
									} else {

										logsDataset.add(new String[] { molFolder.getName(), molecule,
												line.substring(1, line.indexOf("Tanimoto")).trim(),
												line.substring(line.lastIndexOf("=") + 1, line.length()).trim() });
									}
								}
							}
						} // try

					} // for logs
				} // if log exists
			}
		} // for molecules in /pdb

		TsvParserSettings settings = new TsvParserSettings();
		settings.getFormat().setLineSeparator("\n");

		TsvParser parser = new TsvParser(settings);

		List<String[]> rows = parser.parseAll(new File(idsFile));

		Seq<String[]> seqLogs = Seq.seq(logsDataset);

		Seq<String[]> seqIds = Seq.seq(rows);

		List<Tuple2<String[], String[]>> joinDataset = seqIds
				.leftOuterJoin(seqLogs, (t, u) -> t[0].equals(u[0]) && t[2].equals(u[2])).collect(Collectors.toList());

		for (Tuple2<String[], String[]> tuple : joinDataset) {
			ligandsDataset.add(new String[] { tuple.v1[0], tuple.v1[1], tuple.v1[2], tuple.v2[3] });
		}

		ligandsDataset = Ligand3D.getLigandsIDsFiltered(ligandsDataset, keepall);

		return ligandsDataset;

	}
	
	/**
	 * A method to combine ligands information with Tanimoto similarity and SMILES
	 * @param ligandsPath the OpenBabel-selected ligands folder path
	 * @param idsFile ligands IDs file
	 * @throws IOException in case of error in IO operations
	 * @throws FileNotFoundException in case file not found
	 * @return a list of string arrays holding the filtered ligands names, IDs, Tanimoto similarity and SMILES 
	 */
	public static List<String[]> combineIDsAndTanimotoAndSmilesOfLigands(String ligandsPath, String idsFile, boolean keepall)
			throws FileNotFoundException, IOException {

		List<String[]> ligandsDataset = new ArrayList<String[]>();
		List<String[]> ligandsDatasetFinal = new ArrayList<String[]>();

		List<String[]> logsDataset = new ArrayList<String[]>();

		List<String[]> smilesDataset = new ArrayList<String[]>();

		
		File casf = new File(ligandsPath);
		File[] mols = casf.listFiles();

		for (File molFolder : mols) {
			if (molFolder.isDirectory()) {

				File logsFolder = new File(molFolder.getAbsolutePath() + "/logs");

				if (logsFolder.exists()) {

					File[] logs = logsFolder.listFiles();

					for (File log : logs) {

						try (BufferedReader reader = new BufferedReader(new FileReader(log))) {
							String line;
							boolean moleculeDefined = false;
							String molecule = "";

							while ((line = reader.readLine()) != null) {

								if (line.startsWith(">")) {

									if (!moleculeDefined) {
										molecule = line.substring(1, 5);
										moleculeDefined = true;
									} else {

										logsDataset.add(new String[] { molFolder.getName(), molecule,
												line.substring(1, line.indexOf("Tanimoto")).trim(),
												line.substring(line.lastIndexOf("=") + 1, line.length()).trim() });
									}
								}
							}
						} // try

					} // for logs
				} // if log exists
				
				
				File smiFolder = new File(molFolder.getAbsolutePath() + "/results-smi");

				if (smiFolder.exists()) {

					File[] smiles = smiFolder.listFiles();

					for (File smile : smiles) {

						try (BufferedReader reader = new BufferedReader(new FileReader(smile))) {
							String line;

							while ((line = reader.readLine()) != null) {

								if (line.split("\t").length == 2) {

									smilesDataset.add(new String[] { molFolder.getName(), smile.getName().substring(0,4),
												line.split("\t")[1].trim(),
												line.split("\t")[0].trim() });
								}
							}
						}
					} // for smiles
				} // if smile exists
			} // if (molFolder.isDirectory())
		} // for molecules in /pdb

		TsvParserSettings settings = new TsvParserSettings();
		settings.getFormat().setLineSeparator("\n");

		TsvParser parser = new TsvParser(settings);

		List<String[]> rows = parser.parseAll(new File(idsFile));

		Seq<String[]> seqLogs = Seq.seq(logsDataset);

		Seq<String[]> seqSmiles = Seq.seq(smilesDataset);

		Seq<String[]> seqIds = Seq.seq(rows);

		List<Tuple2<String[], String[]>> joinDataset = seqIds
				.leftOuterJoin(seqLogs, (t, u) -> t[0].equals(u[0]) && t[2].equals(u[2])).collect(Collectors.toList());

		for (Tuple2<String[], String[]> tuple : joinDataset) {
			ligandsDataset.add(new String[] { tuple.v1[0], tuple.v1[1], tuple.v1[2], tuple.v2[3] });
		}

		// A second join
		Seq<String[]> seqligands = Seq.seq(ligandsDataset);

		List<Tuple2<String[], String[]>> joinDataset2 = seqligands
				.leftOuterJoin(seqSmiles, (t, u) -> t[0].equals(u[0]) && t[2].equals(u[2])).collect(Collectors.toList());

		for (Tuple2<String[], String[]> tuple : joinDataset2) {
			ligandsDatasetFinal.add(new String[] { tuple.v1[0], tuple.v1[1], tuple.v1[2], tuple.v1[3], tuple.v2[3] });
		}
		
		
		ligandsDatasetFinal = Ligand3D.getLigandsIDsFiltered(ligandsDatasetFinal, keepall);

		return ligandsDatasetFinal;

	}

}
