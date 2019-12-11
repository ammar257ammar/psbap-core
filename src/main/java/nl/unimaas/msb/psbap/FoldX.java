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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.apache.commons.io.FileUtils;

import nl.unimaas.msb.psbap.utils.DataHandler;

/**
 * A class that prepare mutation files for FoldX and generate reports after introducing mutation
 * to PDB structures
 * 
 * @author Ammar Ammar
 *
 */
public class FoldX {
	
	/**
	 * A method to generate folder structure and mutation files for FoldX
	 * @param pdbbindPocketVariants a list of pocket variants data
	 * @param inputPath the PdbBind entries folder path
	 * @param outputPath the FoldX processing path
	 * @throws IOException in case of error in IO operations
	 */
	public static void createMutationConfigFiles(List<String[]> pdbbindPocketVariants, String inputPath, String outputPath) throws IOException {
		
		Map<String, List<String[]>> mutationMap = new HashMap<String, List<String[]>>();
		
		for(String[] line: pdbbindPocketVariants) {
			
			String pdb = line[4];
			
			if(mutationMap.containsKey(pdb)) {
				mutationMap.get(pdb).add(new String[] {line[7]+line[10]+line[12]+line[8]+";"});
			}else {
				mutationMap.put(pdb, new ArrayList<String[]>());
				mutationMap.get(pdb).add(new String[] {line[7]+line[10]+line[12]+line[8]+";"});
			}
		}
		
		int count = 1;
		
		for (Entry<String, List<String[]>> entry : mutationMap.entrySet()){	
			
			File molFolder = new File(inputPath+entry.getKey());
			
			if(molFolder.isDirectory()) {
				
				File dir = new File(outputPath+molFolder.getName());
				
				boolean dirCreated = dir.mkdir();
				
				if(dirCreated) {

					File src = new File(inputPath+molFolder.getName()+"/"+molFolder.getName()+"_protein.pdb");
					File dst = new File(outputPath+molFolder.getName()+"/"+molFolder.getName()+"_protein.pdb");

					if(src.exists()) {
						
						FileUtils.copyFile(src, dst);	
						new File(outputPath+molFolder.getName()+"/input").mkdir();
						new File(outputPath+molFolder.getName()+"/output").mkdir();

						System.out.println(molFolder.getName()+ " copied");
					}else {
						System.out.println(molFolder.getName()+ " not copied");						
					}
					
				}
				System.out.println(entry.getKey() + " mutation list is ready " + count++);
				DataHandler.writeDatasetToTSV(entry.getValue(), outputPath + entry.getKey() + "/input/individual_list.txt");
			}
			
		}
	}
	

	/**
	 * A method to check FoldX results for the twl commands (RepairPDB and BuildModel) and report
	 * successes and failures
	 * 
	 * @param entriesPath of the FoldX processing folder
	 * @return a list of string arrays holding the results of FoldX
	 */
	public static List<String[]> getFoldxResults(String entriesPath) {
		
		File casf = new File(entriesPath);
		File[] mols = casf.listFiles();
		
		List<String[]> logResults = new ArrayList<String[]>();
		
		System.out.println(mols.length+ " files");
				
		String doneRepair = "failure";
		String doneResult = "failure";
		
		for(File molFolder: mols) {
						
			if(molFolder.isDirectory()) {
			
				doneRepair = "failure";
				doneResult = "failure";
								
				String repairPdbPath = entriesPath+"/"+molFolder.getName()+"/input"+"/log.txt";
				String buildModelPath = entriesPath+"/"+molFolder.getName()+"/output"+"/log.txt";					
				
				BufferedReader reader;
				try {
					reader = new BufferedReader(new FileReader(repairPdbPath));
					String line ;
					while ((line = reader.readLine()) != null) {
						if(line.equals("Cleaning RepairPDB...DONE")) {
							doneRepair = "success";
						}
					}
					reader.close();
					
					reader = new BufferedReader(new FileReader(buildModelPath));
					line = "";
					while ((line = reader.readLine()) != null) {
						if(line.equals("Cleaning BuildModel...DONE")) {
							doneResult = "success";
						}
					}
					reader.close();
					
					
					logResults.add(new String[] { molFolder.getName(), doneRepair, doneResult });						
					
				} catch (IOException e) {
					e.printStackTrace();
				}
				
			}
			
		}
		return logResults;
	}

}
