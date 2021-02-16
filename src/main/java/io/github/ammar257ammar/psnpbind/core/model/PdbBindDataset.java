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

package io.github.ammar257ammar.psnpbind.core.model;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import io.github.ammar257ammar.psnpbind.core.Config;

/**
 * A class represents a PdbBind dataset object with methods to manipulate the dataset (filtering,
 * sorting, grouping and generating URLs)
 * 
 * @author Ammar Ammar
 *
 */
public class PdbBindDataset {
	
	private String pathGeneralFile;
	private String pathGeneralNamesFile;
	private String entriesPath;
	
	private List<String[]> pdbbindData = new ArrayList<String[]>();
	
	public enum PdbbindAttribute 
	{ 
	    PDB, RESOLUTION, LIGAND, UNIPROT, POCKET_RES_COUNT
	} 
	
	/**
	 * No-argument constructor initializes instance variables
	 */
	private PdbBindDataset() {
		this.pathGeneralFile = Config.getProperty("PDBBIND_DATA_PATH_1");
		this.pathGeneralNamesFile = Config.getProperty("PDBBIND_DATA_PATH_2");
		this.entriesPath = Config.getProperty("PDBBIND_ENTRIES_PATH");
	}

	/**
	 * PdbBindDataset Constructor
	 * 
	 * @param pathGeneralFile from the PdbBind downloaded file provides (PDB, Ligand, resolution) information 
	 * @param pathGeneralNamesFile from the PdbBind downloaded file provides Uniprot IDs
	 * @param entriesPath path of the PdbBind entries folders as extracted from the dataset downloaded file
	 */
	private PdbBindDataset(String pathGeneralFile, String pathGeneralNamesFile, String entriesPath) {
	    this.pathGeneralFile = pathGeneralFile;
	    this.pathGeneralNamesFile = pathGeneralNamesFile;
	    this.entriesPath = entriesPath;
	}
	
	/**
	 * Create a new PdbBindDataset instance
	 * 
	 * @return a new PdbBindDataset with initialized variables
	 */
	public static PdbBindDataset create() {
	    return new PdbBindDataset();
	}
	
	/**
	 * Create a new PdbBindDataset instance with arguments provided
	 * 
	 * @param pathGeneralFile from the PdbBind downloaded file provides (PDB, Ligand, resolution) information 
	 * @param pathGeneralNamesFile from the PdbBind downloaded file provides Uniprot IDs
	 * @param entriesPath path of the PdbBind entries folders as extracted from the dataset downloaded file
	 * @return a new PdbBindDataset
	 */
	public static PdbBindDataset create(String pathGeneralFile, String pathGeneralNamesFile, String entriesPath) {
	    return new PdbBindDataset(pathGeneralFile, pathGeneralNamesFile, entriesPath);
	}
	
	/**
	 * Load the data from the two files provided by PdbBind downloaded file (http://www.pdbbind.org.cn)
	 * @return the current PdbBindDataset after loading "pdbbindData" List with data
	 * from both files of PdbBind
	 */
	public PdbBindDataset loadData() {
		
		List<String[]> generalFileList = new ArrayList<String[]>(); // to hold the data from first file
		List<String[]> generalNamesFileList = new ArrayList<String[]>(); // to hold the data from second file
		
		try (Stream<String> stream = Files.lines(Paths.get(this.pathGeneralFile))) {
			
			generalFileList = stream.skip(6). // the first six lines are header, so skipped
					map(line -> { // map each line to a function that returns a String array with 3 values
						
						String[] lineArr = new String[3];
						
						// ligand ID is wrapped inside parenthesis so we need to extract it
						// by using substring between the index of the open and close parenthesis 
						int ligandStart = line.indexOf("(")+1;
						int ligandEnd = line.indexOf(")");
						
						lineArr[0] = line.substring(0,4); 					// PDB
						lineArr[1] = line.substring(6,10); 					// Resolution
						lineArr[2] = line.substring(ligandStart,ligandEnd); // Ligand
						
						return  lineArr;
					}).
					collect(Collectors.toList()); // return the results as List<String>
	
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		try (Stream<String> stream = Files.lines(Paths.get(this.pathGeneralNamesFile))) {
			
			generalNamesFileList = stream.skip(6). // the first six lines are header, so skipped
					map(line -> { // map each line to a function that returns a String array with 2 values
						
						String[] lineArr = new String[2];
						
						lineArr[0] = line.substring(0,4);   // PDB
						lineArr[1] = line.substring(12,18); // Uniprot
						
						return  lineArr;
					}).
					collect(Collectors.toList()); // return the results as List<String>
	
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		// Left join the Uniprot IDs list (2nd list) to the first one on PDB ID 
		// to merge the two lists in one
		for(String[] general : generalFileList) {
			for(String[] generalNames : generalNamesFileList) {
				
				// if the PDB ID matches, add a new row to the main List containing data from both lists 
				if(general[0].equals(generalNames[0])) {
					
					this.pdbbindData.add(new String[]{general[0], general[1], general[2], generalNames[1]});
					break;
				}
			}
		}
			
		return this;
	}
		
	/**
	 * Get the data of the PdbBindDataset object as List dataset
	 * @return the dataset list of the current PdbBindDataset object
	 */
	public List<String[]> getData(){
		
		return this.pdbbindData;		
	}	

	/**
	 * Filter a dataset to execlude instances that has an item that matches a provided string
	 * @param attr which is the column to be filtered in a PdbBindDataset object
	 * @param filterValue string value to filter the column against
	 * @return the PdbBindDataset object after applying the filter
	 */
	public PdbBindDataset filterStringNotEqual(PdbbindAttribute attr, String filterValue) {
		
		this.pdbbindData = this.pdbbindData.stream().
					filter(line -> !filterValue.equals(line[attr.ordinal()].trim())).
					collect(Collectors.toList());

		return this;
	}
	

	/**
	 *  Filter a dataset to execlude instances that has an item with a value lower than a provided double 
	 * @param attr which is the column to be filtered in a PdbBindDataset object
	 * @param cutoff a double value represent a threshold to filter values samller than it
	 * @return the PdbBindDataset object after applying the filter
	 */
	public PdbBindDataset filterDoubleCutoff(PdbbindAttribute attr, Double cutoff) {
		
		this.pdbbindData = this.pdbbindData.stream().
					filter(line -> Double.parseDouble(line[attr.ordinal()]) < cutoff).
					collect(Collectors.toList());
		return this;
	}
	
	/**
	 * Sort a PdbBindDataset object by a specified column
	 * @param attr the column to be sorted by
	 * @return the PdbBindDataset object after sorting
	 */
	public PdbBindDataset sortBy(PdbbindAttribute attr) {
		
		this.pdbbindData.sort(Comparator.comparing(line -> line[attr.ordinal()]));
		
		return this;
	}

	/**
	 * A method to keep one PDB for each Uniprot ID by grouping by Uniport ID and select the PDB
	 * with the minimum resolution for that Uniprot ID
	 * @param preserveLigandData boolean value to specify if ligands IDs 
	 * for all PDBs for a Uniprot ID should be saved
	 * @return the PdbBindDataset object after grouping
	 */
	public PdbBindDataset groupByUniProtAndKeepMinResolution(boolean preserveLigandData) {
		
		Map<Object, List<String[]>> map = this.pdbbindData.stream().
				collect(Collectors.groupingBy(line -> line[PdbbindAttribute.UNIPROT.ordinal()]));
		
		this.pdbbindData.clear();
				
		Iterator<?> it = map.entrySet().iterator();
		
	    while (it.hasNext()) {
	    	
	        @SuppressWarnings("unchecked")
			Map.Entry<String, List<String[]>> pair = (Map.Entry<String, List<String[]>>)it.next();
	       
	        List<String[]> ligandAndRes = (List<String[]>) pair.getValue();
	        List<String[]> ligandAndResTemp = ligandAndRes.stream().filter(line -> (!"3e5a".equals(line[0].trim()) && !"1z95".equals(line[0].trim()))).collect(Collectors.toList());     
	        String[] ligandAndResMin = Collections.min(ligandAndResTemp, Comparator.comparing(c -> c[1]));

	        if(preserveLigandData) {
	        	String ligandData = "";
	        	for(String[] row: ligandAndRes) {
	        		ligandData += row[PdbbindAttribute.PDB.ordinal()]+":"+row[PdbbindAttribute.LIGAND.ordinal()]+";";
	        	}
	        	ligandAndResMin[PdbbindAttribute.LIGAND.ordinal()] = ligandData;
	        }
	        
			this.pdbbindData.add(ligandAndResMin);

	        it.remove(); // avoids a ConcurrentModificationException
	    }
		
		return this;
	}
	

	/**
	 * A method to filter the PdbBindDataset to a folder structure like (coreset, refinedset, generalset)
	 * from the PdbBind dataset or any custom subset.
	 * @return the PdbBindDataset object after grouping
	 */
	public PdbBindDataset keepAsFolderMatch(){
		
		List<String[]> tempList = new ArrayList<String[]>();

		File casf = new File(this.entriesPath);
		File[] mols = casf.listFiles();
				
		for(File molFolder: mols) {
			if(molFolder.isDirectory()) {
					
				for(int i = 0; i < this.pdbbindData.size(); i++) {
					
					String[] row = this.pdbbindData.get(i);

					if(molFolder.getName().trim().toLowerCase().equals(row[0].trim().toLowerCase())) {
						tempList.add(row);
					}
				}
				
			}
		}
		
		this.pdbbindData = new ArrayList<String[]>(tempList);
		
		return this;
	}
	

	/**
	 * A method to create a dataset of SIFTS download URLs for PdbBindDataset PDB IDs
	 * @return  a list of SIFTS download URLs
	 */
	public List<String[]> asSiftsDownloadUrlsList(){
		
		List<String[]> newPdbbindData = new ArrayList<String[]>();

		for(int i = 0; i < this.pdbbindData.size(); i++) {

			String url = "http://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/" + 
						this.pdbbindData.get(i)[PdbbindAttribute.PDB.ordinal()] + 
						".xml.gz";
			
			newPdbbindData.add(new String[] {url});
			
		}
				
		return newPdbbindData;
	}
	

	/**
	 * A method to create a dataset of FASTA download URLs for PdbBindDataset UniProt IDs
	 * @return  a list of FASTA download URLs
	 */
	public List<String[]> asFastaDownloadUrlsList(){
		
		List<String[]> newPdbbindData = new ArrayList<String[]>();

		for(int i = 0; i < this.pdbbindData.size(); i++) {

			String url = "http://www.uniprot.org/uniprot/" + 
						this.pdbbindData.get(i)[PdbbindAttribute.UNIPROT.ordinal()] + 
						".fasta";
			
			newPdbbindData.add(new String[] {url});
			
		}
				
		return newPdbbindData;
	}
	

	/**
	 * A method to create a dataset of DSSP download URLs for PdbBindDataset PDB IDs
	 * @return  a list of DSSP download URLs
	 */
	public List<String[]> asDsspDownloadUrlsList(){
		
		List<String[]> newPdbbindData = new ArrayList<String[]>();

		for(int i = 0; i < this.pdbbindData.size(); i++) {

			String url = "http://files.rcsb.org/dssp/" +
					this.pdbbindData.get(i)[PdbbindAttribute.PDB.ordinal()].toLowerCase().substring(1, 3) + "/" +
					this.pdbbindData.get(i)[PdbbindAttribute.PDB.ordinal()].toLowerCase() + "/" +
					this.pdbbindData.get(i)[PdbbindAttribute.PDB.ordinal()].toLowerCase() + ".dssp.gz";
			
			newPdbbindData.add(new String[] {url});
			
		}
				
		return newPdbbindData;
	}
	
}
