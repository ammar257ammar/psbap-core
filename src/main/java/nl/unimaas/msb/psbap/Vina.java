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

import org.apache.commons.io.FileUtils;


/**
 * A class to prepare AutoDock Vina folder structure and prepare the grid box and extract binding affinities from docking results
 * 
 * @author Ammar Ammar
 *
 */
public class Vina {
	
	
	/**
	 * A method to generate folder structure for AutoDock Vina from Gromacs EM results folders for one PDB bind entry
	 * @param gromacsPath the path for Gromacs folder where EM took place
	 * @param pdb the PdbBind entry
	 * @param vinaPath the AutoDock Vina folder being prepared docking
	 * @throws IOException in case of error in IO operations
	 */
	public static void createFolderStructureFromGromacsOneFolder(String gromacsPath, String pdb, String vinaPath) throws IOException {

		int index = 1;

		File molVarsFolder = new File(gromacsPath+pdb+"/output");
		File[] molVars = molVarsFolder.listFiles();
		
		File newMolFolder = new File(vinaPath+pdb);
		
		if(!newMolFolder.exists()) {
			newMolFolder.mkdir();
		}
		
		for(File molVarFolder: molVars) {
			
			if(molVarFolder.isDirectory()) {
				
				File newMolVarFolder = new File(vinaPath+pdb+"/proteins/"+molVarFolder.getName());

				if(!newMolVarFolder.exists()) {
					newMolVarFolder.mkdir();
				}

				String oldMolVarFolder = gromacsPath+pdb+"/output/"+molVarFolder.getName();
				String newMolVarFolderStr = vinaPath+pdb+"/proteins/"+molVarFolder.getName();

				File src1 = new File(oldMolVarFolder+"/"+molVarFolder.getName()+"_final.pdb");
				File dst1 = new File(newMolVarFolderStr+"/"+molVarFolder.getName()+"_final.pdb");
				
				FileUtils.copyFile(src1, dst1);	
				
				System.out.println("Varinat folder num: "+index);
				index++;
				
			}
		}	
	}
	

	/**
	 * A method to generate folder structure for AutoDock Vina from Gromacs EM results folders
	 * @param gromacsPath the path for Gromacs folder where EM took place
	 * @param vinaPath the AutoDock Vina folder being prepared docking
	 * @throws IOException in case of error in IO operations
	 */
	public static void createFolderStructureFromGromacsFolder(String gromacsPath, String vinaPath) throws IOException {
		
		File entries = new File(gromacsPath);
		File[] mols = entries.listFiles();
		
		System.out.println(mols.length+ " files");
				
		for(File molFolder: mols) {
			if(molFolder.isDirectory()) {
				
				Vina.createFolderStructureFromGromacsOneFolder(gromacsPath, molFolder.getName(), vinaPath);
				
			}
		}
	}
	
	
	/**
	 * A method to generate ligands folder structure for AutoDock Vina for single PDBbind entry
	 * @param vinaPath the AutoDock Vina folder being prepared docking
	 * @param pdb the PdbBind entry
	 * @param ligandsPath the path for ligands prepared in a previous step
	 * @param minimized a boolean value to choose if the minimized ligands should be used or the original ones
	 * @throws IOException in case of error in IO operations
	 */
	public static void createLigandsStructureAfterVinaFolderSingle(String vinaPath, String pdb, String ligandsPath, boolean minimized) throws IOException {

		File ligandsDir = new File(vinaPath+pdb+"/ligands");
		
		if(!ligandsDir.exists()) {
			ligandsDir.mkdir();
		}
		
		File molVarsFolder = new File(vinaPath+pdb+"/proteins");
		File[] molVars = molVarsFolder.listFiles();
						
		for(File molVarFolder: molVars) {
			
			if(molVarFolder.isDirectory()) {
				
				String varVinaFolder = vinaPath+pdb+"/proteins/"+molVarFolder.getName();
				
				File dir = new File(varVinaFolder+"/vina");
				
				boolean dirCreated = dir.mkdir();
				
				if(dirCreated || dir.exists()) {
					
					File ligandsFoldersDir = new File(ligandsPath+pdb+"/splitted");
					File[] ligandsFiles = ligandsFoldersDir.listFiles();
					
					for(File ligandsFile: ligandsFiles) {
						
						if(minimized) {
							
							if(ligandsFile.getName().contains("_min")) {
								
								String ligandDockingFolderName = ligandsFile.getName().substring(0, ligandsFile.getName().length()-9); 
							
								File ligandDockingFolder = new File(varVinaFolder+"/vina/"+ligandDockingFolderName);
								
								boolean ligandDockingFolderCreated = ligandDockingFolder.mkdir();
								
								if(ligandDockingFolderCreated || ligandDockingFolder.exists()) {
								
									File srcLigand = new File(ligandsPath+pdb+"/splitted/"+ligandsFile.getName());
									File dstLigand = new File(vinaPath+pdb+"/ligands/"+ligandsFile.getName());
									
									if(srcLigand.exists()) {
										
										FileUtils.copyFile(srcLigand, dstLigand);
										
										System.out.println(ligandsFile.getName() + " Ligand copied");
									
									}else {
										System.out.println(srcLigand.getName()+ " ligand file not copied #####################################");						
									}
									
								}else {
									System.out.println(ligandDockingFolderName+ " Ligand folder cannot be created!!");						
								}
									
							} // if(ligandsFile.getName().contains("_min")) {
							
						}else {
							
							if(!ligandsFile.getName().contains("_min")) {
								
								String ligandDockingFolderName = ligandsFile.getName().substring(0, ligandsFile.getName().length()-5); 
							
								File ligandDockingFolder = new File(varVinaFolder+"/vina/"+ligandDockingFolderName);
								
								boolean ligandDockingFolderCreated = ligandDockingFolder.mkdir();
								
								if(ligandDockingFolderCreated || ligandDockingFolder.exists()) {
								
									File srcLigand = new File(ligandsPath+pdb+"/splitted/"+ligandsFile.getName());
									File dstLigand = new File(vinaPath+pdb+"/ligands/"+ligandDockingFolderName+"_min.mol2");
									
									if(srcLigand.exists()) {
										
										FileUtils.copyFile(srcLigand, dstLigand);
										
										System.out.println(ligandsFile.getName() + " Ligand copied");
									
									}else {
										System.out.println(srcLigand.getName()+ " ligand file not copied #####################################");						
									}
									
								}else {
									System.out.println(ligandDockingFolderName+ " Ligand folder cannot be created!!");						
								}
									
							} // if(ligandsFile.getName().contains("_min")) {
							
						}
					} // for(File ligandsFile: ligandsFiles) {			
				}
			}
		}
	}
	

	/**
	 * A method to generate ligands folder structure for AutoDock Vina for all PDBbind entry
	 * @param vinaPath the AutoDock Vina folder being prepared docking
	 * @param ligandsPath the path for ligands prepared in a previous step
	 * @param minimized a boolean value to choose if the minimized ligands should be used or the original ones
	 * @throws IOException in case of error in IO operations
	 */
	public static void createLigandsStructureAfterVinaFolder(String vinaPath, String ligandsPath, boolean minimized) throws IOException {
		
		File entries = new File(vinaPath);
		File[] mols = entries.listFiles();
		
		for(File molFolder: mols) {
			if(molFolder.isDirectory()) {
				
				Vina.createLigandsStructureAfterVinaFolderSingle(vinaPath, molFolder.getName(), ligandsPath, minimized);	
			}
		}
	}

}
