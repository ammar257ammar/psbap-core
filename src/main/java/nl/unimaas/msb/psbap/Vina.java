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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.vecmath.Point3d;

import org.apache.commons.io.FileUtils;
import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.ChainImpl;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureImpl;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.PDBFileReader;

import nl.unimaas.msb.PockerSnpBindingAffinity.Vina;
import nl.unimaas.msb.psbap.model.PDBbindEntry;
import nl.unimaas.msb.psbap.utils.PdbTools;

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

	/**
	 * A method to generate a BioJava Structure from the binding pocket of PDBbind entry protein after minimization
	 * @param pdb the PdbBind entry protein
	 * @param finalPDBPath the path of the minimzed PDB (wild-type or mutation)
	 * @throws IOException in case of error in IO operations
	 */
	public static Structure getPocketForModifiedPDBbindStructure(String pdb, String finalPDBPath) throws IOException{
		
		PDBbindEntry pdbEntry = new PDBbindEntry(pdb, false, false);

		PDBFileReader reader = PdbTools.configureReader(false);

		Structure proteinStructure = reader.getStructure(finalPDBPath);
		
		List<AminoAcid> proteinAAlist = PdbTools.getAminoAcidsFromStructure(proteinStructure);

		Structure s = new StructureImpl();

		Chain c = new ChainImpl();
    	
    	for(AminoAcid aa : pdbEntry.getPocketAminoAcids()){
    		
    		for (AminoAcid aap : proteinAAlist) {
    			
    			if(aa.getResidueNumber().toString().equals(aap.getResidueNumber().toString())) {

    				c.addGroup(aap);

    				break;
    			}		
    		}
    	}
    	
    	s.addChain(c);
    	
		return s;
	}
	
	/**
	 * A method to generate Vina grid box config from the binding pocket of PDBbind entry protein after minimization
	 * @param pdb the PdbBind entry protein
	 * @param finalPDBPath the path of the minimzed PDB (wild-type or mutation)
	 * @throws IOException in case of error in IO operations
	 */
	public static List<double[]> calculateVinaGridEnhanced(String pdb, String finalPDBPath) throws IOException{

		Structure pocketStructure = Vina.getPocketForModifiedPDBbindStructure(pdb, finalPDBPath);

    	double minX = 1000.0;
    	double minY = 1000.0;
    	double minZ = 1000.0;
    	
    	double maxX = 0.0;
    	double maxY = 0.0;
    	double maxZ = 0.0;
		
    	Atom[] pocketAtoms = StructureTools.getAllAtomArray(pocketStructure);
		
		for(Atom atom: pocketAtoms){
			
			Point3d p = atom.getCoordsAsPoint3d();

			if(p.x < minX)	minX = p.x;
			
			if(p.y < minY)  minY = p.y;
			
			if(p.z < minZ)	minZ = p.z;
			
			if(p.x > maxX)	maxX = p.x;

			if(p.y > maxY)	maxY = p.y;

			if(p.z > maxZ)	maxZ = p.z;
		}
		
		List<double[]> grid = new ArrayList<double[]>();
    	
    	grid.add(new double[]{minX , minY , minZ});
    	grid.add(new double[]{maxX , maxY , maxZ});
    	
    	grid.add(new double[]{(Math.round((((maxX + minX)/2)*1000.0)))/1000.0,
    						  (Math.round((((maxY + minY)/2)*1000.0)))/1000.0,
    						  (Math.round((((maxZ + minZ)/2)*1000.0)))/1000.0});
    	
    	grid.add(new double[]{(int)(Math.ceil(maxX - minX)),
			    			  (int)(Math.ceil(maxY - minY)),
			    			  (int)(Math.ceil(maxZ - minZ))});
    	
    	return grid;
	}
	
	

	/**
	 * A method to generate Vina docking config for PDBbind entries
	 * @param entriesPath the path of the selected PDBbind entries
	 * @param replace a boolean to choose if the config file should be replaced if exists
	 * @throws IOException in case of error in IO operations
	 */
	public static void createConfigFilesAfterLigandsStructure(String entriesPath, boolean replace) throws IOException {
		
		File entries = new File(entriesPath);
		File[] mols = entries.listFiles();
		
		System.out.println(mols.length+ " files");
		
		for(File molFolder: mols) {
			if(molFolder.isDirectory()) {
				
				File molVarsFolder = new File(entriesPath + molFolder.getName()+"/proteins");
				File[] molVars = molVarsFolder.listFiles();
				
				for(File molVarFolder: molVars) {
					
					if(molVarFolder.isDirectory()) {
						
						File configFile = new File(entriesPath+"/"+molFolder.getName()+"/proteins/" + molVarFolder.getName() +"/"+ molVarFolder.getName()+"_config.txt");
						
						if(replace) {
						
							List<double[]> grid = Vina.calculateVinaGridEnhanced(molFolder.getName(),entriesPath+"/"+molFolder.getName()+"/proteins/" + molVarFolder.getName() +"/"+ molVarFolder.getName()+"_final.pdb");
							
							String gridBox =    "center_x = "+grid.get(2)[0]+"\n" + 
								    			"center_y = "+grid.get(2)[1]+"\n" + 
								    			"center_z = "+grid.get(2)[2]+"\n" + 
								    			"\n" + 
								    			"size_x = "+(int) grid.get(3)[0]+"\n" + 
								    			"size_y = "+(int) grid.get(3)[1]+"\n" + 
								    			"size_z = "+(int) grid.get(3)[2]+"\n\n";
							
							String technical =  "cpu = 12\n" + 
												"num_modes = 3\n" + 
												"energy_range = 2\n" + 
												"exhaustiveness = 12\n";
							
							
							try (BufferedWriter writer = new BufferedWriter(new FileWriter(configFile))) 
							{
								writer.write(gridBox);
								writer.write(technical);
							}
							System.out.println(molVarFolder.getName()+ " config written!!");
							
						}else {
							
							if(!configFile.exists()) {
								
								List<double[]> grid = Vina.calculateVinaGridEnhanced(molFolder.getName(),entriesPath+"/"+molFolder.getName()+"/proteins/" + molVarFolder.getName() +"/"+ molVarFolder.getName()+"_final.pdb");
								
								String gridBox =    "center_x = "+grid.get(2)[0]+"\n" + 
									    			"center_y = "+grid.get(2)[1]+"\n" + 
									    			"center_z = "+grid.get(2)[2]+"\n" + 
									    			"\n" + 
									    			"size_x = "+(int) grid.get(3)[0]+"\n" + 
									    			"size_y = "+(int) grid.get(3)[1]+"\n" + 
									    			"size_z = "+(int) grid.get(3)[2]+"\n\n";
								
								String technical =  "cpu = 12\n" + 
													"num_modes = 3\n" + 
													"energy_range = 2\n" + 
													"exhaustiveness = 12\n";
								
								
								try (BufferedWriter writer = new BufferedWriter(new FileWriter(configFile))) 
								{
									writer.write(gridBox);
									writer.write(technical);
								}
								System.out.println(molVarFolder.getName()+ " config written!!");
							}
						}
						
					}
				}
			}	
		}
	}
	
	

}
