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

}
