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

package nl.unimaas.msb.psbap.utils;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 * A utility class to manipulate ligand files mainly using CDK (Mol2 and SDF formats)
 * 
 * @author Ammar Ammar
 *
 */
public class LigandTools {


	/**
	 * A method to read a molecule file from PdbBind entry folder and parse it
	 * as a CDK IAtomContainer
	 * @param fname the molecule file path
	 * @return a CDK IAtomContainer
	 */
	public static IAtomContainer readMolandAromatizeKekulize(String fname){
		
		String filename = fname;
		
		try{
		
			InputStream ins = new FileInputStream(new File(filename));
	        
			MDLV2000Reader reader = new MDLV2000Reader(ins);
	        
	        ChemFile chemFile = (ChemFile) reader.read((ChemObject) new ChemFile());
	        
	        IChemSequence cs = chemFile.getChemSequence(0);
	                
	        IChemModel cm = cs.getChemModel(0);
	        
	        IAtomContainerSet acs = cm.getMoleculeSet();
	                
	        reader.close();
	        
	        IAtomContainer ac = acs.getAtomContainer(0);
	        
	        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);
	        
	        ElectronDonation ed = ElectronDonation.daylight();
			
			Aromaticity ar = new Aromaticity(ed, Cycles.all());
			
			ar.apply(ac);
				                        
	        return ac;  
	        
		} catch (CDKException e) {
			return null;
		} catch (IOException e) {
			return null;
		}
	}
}
