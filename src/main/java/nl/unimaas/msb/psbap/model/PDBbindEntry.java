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

package nl.unimaas.msb.psbap.model;

import java.util.List;

import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.io.sifts.SiftsEntity;
import org.biojava.nbio.structure.io.sifts.SiftsResidue;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 * A class that represent a PdbBind entry folder including all information (amino acids, ligands and SIFTS)
 * 
 * @author Ammar Ammar
 *
 */
public class PDBbindEntry {
	
	private Structure proteinStructure = null;
	private Structure pocketStructure = null;
	private IAtomContainer ligandStructure = null;
	
	private List<SiftsEntity> siftEntities = null;
	
	private List<SiftsResidue> siftResidues = null;
	
	private List<AminoAcid> pocketAminoAcids = null;

	private List<AminoAcid> proteinAminoAcids = null;

	private PDBFileReader reader;	
	
	private String pdb;

}
