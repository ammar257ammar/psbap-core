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

import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.io.sifts.SiftsEntity;
import org.biojava.nbio.structure.io.sifts.SiftsResidue;
import org.openscience.cdk.interfaces.IAtomContainer;

import nl.unimaas.msb.psbap.Config;
import nl.unimaas.msb.psbap.utils.LigandTools;
import nl.unimaas.msb.psbap.utils.PdbTools;

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
	
	
	/**
	 * Non-argument constructor
	 */
	public PDBbindEntry() {
		this("", false, false);
	}

	/**
	 * A constructor for the PDBbindEntry object class that initializes the object with seveal
	 * parsed information (protein structure, pocket structure, ligand structure, SIFTS data)
	 * 
	 * @param pdb the PDB ID from the PdbBind dataset
	 * @param parseLigand boolean to specifiy if the secondary structure information and sequence 
	 * alignment should be performed by BioJava when parsing the PDB files
	 * @param alignSeqAndPraseSS a boolean to specify if SS and alignment should be performed
	 */
	public PDBbindEntry(String pdb, boolean parseLigand, boolean alignSeqAndPraseSS) {
		
		this.pdb = pdb;
		this.reader = PdbTools.configureReader(alignSeqAndPraseSS);
				
		try {
			this.proteinStructure = reader.getStructure(Config.getProperty("PDBBIND_ENTRIES_PATH")+"/"+pdb+"/"+pdb+"_protein.pdb");
			this.proteinAminoAcids = PdbTools.getAminoAcidsFromStructure(this.proteinStructure);		
		} catch (IOException e) {
			this.proteinStructure = null;
		}
		
		try {
			this.pocketStructure = reader.getStructure(Config.getProperty("PDBBIND_ENTRIES_PATH")+"/"+pdb+"/"+pdb+"_pocket.pdb");
			this.pocketAminoAcids = PdbTools.getAminoAcidsFromStructure(this.pocketStructure);
		} catch (IOException e) {
			this.pocketStructure = null;
		}
		 
		if(parseLigand){
			this.ligandStructure = LigandTools.readMol2andAddHydrogens(
					Config.getProperty("PDBBIND_ENTRIES_PATH") + "/" + pdb + "/" + pdb + "_ligand.mol2", true);
		}
		
		this.siftEntities = PdbTools.getSiftsEntitiesForPDB(Config.getProperty("SIFTS_PATH"), pdb);
		
		if(this.siftEntities != null){
			this.siftResidues = PdbTools.getSiftResiduesFromSiftsMapping(this.siftEntities);		
		}
		
	}
	
	/**
	 * check if the current entry has protein structure (i.e. the PDB files is parsed correctly)
	 * @return a boolean value True or False
	 */
	public boolean hasProteinStructure(){
		return this.proteinStructure != null;
	}

	/**
	 * check if the current entry has pocket structure (i.e. the PDB files is parsed correctly)
	 * @return a boolean value True or False
	 */
	public boolean hasPocketStructure(){
		return this.pocketStructure != null;
	}
	
	/**
	 * check if the current entry has ligand structure (i.e. the Mol files is parsed correctly)
	 * @return a boolean value True or False
	 */
	public boolean hasLigandStructure(){
		return this.ligandStructure != null;
	}
	
	/**
	 * check if the current entry has a SIFTS mapping (i.e. the SIFTS file is parsed correctly)
	 * @return a boolean value True or False
	 */
	public boolean hasSiftsMapping(){
		return this.siftEntities != null;
	}

	/**
	 * Get the protein structure as BioJava Structure object
	 * @return a BioJava Structure object
	 */
	public Structure getProteinStructure() {
		return proteinStructure;
	}

	/**
	 * Set the protein structure as BioJava Structure object
	 * @param proteinStructure a BioJava Structure object
	 */
	public void setProteinStructure(Structure proteinStructure) {
		this.proteinStructure = proteinStructure;
	}

	/**
	 * Get the pocket structure as BioJava Structure object
	 * @return a BioJava Structure object
	 */
	public Structure getPocketStructure() {
		return pocketStructure;
	}

	/**
	 * Set the pocket structure as BioJava Structure object
	 * @param pocketStructure a BioJava Structure object
	 */
	public void setPocketStructure(Structure pocketStructure) {
		this.pocketStructure = pocketStructure;
	}

	/**
	 * Get the ligand structure as CDK IAtomContainer object
	 * @return a BioJava Structure object
	 */
	public IAtomContainer getLigandStructure() {
		return ligandStructure;
	}

	/**
	 * Set the ligand structure as CDK IAtomContainer object
	 * @param ligandStructure a CDK IAtomContainer object
	 */
	public void setLigandStructure(IAtomContainer ligandStructure) {
		this.ligandStructure = ligandStructure;
	}

	/**
	 * Get SIFTS entities for the current PDBbindEntry
	 * @return a List of SiftsEntity
	 */
	public List<SiftsEntity> getSiftEntities() {
		return siftEntities;
	}

	/**
	 * Set a list of SIFTS entities for the current PDBbindEntry
	 * @param siftEntities a List of SiftsEntity
	 */
	public void setSiftEntities(List<SiftsEntity> siftEntities) {
		this.siftEntities = siftEntities;
	}

	
	/**
	 * Get SIFTS residues for the current PDBbindEntry
	 * @return a List of SiftsResidue
	 */
	public List<SiftsResidue> getSiftResidues() {
		return siftResidues;
	}

	/**
	 * Set a list of SIFTS residues for the current PDBbindEntry
	 * @param siftResidues a List of SiftsResidue
	 */
	public void setSiftResidues(List<SiftsResidue> siftResidues) {
		this.siftResidues = siftResidues;
	}

	
	/**
	 * Get a list of AminoAcid from the pocket structure of a PDBbindEntry
	 * @return a list of AminoAcid
	 */
	public List<AminoAcid> getPocketAminoAcids() {
		return pocketAminoAcids;
	}

	/**
	 * Set a list of AminoAcid from the pocket structure to a PDBbindEntry
	 * @param pocketAminoAcids a list of AminoAcid objects
	 */
	public void setPocketAminoAcids(List<AminoAcid> pocketAminoAcids) {
		this.pocketAminoAcids = pocketAminoAcids;
	}

	/**
	 * Get a list of AminoAcid from the protein structure of a PDBbindEntry
	 * @return a list of AminoAcid
	 */
	public List<AminoAcid> getProteinAminoAcids() {
		return proteinAminoAcids;
	}

	/**
	 * Set a list of AminoAcid from the protein structure to a PDBbindEntry
	 * @param proteinAminoAcids a list of AminoAcid objects
	 */
	public void setProteinAminoAcids(List<AminoAcid> proteinAminoAcids) {
		this.proteinAminoAcids = proteinAminoAcids;
	}

	/**
	 * Get a PDB reader object
	 * @return a PDBFileReader
	 */
	public PDBFileReader getReader() {
		return reader;
	}

	/**
	 * Set a PDB reader object
	 * @param reader a PDBFileReader
	 */
	public void setReader(PDBFileReader reader) {
		this.reader = reader;
	}

	/**
	 * Get PDB ID
	 * @return a String of the PDB ID of current PDBbindEntry
	 */
	public String getPdb() {
		return pdb;
	}

	/**
	 * Set PDB ID
	 * @param pdb a String of the PDB ID of current PDBbindEntry
	 */
	public void setPdb(String pdb) {
		this.pdb = pdb;
	}

	@Override
	public String toString() {
		return "PDBbindEntry [proteinStructure=" + proteinStructure.getName() + ", pocketStructure=" + pocketStructure.getName()
				+ ", PDB ID=" + this.pdb
				+ "]";
	}

}
