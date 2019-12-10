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

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.core.util.InputStreamProvider;
import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.Chain;
import org.biojava.nbio.structure.Group;
import org.biojava.nbio.structure.GroupType;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;
import org.biojava.nbio.structure.io.sifts.SiftsEntity;
import org.biojava.nbio.structure.io.sifts.SiftsResidue;
import org.biojava.nbio.structure.io.sifts.SiftsSegment;
import org.biojava.nbio.structure.io.sifts.SiftsXMLParser;
import org.biojava.nbio.structure.io.PDBFileReader;


/**
 * A utility class to perform operations on PDB and SIFTS files related to PdbBind Datasets
 *  
 * @author Ammar Ammar
 * 
 */

public class PdbTools {
	

	/**
	 * A method to configure a PDB file reader 
	 * 
	 * @param alignSeqAndPraseSS boolean value to specify if secondary structure info need to be parsed
	 *        and if SEQRES sequence and ATOM sequence need to be aligned when reading a PDB file
	 * @return a PDB file reader
	 */
	public static PDBFileReader configureReader(boolean alignSeqAndPraseSS){
		
		AtomCache cache = new AtomCache();

    	FileParsingParameters params = cache.getFileParsingParams();
    	
    	params.setAlignSeqRes(alignSeqAndPraseSS);
    	params.setParseSecStruc(alignSeqAndPraseSS);
    	
    	cache.setFileParsingParams(params);
    	
    	PDBFileReader reader = new PDBFileReader();
    	
    	reader.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
    	
    	reader.setFileParsingParameters(params);
    	
    	return reader;
	}
	

	/**
	 * A static method to extract AminoAcids from a BioJava Structure object and store the in a List
	 * 
	 * @param structure as a BioJava Structure object
	 * @return a list of AminoAcids extracted from the structure
	 */
	public static List<AminoAcid> getAminoAcidsFromStructure(Structure structure){
		
		List<AminoAcid> aaList = new ArrayList<AminoAcid>();
		
		for(Chain chain: structure.getChains()) {
        	
			List<Group> groups = chain.getAtomGroups(GroupType.AMINOACID);
			
			for (Group group : groups) {
		    
				AminoAcid aa = (AminoAcid) group;
				
				aaList.add(aa);
				
			}
			
		}
		return aaList;
	}
	
	/**
	 * Get a list of Sifts entities from the SIFTS file corresponding to a PdbBind PDB ID
	 * @param path of the downloaded SIFTS files
	 * @param pdb which is the ID for the PDB needed to get its SIFTS
	 * @return a list of SIFTS entities for a specfic PDB ID
	 */
	public static List<SiftsEntity> getSiftsEntitiesForPDB(String path,String pdb){
		
		InputStreamProvider prov = new InputStreamProvider();
		InputStream is;
		
		try {
			
			is = prov.getInputStream(path + "/" + pdb + ".xml.gz");
		
			SiftsXMLParser siftParser = new SiftsXMLParser();
	
			siftParser.parseXmlFile(is);
	
			List<SiftsEntity> entities = siftParser.getEntities();
			
			return entities;
		} catch (IOException e) {
			return null;
		}		
	} 
	
	/**
	 * A method to extract SIFTS residues from SIFTS entities list and store them in a List
	 * @param siftsEntities a List of entities extracted from a SIFTS file for a certain PDB
	 * @return a list of SIFTS residues
	 */
	public static List<SiftsResidue> getSiftResiduesFromSiftsMapping(List<SiftsEntity> siftsEntities){
		
		List<SiftsResidue> srList = new ArrayList<SiftsResidue>();
		
		for (SiftsEntity entity : siftsEntities){

			if(entity.getType().equals("protein")){
				
	            for ( SiftsSegment seg: entity.getSegments()) {
	                
	                for ( SiftsResidue res: seg.getResidues() ) {
	                	
	                	srList.add(res);
	                
	                }
	            }
			}
		}
		return srList;
	}
	

}
