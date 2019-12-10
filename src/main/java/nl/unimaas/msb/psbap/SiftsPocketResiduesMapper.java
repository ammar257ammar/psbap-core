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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.io.sifts.SiftsResidue;

import com.univocity.parsers.csv.CsvParser;
import com.univocity.parsers.csv.CsvParserSettings;

import nl.unimaas.msb.psbap.model.PDBbindEntry;
import nl.unimaas.msb.psbap.utils.PdbTools;

/**
 * This class maps the pocket residues for PdbBind files to the UniProt variants and keep only
 * those variant which occur in the binding pocker
 *  
 * @author Ammar Ammar
 *
 */
public class SiftsPocketResiduesMapper {
	
	
	/**
	 * A method to map pocket residues of PdbBind entries to UniProt variants and return a list of the results
	 * 
	 * @param path of the variants file
	 * @param pdbbindEntriesFile path of the PdbBind dataset file
	 * @param mappingType wich is wither "protein" or "pocket"
	 * @return a dataset as list of string arrays
	 */
	public static List<String[]> mapPocketResidues(String path, String pdbbindEntriesFile, String mappingType){

    	CsvParserSettings settings = new CsvParserSettings();
		
		settings.getFormat().setLineSeparator("\n");
		settings.getFormat().setDelimiter('\t');

		CsvParser parser = new CsvParser(settings);
		
		List<String[]> rows = parser.parseAll(new File(path));
		
		List<String[]> pdbVariantsRows = new ArrayList<String[]>();
						
		Map<String, PDBbindEntry> pdbbindEntries = PdbTools.parsePDBbindEntriesFromFile(pdbbindEntriesFile);
						
		for (String[] row: rows){
						
			String pdb = row[4];
			String snp = row[2];
			
			PDBbindEntry pbEntry = pdbbindEntries.get(pdb);
			
			if(pbEntry != null &&
			   pbEntry.hasPocketStructure()  &&
			   pbEntry.hasSiftsMapping()){

				List<AminoAcid> aaList = null;
				
				if(mappingType.equals("pocket")){
					
					aaList = pbEntry.getPocketAminoAcids();
					
				}else if(mappingType.equals("protein")){
					
					aaList = pbEntry.getProteinAminoAcids();
					
				}
										
				String snpResidue = snp.substring(2);
				
				String sourceAminoAcid = snpResidue.substring(0, 3);
				String targetAminoAcid = snpResidue.substring(snpResidue.length()-3, snpResidue.length());
				
				String residueNum = snpResidue.substring(3, snpResidue.length()-3);
		    	
		    					
				for(AminoAcid aa: aaList) {
					
					for(SiftsResidue res: pbEntry.getSiftResidues())
		            	
	                    if(res.getPdbResName().equals(aa.getPDBName()) && 
	                       res.getPdbResNum().equals(aa.getResidueNumber().getSeqNum().toString()) &&
	                       res.getUniProtAccessionId() != null &&
	                       !res.getUniProtAccessionId().equals("null") &&
	                       res.getUniProtPos().toString().equals(residueNum) &&
	                       res.getPdbResName() != null) {
						           			                    	
	                    	List<String> pdbTempRow = new ArrayList<String>();
	                    	
	                    	pdbTempRow.add(row[0]);
	                    	pdbTempRow.add(row[1]);
	                    	pdbTempRow.add(row[2]);
	                    	pdbTempRow.add(row[3]);
	                    	pdbTempRow.add(row[4]);
	                    	pdbTempRow.add(row[5]);
	                    	pdbTempRow.add(row[6]);
	                    	
	                    	pdbTempRow.add((StructureTools.get1LetterCodeAmino(sourceAminoAcid.toUpperCase()) != null) ? StructureTools.get1LetterCodeAmino(sourceAminoAcid.toUpperCase()).toString() : sourceAminoAcid);
	                    	pdbTempRow.add((StructureTools.get1LetterCodeAmino(targetAminoAcid.toUpperCase()) != null) ? StructureTools.get1LetterCodeAmino(targetAminoAcid.toUpperCase()).toString() : "");
	                    	pdbTempRow.add(residueNum);
	                    	pdbTempRow.add(res.getChainId());

	                    	try {
	                        	pdbTempRow.add((res.getPdbResName() != null) ? StructureTools.get1LetterCodeAmino(res.getPdbResName()).toString() : "");
	                    		
	                    	}catch(NullPointerException ex) {
	                    		pdbTempRow.add("");
	                    	}
	                    	
	                    	pdbTempRow.add((res.getPdbResNum() != null) ? res.getPdbResNum() : "");

	                    	try {
	                        	pdbTempRow.add((aa.getPDBName() != null) ? StructureTools.get1LetterCodeAmino(aa.getPDBName()).toString() : "");
	                    		
	                    	}catch(NullPointerException ex) {
	                    		pdbTempRow.add("");
	                    	}
	                    	
	                    	pdbTempRow.add((aa.getResidueNumber() != null) ? aa.getResidueNumber().getSeqNum().toString() : "");
	                    	pdbTempRow.add((res.getUniProtResName() != null) ? res.getUniProtResName() : "");
	                    	pdbTempRow.add((res.getUniProtPos() != null) ? res.getUniProtPos().toString() : "");
	                    	pdbTempRow.add((res.getUniProtAccessionId() != null) ? res.getUniProtAccessionId() : "");
	                    	pdbTempRow.add((res.getPdbId() != null) ? res.getPdbId() : "");
	                    	pdbTempRow.add((res.getSeqResName() != null) ? res.getSeqResName() : "");
	                    	pdbTempRow.add((res.getNaturalPos() != null) ? res.getNaturalPos().toString() : "");
	            			
	                    	if(pdbTempRow.get(7).equals(pdbTempRow.get(13)) && pdbTempRow.get(7).equals(pdbTempRow.get(15)) && !pdbTempRow.get(8).equals("")) {
	                    		pdbVariantsRows.add(pdbTempRow.toArray(new String[pdbTempRow.size()]));
	                    	}
						}
					}
				}
		}
			
    	return pdbVariantsRows;
    }

}
