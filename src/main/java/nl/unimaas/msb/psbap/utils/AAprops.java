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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.biojava.nbio.structure.AminoAcid;
import org.biojava.nbio.structure.Calc;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.io.PDBFileReader;

import com.univocity.parsers.csv.CsvParser;
import com.univocity.parsers.csv.CsvParserSettings;

/**
 * A utility class to use 48 amino acid properties in creating features for AA residues and their surrounding
 * It can calculate the AA properties for all residues surrounding a specific residue with 8A radius
 * and return the results as a list of double values
 * 
 * @author Ammar Ammar
 * 
 */
public class AAprops {
			
	List<String[]> rows = null;
	
	public AAprops() {

		CsvParserSettings settings = new CsvParserSettings();

		settings.getFormat().setLineSeparator("\n");
		settings.getFormat().setDelimiter(',');

		settings.setNumberOfRowsToSkip(1);

		CsvParser parser = new CsvParser(settings);
		
        this.rows = parser.parseAll(getClass().getResourceAsStream("/aaprops.csv"));	   
	}
	
	/**
	 * create a header from the AA properties names in the resource CSV file
	 * @return the header as a list of strings
	 */
	public List<String> getAApropsHeader() {

		List<String> header = new ArrayList<String>();


		for (String[] row : this.rows) {

			header.add(row[0]);
		}

		return header;
		
	}
	
	/**
	 * A method to get a list of AA neighbours of a specific amino acid within 8A radius sphere 
	 * @param aminoAcids list of amino acids from the protein structure
	 * @param residue the amino acid to get its neighbours
	 * @return list of neighbour residues
	 */
	public List<AminoAcid> getResidueNeighbours(List<AminoAcid> aminoAcids, AminoAcid residue) {

		List<AminoAcid> neighbours = new ArrayList<AminoAcid>();

		for (AminoAcid aa : aminoAcids) {

			if (Calc.getDistance(residue.getCA(), aa.getCA()) <= 8.0) {
				neighbours.add(aa);
			}
		}

		return neighbours;
	}
	

	/**
	 * A method to calculate the properties of a specific residue AA and its surrounding environment
	 * @param path a path string for the PDB of the protein
	 * @param residueNumber the targeted residue number
	 * @param WT the wild-type residue name
	 * @return list of 96 double values (48 for residue properties, 48 for surrounding properties)
	 */
	public List<Double> getResidueAndSurroundingProps(String path, String residueNumber, String WT) throws IOException{
		
		List<Double> props = new ArrayList<Double>();
		List<Double> nbProps = new ArrayList<Double>();
		
		HashMap<Character,Integer> map = new HashMap<Character,Integer>();
		map.put('A', 1);
		map.put('R', 15);
		map.put('N', 12);
		map.put('D', 2);
		map.put('C', 3);
		map.put('Q', 14);
		map.put('E', 4);
		map.put('G', 6);
		map.put('H', 7);
		map.put('I', 8);
		map.put('L', 10);
		map.put('K', 9);
		map.put('M', 11);
		map.put('F', 5);
		map.put('P', 13);
		map.put('S', 16);
		map.put('T', 17);
		map.put('W', 19);
		map.put('Y', 20);
		map.put('V', 18);
		
		PDBFileReader reader = PdbTools.configureReader(false);

		Structure proteinStructure = reader.getStructure(path);

		List<AminoAcid> aminoAcids = PdbTools.getAminoAcidsFromStructure(proteinStructure);
		
		for (int i = 0; i < aminoAcids.size(); i++) {

			if (aminoAcids.get(i).getResidueNumber().toString().equals(residueNumber)) {
				
				AminoAcid residue = aminoAcids.get(i);
				
				List<AminoAcid> neighbours = this.getResidueNeighbours(aminoAcids, residue);

				for(String[] row: this.rows) {

					double aaPropValue = Double.parseDouble(row[map.get(residue.getAminoType())]);
					double aaWTPropValue = Double.parseDouble(row[map.get(WT.trim().charAt(0))]);

					double neighboursPropValue = 0.0;

					props.add(aaPropValue-aaWTPropValue);
					
					if(neighbours.size() > 0) {
						
						for(AminoAcid nb: neighbours) {
							neighboursPropValue += Double.parseDouble(row[map.get(nb.getAminoType())]);
						}
						
						nbProps.add(neighboursPropValue - aaPropValue);

					}else {
						nbProps.add(0.0);
					}
					
				}
				
				props.addAll(nbProps);
				
				break;

			}
		}
		
		return props;
	}

}