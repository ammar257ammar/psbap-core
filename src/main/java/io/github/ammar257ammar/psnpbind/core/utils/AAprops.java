package io.github.ammar257ammar.psnpbind.core.utils;

import java.io.File;
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

public class AAprops {

	public static List<String> getAApropsHeader() {

		List<String> header = new ArrayList<String>();

		CsvParserSettings settings = new CsvParserSettings();

		settings.getFormat().setLineSeparator("\n");
		settings.getFormat().setDelimiter(',');

		settings.setNumberOfRowsToSkip(1);

		CsvParser parser = new CsvParser(settings);

		List<String[]> rows = parser.parseAll(new File("/config/AAprops.csv"));

		for (String[] row : rows) {

			header.add(row[0]);
		}

		return header;
	}

	public static List<AminoAcid> getResidueNeighbours(List<AminoAcid> aminoAcids, AminoAcid residue) {

		List<AminoAcid> neighbours = new ArrayList<AminoAcid>();

		for (AminoAcid aa : aminoAcids) {

			if (Calc.getDistance(residue.getCA(), aa.getCA()) <= 8.0) {
				neighbours.add(aa);
			}
		}

		return neighbours;
	}

	public static List<Double> getResidueAndSurroundingProps(String path, String residueNumber, String WT) throws IOException{
		
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
				
				CsvParserSettings settings = new CsvParserSettings();

				settings.getFormat().setLineSeparator("\n");
				settings.getFormat().setDelimiter(',');

				settings.setNumberOfRowsToSkip(1);
				
				CsvParser parser = new CsvParser(settings);

				List<String[]> rows = parser.parseAll(new File("config/AAprops.csv"));

				List<AminoAcid> neighbours = AAprops.getResidueNeighbours(aminoAcids, residue);

				for(String[] row: rows) {

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
