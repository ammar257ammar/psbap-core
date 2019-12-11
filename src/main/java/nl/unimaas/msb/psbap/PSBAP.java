/**
* binding Pocket's SNPS effect on Binding Affinity Project (PSBAP) 
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

import java.io.IOException;
import java.util.List;

import nl.unimaas.msb.psbap.model.PdbBindDataset.PdbbindAttribute;
import nl.unimaas.msb.psbap.Config;
import nl.unimaas.msb.psbap.FoldX;
import nl.unimaas.msb.psbap.SiftsPocketResiduesMapper;
import nl.unimaas.msb.psbap.UniProtVariantsMapper;
import nl.unimaas.msb.psbap.model.PdbBindDataset;
import nl.unimaas.msb.psbap.utils.DataHandler;

/**
 * This class perform the tasks passed through command line
 * 
 * @author Ammar Ammar
 *
 */
public class PSBAP 
{
    public static void main( String[] args )
    {

    	CliOptions cli = new CliOptions(args);
    	
    	switch(cli.operation){
    	
    	case "generate-pdbbind-dataset":
    		
    		PdbBindDataset pdbbindData = PdbBindDataset.create().
											loadData().
											filterStringNotEqual(PdbbindAttribute.UNIPROT, "------").
											filterStringNotEqual(PdbbindAttribute.RESOLUTION, "NMR").
											sortBy(PdbbindAttribute.RESOLUTION).
											keepAsFolderMatch().
											filterDoubleCutoff(PdbbindAttribute.RESOLUTION, 2.51).
											groupByUniProtAndKeepMinResolution(true);
    		
    	    DataHandler.writeDatasetToTSV(pdbbindData.getData(), Config.getProperty("DATASETS_PATH") + "/pdbbind_entries_data.tsv");

    	    DataHandler.printDatasetStats(pdbbindData.getData());

    	    DataHandler.printDatasetHead(pdbbindData.getData());
    	    
			break;
			
    	case "generate-sifts-urls":
    				
    		PdbBindDataset pdbbindDataForSifts = PdbBindDataset.create().
											loadData().
											filterStringNotEqual(PdbbindAttribute.UNIPROT, "------").
											filterStringNotEqual(PdbbindAttribute.RESOLUTION, "NMR").
											sortBy(PdbbindAttribute.RESOLUTION).
											keepAsFolderMatch().
											filterDoubleCutoff(PdbbindAttribute.RESOLUTION, 2.51).
											groupByUniProtAndKeepMinResolution(true);

        	List<String[]> pdbbindDataSiftsUrls = pdbbindDataForSifts.asSiftsDownloadUrlsList();    	
        	
        	DataHandler.writeDatasetToTSV(pdbbindDataSiftsUrls, Config.getProperty("DATASETS_PATH") + "/pdbbind_sifts_urls.tsv");
    		
    		break;
    		
    	case "map-uniprot-variant-to-proteins":
    		
    		PdbBindDataset pdbbindDataForVariants = PdbBindDataset.create().
													loadData().
													filterStringNotEqual(PdbbindAttribute.UNIPROT, "------").
													filterStringNotEqual(PdbbindAttribute.RESOLUTION, "NMR").
													sortBy(PdbbindAttribute.RESOLUTION).
													keepAsFolderMatch().
													filterDoubleCutoff(PdbbindAttribute.RESOLUTION, 2.51).
													groupByUniProtAndKeepMinResolution(true);
    		
    		List<String[]> pdbbindVariants = UniProtVariantsMapper.mapMissenseVariantsToPdbbindDataset(pdbbindDataForVariants.getData());
	    	
    		DataHandler.writeDatasetToTSV(pdbbindVariants, Config.getProperty("DATASETS_PATH") + "/pdbbind_protein_variants.tsv");  	
	    	
    		System.gc();
    		
    		break;
    		
    	case "map-pocket-residues-to-uniprot-variants":

    		List<String[]> pdbbindPocketVariants = SiftsPocketResiduesMapper.
			mapPocketResidues(Config.getProperty("DATASETS_PATH") + "/pdbbind_protein_variants.tsv",
							  Config.getProperty("DATASETS_PATH") + "/pdbbind_entries_data.tsv",
							  "pocket");

			String[] header = new String[]{"gene_name",
											"uniprot",
											"snp",
											"rs_id",
											"pdb",
											"min_res",
											"ligand",
											"sourceAminoAcid",
											"targetAminoAcid",
											"residueNum",
											"chain",
											"PdbResName",
											"PdbResNum",
											"aaPDBName",
											"aaResidueNumber",
											"UniProtResName",
											"UniProtPos",
											"UniProtAccessionId",
											"PdbId",
											"SeqResName",
											"NaturalPos"};
			
			DataHandler.writeDatasetToTSV(pdbbindPocketVariants, Config.getProperty("DATASETS_PATH") + "/pdbbind_pocket_variants.tsv", header);

    		break;
    		
    	case "generate-foldx-mutations-and-structure":
    		
	    	try {
	    		
	    		List<String[]> pdbbindPocketVariantsForFoldX = SiftsPocketResiduesMapper.
	    				mapPocketResidues(Config.getProperty("DATASETS_PATH") + "/pdbbind_protein_variants.tsv",
	    								  Config.getProperty("DATASETS_PATH") + "/pdbbind_entries_data.tsv",
	    								  "pocket");
	    		
				FoldX.createMutationConfigFiles(pdbbindPocketVariantsForFoldX, Config.getProperty("PDBBIND_ENTRIES_PATH"), Config.getProperty("FOLDX_PDB_DIR"));
			
	    	} catch (IOException e) {
				e.printStackTrace();
			}

    		break;
    		
    	case "foldx-success-report":
    		
    		List<String[]> repairResults2 = FoldX.getFoldxResults(Config.getProperty("FOLDX_PDB_DIR"));
        	DataHandler.writeDatasetToTSV(repairResults2, Config.getProperty("DATASETS_PATH") + "/foldx_results.tsv");
        	
        	break;
        	
    	case "foldx-energies-report":
    		
    		List<String[]> mutationResults = FoldX.buildFoldxReport(Config.getProperty("FOLDX_PDB_DIR"));
        	
    		DataHandler.writeDatasetToTSV(mutationResults, Config.getProperty("DATASETS_PATH") + "/foldx_mutation_results.tsv", 
        			new String[] {"PDB", "Mutation", "Energy", "SD"});

    		break;
    		
    	case "foldx-energies-html-report":
    		
    
    		break;
    	}	
    }
}
