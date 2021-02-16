/**
* Binding Pocket SNPs' effect on Binding Affinity Database Project (PSnpBind)
* 
*Copyright (C) 2019-2021  Ammar Ammar <ammar257ammar@gmail.com> ORCID:0000-0002-8399-8990
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

package io.github.ammar257ammar.psnpbind.core;

import java.io.IOException;
import java.util.List;

import io.github.ammar257ammar.psnpbind.core.Config;
import io.github.ammar257ammar.psnpbind.core.FoldX;
import io.github.ammar257ammar.psnpbind.core.SiftsPocketResiduesMapper;
import io.github.ammar257ammar.psnpbind.core.UniProtVariantsMapper;
import io.github.ammar257ammar.psnpbind.core.db.DbDataFabricator;
import io.github.ammar257ammar.psnpbind.core.model.PdbBindDataset;
import io.github.ammar257ammar.psnpbind.core.model.PdbBindDataset.PdbbindAttribute;
import io.github.ammar257ammar.psnpbind.core.utils.DataHandler;
import io.github.ammar257ammar.psnpbind.core.utils.PdbTools;

/**
 * This class perform the tasks passed through command line
 * 
 * @author Ammar Ammar
 *
 */
public class PSnpBindCore 
{
    public static void main( String[] args )
    {

    	CliOptions cli = new CliOptions(args);
    	
    	switch(cli.operation){
    	
    	case "init":
    		
    		PdbBindDataset pdbbindData = PdbBindDataset.create().
											loadData().
											keepAsFolderMatch().
											filterStringNotEqual(PdbbindAttribute.UNIPROT, "------").
											filterStringNotEqual(PdbbindAttribute.RESOLUTION, "NMR").
											sortBy(PdbbindAttribute.RESOLUTION).
											groupByUniProtAndKeepMinResolution(true).
											filterDoubleCutoff(PdbbindAttribute.RESOLUTION, 2.51);
    		
    	    DataHandler.writeDatasetToTSV(pdbbindData.getData(), Config.getProperty("TSV_PATH") + "/pdbbind_entries_data.tsv");

    	    DataHandler.printDatasetStats(pdbbindData.getData());
    	    DataHandler.printDatasetHead(pdbbindData.getData());

    	    List<String[]> pdbbindDataSiftsUrls = pdbbindData.asSiftsDownloadUrlsList();    	  	
        	DataHandler.writeDatasetToTSV(pdbbindDataSiftsUrls, Config.getProperty("TSV_PATH") + "/pdbbind_sifts_urls.tsv");
    		
        	try {
				String siftsDownloadStatus = PdbTools.downloadSifts(Config.getProperty("TSV_PATH") + "/pdbbind_sifts_urls.tsv",
						Config.getProperty("SIFTS_PATH"));

				System.out.println(siftsDownloadStatus);

			} catch (IOException | InterruptedException e1) {
				e1.printStackTrace();
			}

        	List<String[]> pdbbindVariants = UniProtVariantsMapper.mapMissenseVariantsToPdbbindDataset(pdbbindData.getData());
	    	
    		DataHandler.writeDatasetToTSV(pdbbindVariants, Config.getProperty("TSV_PATH") + "/pdbbind_protein_variants.tsv");  	
	    	
    		System.gc(); 
    		
    		break;
    		
    	case "pocket-snps-mapping-and-foldx-prep":
    		
    		List<String[]> pdbbindPocketVariants = SiftsPocketResiduesMapper.
    				mapPocketResidues(Config.getProperty("TSV_PATH") + "/pdbbind_protein_variants.tsv",
    								  Config.getProperty("TSV_PATH") + "/pdbbind_entries_data.tsv",
    								  "pocket", true);

    		String[] header = new String[]{	"uniprot",
											"snp",
											"rs_id",
											"mutation_type",
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
    				
    		DataHandler.writeDatasetToTSV(pdbbindPocketVariants, Config.getProperty("TSV_PATH") + "/pdbbind_pocket_variants.tsv", header);

    		try {
	    		
    			FoldX.createMutationConfigFiles(pdbbindPocketVariants, Config.getProperty("PDBBIND_ENTRIES_PATH"), Config.getProperty("FOLDX_PDB_DIR"));
			
	    	} catch (IOException e) {
				e.printStackTrace();
			}

			break;
					
    	case "foldx-report":
    		
    		List<String[]> repairResults2 = FoldX.getFoldxResults(Config.getProperty("FOLDX_PDB_DIR"));
        	DataHandler.writeDatasetToTSV(repairResults2, Config.getProperty("TSV_PATH") + "/foldx_results.tsv");
        	
    		List<String[]> mutationResults = FoldX.buildFoldxReport(Config.getProperty("FOLDX_PDB_DIR"));
        	
    		DataHandler.writeDatasetToTSV(mutationResults, Config.getProperty("TSV_PATH") + "/foldx_mutation_results.tsv", 
        			new String[] {"PDB", "Mutation", "Energy", "SD"});
    		
        	break;
        	
    	case "prepare-ligands-folders":
    		
        	try {
				Ligand3D.prepareLigandsFolder(Config.getProperty("FOLDX_PDB_DIR"), 
						Config.getProperty("PDBBIND_ENTRIES_PATH"), 
						Config.getProperty("LIGANDS_PATH"));
			} catch (IOException e) {
				e.printStackTrace();
			}
        	
        	break;
        	
    	case "ligands-tanimoto-dataset":
    		    		
			try {
				List<String[]> similarLigands = Ligand3D.getLigandsIDsFiltered(Config.getProperty("LIGANDS_PATH"), true);
	        	
				DataHandler.writeDatasetToTSV(similarLigands, 
	        			Config.getProperty("TSV_PATH") + "/chembl_ligands_filtered.tsv");

	        	List<String[]> ligandsWithIDsAndTanimoto = Ligand3D.combineIDsAndTanimotoOfLigands(Config.getProperty("LIGANDS_PATH"), 
	        			Config.getProperty("TSV_PATH") + "/chembl_ligands_filtered.tsv", true);

	        	DataHandler.writeDatasetToTSV(ligandsWithIDsAndTanimoto, 
	        			Config.getProperty("TSV_PATH") + "/chembl_ligands_filtered_combined_tanimoto.tsv");
	        	
	        	
	        	List<String[]> ligandsWithIDsAndTanimotoAndSmiles = Ligand3D.combineIDsAndTanimotoAndSmilesOfLigands(Config.getProperty("LIGANDS_PATH"), 
	        			Config.getProperty("TSV_PATH") + "/chembl_ligands_filtered.tsv", true);
	        	
	        	DataHandler.writeDatasetToTSV(ligandsWithIDsAndTanimotoAndSmiles, 
	        			Config.getProperty("TSV_PATH") + "/chembl_ligands_filtered_combined_tanimoto_smiles.tsv");
        	
			} catch (IOException e) {
				e.printStackTrace();
			}

			break;
        	
    	case "prepare-vina-folders-config":
    		
    		try {
    			Vina.createFolderStructureFromGromacsFolder(Config.getProperty("FOLDX_PDB_DIR"), Config.getProperty("VINA_DOCKING_DIR"));

    			Vina.createLigandsStructureAfterVinaFolder(Config.getProperty("VINA_DOCKING_DIR"), Config.getProperty("LIGANDS_PATH"), true);

    			Vina.createSeedFilesAfterLigandsStructure(Config.getProperty("VINA_DOCKING_DIR"), "1264647227", false);

    			Vina.createConfigFilesAfterLigandsStructure(Config.getProperty("VINA_DOCKING_DIR"),false);       		
       		
			} catch (IOException e) {
				e.printStackTrace();
			}
    		
    		
    		break;
  
    	case "generate-dockings-results":
    		
    		try {
    	    	Vina.generateVinaReportAll(Config.getProperty("VINA_DOCKING_DIR"),Config.getProperty("TSV_PATH")+"/docking-results");
			} catch (IOException e) {
				e.printStackTrace();
			}
    		
    		break;
    		
    	case "build-database":
    		
    		DbDataFabricator.sanityChecks();
    		DbDataFabricator.populateDBWithProteinsData();
    		DbDataFabricator.populateDBWithVariantsData();
    		DbDataFabricator.populateDBWithLigandsData();
    		
    		break;
    	}
    }
}
