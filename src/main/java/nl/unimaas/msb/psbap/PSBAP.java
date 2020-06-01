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

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.biojava.nbio.structure.StructureException;
import org.openscience.cdk.exception.CDKException;

import nl.unimaas.msb.psbap.model.PdbBindDataset.PdbbindAttribute;
import nl.unimaas.msb.psbap.Config;
import nl.unimaas.msb.psbap.FoldX;
import nl.unimaas.msb.psbap.SiftsPocketResiduesMapper;
import nl.unimaas.msb.psbap.UniProtVariantsMapper;
import nl.unimaas.msb.psbap.model.PdbBindDataset;
import nl.unimaas.msb.psbap.utils.DataHandler;
import nl.unimaas.msb.psbap.utils.PdbTools;

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
    	
    	case "init":
    		
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

    	    List<String[]> pdbbindDataSiftsUrls = pdbbindData.asSiftsDownloadUrlsList();    	  	
        	DataHandler.writeDatasetToTSV(pdbbindDataSiftsUrls, Config.getProperty("DATASETS_PATH") + "/pdbbind_sifts_urls.tsv");

    	    List<String[]> pdbbindDataFastaUrls = pdbbindData.asFastaDownloadUrlsList();    	  	
        	DataHandler.writeDatasetToTSV(pdbbindDataFastaUrls, Config.getProperty("DATASETS_PATH") + "/pdbbind_fasta_urls.tsv");
    		
    	    List<String[]> pdbbindDataDsspUrls = pdbbindData.asDsspDownloadUrlsList();    	  	
        	DataHandler.writeDatasetToTSV(pdbbindDataDsspUrls, Config.getProperty("DATASETS_PATH") + "/pdbbind_dssp_urls.tsv");
    		
        	try {
				PdbTools.downloadSifts(Config.getProperty("DATASETS_PATH") + "/pdbbind_sifts_urls.tsv",
						Config.getProperty("SIFTS_PATH") + "/sifts");
				
				PdbTools.downloadSifts(Config.getProperty("DATASETS_PATH") + "/pdbbind_fasta_urls.tsv",
						Config.getProperty("FASTA_PATH") + "/FASTA");
				
				PdbTools.downloadSifts(Config.getProperty("DATASETS_PATH") + "/pdbbind_dssp_urls.tsv",
						Config.getProperty("DSSP_PATH") + "/DSSP");
				
			} catch (IOException | InterruptedException e1) {
				e1.printStackTrace();
			}

        	List<String[]> pdbbindVariants = UniProtVariantsMapper.mapMissenseVariantsToPdbbindDataset(pdbbindData.getData());
	    	
    		DataHandler.writeDatasetToTSV(pdbbindVariants, Config.getProperty("DATASETS_PATH") + "/pdbbind_protein_variants.tsv");  	
	    	
    		System.gc(); 
    		
    		break;
    		
    	case "pocket-snps-mapping-and-foldx-prep":
    		
    		List<String[]> pdbbindPocketVariants = SiftsPocketResiduesMapper.
    				mapPocketResidues(Config.getProperty("DATASETS_PATH") + "/pdbbind_protein_variants.tsv",
    								  Config.getProperty("DATASETS_PATH") + "/pdbbind_entries_data.tsv",
    								  "pocket", true);

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


    		try {
	    		
    			FoldX.createMutationConfigFiles(pdbbindPocketVariants, Config.getProperty("PDBBIND_ENTRIES_PATH"), Config.getProperty("FOLDX_PDB_DIR"));
			
	    	} catch (IOException e) {
				e.printStackTrace();
			}

			break;
					
    	case "foldx-report":
    		
    		List<String[]> repairResults2 = FoldX.getFoldxResults(Config.getProperty("FOLDX_PDB_DIR"));
        	DataHandler.writeDatasetToTSV(repairResults2, Config.getProperty("DATASETS_PATH") + "/foldx_results.tsv");
        	
    		List<String[]> mutationResults = FoldX.buildFoldxReport(Config.getProperty("FOLDX_PDB_DIR"));
        	
    		DataHandler.writeDatasetToTSV(mutationResults, Config.getProperty("DATASETS_PATH") + "/foldx_mutation_results.tsv", 
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
				List<String[]> similarLigands = Ligand3D.getLigandsIDsFiltered(Config.getProperty("LIGANDS_PATH"));
	        	DataHandler.writeDatasetToTSV(similarLigands, 
	        			Config.getProperty("DATASETS_PATH") + "/chembl_ligands_filtered.tsv");

	        	List<String[]> ligandsWithIDsAndTanimoto = Ligand3D.combineIDsAndTanimotoOfLigands(Config.getProperty("LIGANDS_PATH"), 
	        			Config.getProperty("DATASETS_PATH") + "/chembl_ligands_filtered.tsv");
	        	
	        	DataHandler.writeDatasetToTSV(ligandsWithIDsAndTanimoto, 
	        			Config.getProperty("DATASETS_PATH") + "/chembl_ligands_filtered_combined_tanimoto.tsv");
	        	
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
    		
    	case "generate-features":
    		
    		try {
    			
    			File vinaFolder = new File("/processing/vina-docking");
    			
    			File[] pdbs = vinaFolder.listFiles();
    			
    			for(File pdbFile: pdbs) {
    	    	  	
    				String pdb = pdbFile.getName();
    				
    	    		if(!new File("/features/"+pdb).exists()) {
    	    			new File("/features/"+pdb).mkdir();
    	    		}
    	    		
    	    		List<String[]> featuresOnePDB = Featurizer.getSnpsFeatures("/tsv/pdbbind_pocket_variants.tsv",Config.getProperty("FOLDX_PDB_DIR"),pdb);
    	    		DataHandler.writeDatasetToTSV(featuresOnePDB, "/features/"+pdb+"/pdbbind_pocket_variants_features_"+pdb+".tsv");
    	        	  
    	        	List<String[]> featuresOnePDBL = Featurizer.getLigandsFeatures(Config.getProperty("LIGANDS_PATH"),"/tsv/chembl_ligands_filtered_combined_tanimoto.tsv",pdb);
    				DataHandler.writeDatasetToTSV(featuresOnePDBL, "/features/"+pdb+"/chembl_ligands_features_"+pdb+".tsv");
    	    	
    	        	List<String[]> featuresOnePocket = Featurizer.getPocketsFeatures(pdb);
    				DataHandler.writeDatasetToTSV(featuresOnePocket, "/features/"+pdb+"/pdbbind_pocket_features_"+pdb+".tsv");	
    	    		
    	    		Vina.generateVinaReport(Config.getProperty("VINA_DOCKING_DIR"),"/features/"+pdb+"/bindingAffinity-official-"+pdb, pdb);
    	    	}
    			
			} catch (IOException e) {
				e.printStackTrace();
			} catch (StructureException e) {
				e.printStackTrace();
			} catch (ClassNotFoundException e) {
				e.printStackTrace();
			} catch (CDKException e) {
				e.printStackTrace();
			} catch (CloneNotSupportedException e) {
				e.printStackTrace();
			}
    		
    		
    		break;
    	}	
    }
}
