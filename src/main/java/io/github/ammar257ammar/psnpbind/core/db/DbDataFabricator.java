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

package io.github.ammar257ammar.psnpbind.core.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.URL;
import java.net.URLConnection;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.UUID;

import org.apache.ibatis.jdbc.ScriptRunner;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.qsar.result.DoubleResult;

import com.univocity.parsers.csv.CsvParser;
import com.univocity.parsers.csv.CsvParserSettings;

import io.github.ammar257ammar.psnpbind.core.Config;
import io.github.ammar257ammar.psnpbind.core.utils.LigandTools;

/**
 * A class to prepare populate the database with PSnPBind Data.
 * 
 * @author Ammar Ammar
 *
 */
public class DbDataFabricator {
	
  	
    /**
     * Create the table structure of the database
     * 
     */
	public static void sanityChecks() {
		
		Connection con = null;
		
		try {
			
			con = DatabaseConnection.getConnection();

			ScriptRunner sr = new ScriptRunner(con);

			Reader reader = new BufferedReader(new FileReader("/config/schema.sql"));
			sr.runScript(reader);
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}finally{
			try {
				con.close();
			} catch (SQLException e) {
				e.printStackTrace();
			}
		} 
	}
	
	/**
     * Insert Ligand data and docking results into the database. In the process, 
     * 5 descriptors of each ligand are computed (number of rotatable bonds,
     * hydrogen donor count, hydrogen acceptor count, molecular weight, XLogP)
     * 
     */
	public static void populateDBWithLigandsData() {
		
		CsvParserSettings settings = new CsvParserSettings();
		
		settings.getFormat().setLineSeparator("\n");
		settings.getFormat().setDelimiter('\t');
		settings.setNumberOfRowsToSkip(1);

		CsvParser parser = new CsvParser(settings);
		
		List<String[]> ligandsRows = parser.parseAll(new File(Config.getProperty("TSV_PATH")+"/chembl_ligands_filtered_combined_tanimoto_smiles.tsv"));
		List<String[]> dockingsRows = parser.parseAll(new File(Config.getProperty("TSV_PATH")+"/docking-results/docking-results-all.tsv"));
				
		Connection con = null;
		PreparedStatement psL = null;
		PreparedStatement psPL = null;
		PreparedStatement psVL = null;
		Statement stmt = null;
		String queryL = "INSERT INTO psnpbind_ligand values (?,?,?,?,?,?,?,?)";
		String queryPL = "INSERT INTO psnpbind_protein_ligand values (?,?,?,?)";
		String queryVL = "INSERT INTO psnpbind_variant_ligand values (?,?,?,?,?,?,?)";
		String selectVariants = "SELECT variant_id, variant_folder, pdb_id, variant_type from psnpbind_variant";
		
		RotatableBondsCountDescriptor descRot = new RotatableBondsCountDescriptor();
		HBondDonorCountDescriptor descHD = new HBondDonorCountDescriptor();
		HBondAcceptorCountDescriptor descHA = new HBondAcceptorCountDescriptor();
		WeightDescriptor descW = new WeightDescriptor();
		XLogPDescriptor descXlogp = new XLogPDescriptor();
		
		Map<String, String[]> ligands = new HashMap<String, String[]>();
		Map<String, String[]> variants = new HashMap<String, String[]>();
		Set<String> ligandsUnique = new HashSet<String>();
		Set<String> proteinLigandsUnique = new HashSet<String>();
		
		for(String[] ligandRow: ligandsRows){
			
			ligands.put(ligandRow[1], new String[] {ligandRow[2], ligandRow[4]});
		}
		
		
		try {
			
			con = DatabaseConnection.getConnection();
			psL = con.prepareStatement(queryL);
			psPL = con.prepareStatement(queryPL);
			psVL = con.prepareStatement(queryVL);
			stmt = con.createStatement();
		      
			con.setAutoCommit(false);
			
		    ResultSet rs = stmt.executeQuery(selectVariants);
		    
		    while(rs.next()){
		    	String variantFolder = rs.getString("pdb_id") + "_protein_Repair_" + rs.getString("variant_folder");
		    	String uuidBaseString = rs.getString("pdb_id") + rs.getString("variant_type");
		    	
		    	variants.put(variantFolder, new String[] {rs.getString("variant_id") , uuidBaseString});
		    }  
		    
		    rs.close();
		    

			int proteinLigandId = 1;
			int variantLigandId = 1;
			int variantLigandTotal = dockingsRows.size();
			
			for(String[] dockingRow: dockingsRows) {
				
				String[] variantData = variants.get(dockingRow[1].trim());
				String[] ligandData = ligands.get(dockingRow[2].trim());
				
				if(ligandData == null) {
					variantLigandTotal -= 1;
					continue;
				}
				
				if(!ligandsUnique.contains(ligandData[0])) {
					
					ligandsUnique.add(ligandData[0]);
				
					IAtomContainer ac = LigandTools.readSmilesStringandAddHydrogens(ligandData[1].trim(), true);
					
					psL.setString(1, ligandData[0]);
					psL.setString(2, UUID.nameUUIDFromBytes(ligandData[0].getBytes()).toString());
					psL.setString(3, ligandData[1]);
					psL.setString(4, descRot.calculate(ac).getValue().toString());
					psL.setString(5, descHD.calculate(ac).getValue().toString());
					psL.setString(6, descHA.calculate(ac).getValue().toString());
					psL.setString(7, String.valueOf(r(((DoubleResult) descW.calculate(ac).getValue()).doubleValue())));
					psL.setString(8, descXlogp.calculate(ac).getValue().toString());

					psL.addBatch();
				}
				
				if(!proteinLigandsUnique.contains(dockingRow[0] + "_" + dockingRow[2])) {
				
					proteinLigandsUnique.add(dockingRow[0] + "_" + dockingRow[2]);
					
					psPL.setInt(1, proteinLigandId);
					psPL.setString(2, dockingRow[0]);
					psPL.setString(3, ligandData[0]);
					psPL.setString(4, dockingRow[2]);
					
					psPL.addBatch();
					
					proteinLigandId++;
				}
				
				String uuidString = variantData[1]+ "_" + ligandData[0];
				psVL.setInt(1, variantLigandId);
				psVL.setString(2, UUID.nameUUIDFromBytes(uuidString.getBytes()).toString());
				psVL.setInt(3, Integer.parseInt(variantData[0]));
				psVL.setString(4, dockingRow[1]); // variant_folder
				psVL.setString(5, ligandData[0]); //chembl_id
				psVL.setString(6, dockingRow[2]); //ligand_folder
				psVL.setString(7, dockingRow[3]); //binding_affinity
				
				psVL.addBatch();
				
				variantLigandId++;

				if(variantLigandId % 10000 == 0) {
					psL.executeBatch();
					psPL.executeBatch();
					psVL.executeBatch();
					con.commit();
				}
				
				if(variantLigandId % 50000 == 0 || variantLigandId ==  variantLigandTotal) {
					System.out.println(Math.round(((variantLigandId / (double) variantLigandTotal) * 100)) + "% of dockings results data were inserted into DB");
				}
			}
			
			psL.executeBatch();
			psPL.executeBatch();
			psVL.executeBatch();
			con.commit();
			
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (CDKException e) {
			e.printStackTrace();
		}finally{
			try {
				psL.close();
				psPL.close();
				psVL.close();
				con.close();
			} catch (SQLException e) {
				e.printStackTrace();
			}
		} 
	}
	
	/**
     * Insert variants data into the database
     * 
     */
	public static void populateDBWithVariantsData() {
				
		Connection con = null;
		PreparedStatement ps = null;
		String query = "INSERT INTO psnpbind_variant values (?,?,?,?,?,?,?,?,?)";
		
		Set<String> pdbs = new HashSet<String>();

		CsvParserSettings settings = new CsvParserSettings();
		
		settings.getFormat().setLineSeparator("\n");
		settings.getFormat().setDelimiter('\t');
		settings.setNumberOfRowsToSkip(1);

		CsvParser parser = new CsvParser(settings);
		
		List<String[]> rows = parser.parseAll(new File(Config.getProperty("TSV_PATH")+"/pdbbind_pocket_variants.tsv"));
		
		int variantId = 1;
		int variantFolder = 1;
		
		String currentPdb = "";
		
		// 1: variant_id
		// 2: variant_uuid
		// 3: variant_folder
		// 4: pdb_id
		// 5: variant_type
		// 6: source_aa
		// 7: target_aa
		// 8: residue_num
		// 9: chain
		
		try {
			
			con = DatabaseConnection.getConnection();
			ps = con.prepareStatement(query);
			
			con.setAutoCommit(false);
		
			for(String[] rowArr: rows){
				
				if(!pdbs.contains(rowArr[4])){
					pdbs.add(rowArr[4]);
				}
			
				if(!currentPdb.equals(rowArr[4])){
					variantFolder = 1;
					currentPdb = rowArr[4];
				}
				
				ps.setInt(1,variantId);
				ps.setString(2,UUID.nameUUIDFromBytes(String.valueOf(variantId).getBytes()).toString());
				ps.setString(3,String.valueOf(variantFolder));
				ps.setString(4,rowArr[4]);
				ps.setString(5,rowArr[7]+rowArr[10]+rowArr[12]+rowArr[8]);
				ps.setString(6,rowArr[7]);
				ps.setString(7,rowArr[8]);
				ps.setString(8,rowArr[12]);
				ps.setString(9,rowArr[10]);
						
				ps.addBatch();
				variantFolder++;
				variantId++;
			}
			
			for(String pdb: pdbs){
						
				ps.setInt(1,variantId);
				ps.setString(2,UUID.nameUUIDFromBytes(String.valueOf(variantId).getBytes()).toString());
				ps.setString(3,"WT");
				ps.setString(4,pdb);
				ps.setString(5,"WT");
				ps.setString(6,"-");
				ps.setString(7,"-");
				ps.setString(8,"-");
				ps.setString(9,"-");

				ps.addBatch();
				variantId++;
			}
			
			ps.executeBatch();
			con.commit();
			
		} catch (SQLException e) {
			e.printStackTrace();
		}finally{
			try {
				ps.close();
				con.close();
			} catch (SQLException e) {
				e.printStackTrace();
			}
		} 
	}
	
	/**
     * Insert protein data into the database. In the process, an API call to RCSB protein data bank is made to obtain information
     * about each protein.
     * 
     */
	public static void populateDBWithProteinsData() {
		
		Connection con = null;
		PreparedStatement ps = null;
		String query = "INSERT INTO psnpbind_protein values (?,?,?,?,?,?,?,?,?,?,?,?)";
				
		String[] pdbs = {"1owh", "2c3i", "2hb1", "2pog", "2weg", "2y5h","3b27", 
				"3b5r", "3fv1","3jvr","3pxf", "3u9q", "3udh", "3up2", "3utu","4crc", 
				"4dli","4e5w", "4gr0","4j21","4jia","4m0y", "4wiv","4twp","5a7b", "5c28"};
		
		String[] uniprot = {"P00749", "P11309", "P18031", "P03372", "P00918", "P00742","P07900", 
				"P10275", "P39086","O14757","P24941", "P37231", "P56817", "O14965", "P00734","P03951", 
				"Q16539","P23458", "P39900","Q9H2K2","O60674","Q08881", "O60885","P00519","P04637", "Q9Y233"};
		
		JSONParser parser = new JSONParser();

        int index = 0;

        try {
        	
			con = DatabaseConnection.getConnection();
			ps = con.prepareStatement(query);
			
			con.setAutoCommit(false);
	
			for(String pdb: pdbs){
			
	            URL rcsb = new URL("https://data.rcsb.org/rest/v1/core/entry/"+pdb);
	            URLConnection yc = rcsb.openConnection();
	            BufferedReader in = new BufferedReader(new InputStreamReader(yc.getInputStream()));
	           	            
	            // 1: PDB ID
	            // 2: Protein UUID
	            // 3: Protein Name
	            // 4: UniProt ID
	            // 5: Exp method
	            // 6: resolution
	            // 7: deposit date
	            // 8: revised date
	            // 9: chains
	            // 10: sequence length
	            // 11: molecular weight
	            // 12: atom count

	            String inputLine;
	            
	            while ((inputLine = in.readLine()) != null) {              
	            	
	            	JSONObject root = (JSONObject) parser.parse(inputLine);
	               
	                JSONObject entry = (JSONObject) root.get("entry");
	                
	                String depositDate = ((JSONObject) root.get("rcsb_accession_info")).get("deposit_date").toString();
	                String revisionDate = ((JSONObject) root.get("rcsb_accession_info")).get("revision_date").toString();
	                
	                SimpleDateFormat inFormat, outFormat;
	                inFormat = new SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ss'+'SSS");
	                outFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");

	                Date date = inFormat.parse(depositDate);
	                depositDate = outFormat.format(date).toString();
	                
	                date = inFormat.parse(revisionDate);
	                revisionDate = outFormat.format(date).toString();

	                ps.setString(1,pdb);
	                ps.setString(2,UUID.nameUUIDFromBytes(entry.get("id").toString().toLowerCase().getBytes()).toString());
	                ps.setString(3,((JSONObject) root.get("struct")).get("pdbx_descriptor").toString());
	                ps.setString(4,uniprot[index]);
	                ps.setString(5,((JSONObject) root.get("rcsb_entry_info")).get("experimental_method").toString());
	                ps.setFloat(6,Float.parseFloat(((JSONObject) root.get("pdbx_vrpt_summary")).get("pdbresolution").toString()));
	                ps.setTimestamp(7,Timestamp.valueOf(depositDate));
	                ps.setTimestamp(8,Timestamp.valueOf(revisionDate));
	                ps.setInt(9,Integer.parseInt(((JSONObject) root.get("rcsb_entry_info")).get("deposited_polymer_entity_instance_count").toString()));
	                ps.setInt(10,Integer.parseInt(((JSONObject) root.get("rcsb_entry_info")).get("deposited_modeled_polymer_monomer_count").toString()));
	                ps.setFloat(11,Float.parseFloat(((JSONObject) root.get("rcsb_entry_info")).get("molecular_weight").toString()));
	                ps.setInt(12,Integer.parseInt(((JSONObject) root.get("rcsb_entry_info")).get("deposited_atom_count").toString()));
          
	                ps.addBatch();
	            }
		            
		        in.close();
				index++;

			} // for loop
			
			ps.executeBatch();
			con.commit();
			
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ParseException e) {
            e.printStackTrace();
        } catch (java.text.ParseException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		}finally{
			try {
				ps.close();
				con.close();
			} catch (SQLException e) {
				e.printStackTrace();
			}
		} 
	}
	
	/**
	 * A method to round a double number to four digits after the floating point
	 * @param value the value to be rounded
	 * @return a double value of the rounded number
	 */
	public static double r(double value) {

		return (double) Math.round(value * Math.pow(10, 4)) / Math.pow(10, 4);
	}
}
