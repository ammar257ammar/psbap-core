package io.github.ammar257ammar.psnpbind.core;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.biojava.nbio.aaproperties.profeat.convertor.Convert2Charge;
import org.biojava.nbio.aaproperties.profeat.convertor.Convert2Hydrophobicity;
import org.biojava.nbio.aaproperties.profeat.convertor.Convert2NormalizedVanDerWaalsVolume;
import org.biojava.nbio.aaproperties.profeat.convertor.Convert2Polarity;
import org.biojava.nbio.aaproperties.profeat.convertor.Convert2Polarizability;
import org.biojava.nbio.aaproperties.profeat.convertor.Convert2SecondaryStructure;
import org.biojava.nbio.aaproperties.profeat.convertor.Convert2SolventAccessibility;
import org.biojava.nbio.structure.StructureException;
import org.biojava.nbio.structure.secstruc.SecStrucState;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.config.Elements;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.geometry.surface.NumericalSurface;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FractionalPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.JPlogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.SmallRingDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WHIMDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IntegerArrayResult;

import com.univocity.parsers.tsv.TsvParser;
import com.univocity.parsers.tsv.TsvParserSettings;

import io.github.ammar257ammar.psnpbind.core.utils.AAprops;
import io.github.ammar257ammar.psnpbind.core.utils.LigandTools;
import io.github.ammar257ammar.psnpbind.core.utils.PdbTools;


public class Featurizer {

	public static double r(double value) {

		return (double) Math.round(value * Math.pow(10, 4)) / Math.pow(10, 4);
	}

	public static List<String[]> getSnpsFeatures(String path, String foldxPath, String singlePDB)
			throws IOException, StructureException {

		TsvParserSettings settings = new TsvParserSettings();
		settings.getFormat().setLineSeparator("\n");

		TsvParser parser = new TsvParser(settings);

		List<String[]> rows = parser.parseAll(new FileReader(path));
		rows.remove(0);

		List<String[]> annotatedSnps = new ArrayList<String[]>(rows.size());

		List<String> header = new ArrayList<String>();

		header.add("uniprot");
		header.add("snp");
		header.add("rs_id");
    header.add("mutation_type");
		header.add("pdb");
		header.add("min_res");
		header.add("ligand");
		header.add("sourceAminoAcid");
		header.add("targetAminoAcid");
		header.add("residueNum");
		header.add("chain");
		header.add("PdbResName");
		header.add("PdbResNum");
		header.add("aaPDBName");
		header.add("aaResidueNumber");
		header.add("UniProtResName");
		header.add("UniProtPos");
		header.add("UniProtAccessionId");
		header.add("PdbId");
		header.add("SeqResName");
		header.add("NaturalPos");

		header.add("FoldXname");
		header.add("FoldXmutation");

		header.add("secStruct");
		header.add("secStructSimple");

		header.add("CysteineMutation");
		header.add("GlycineMutation");
		header.add("ProlineMutation");

		header.add("ChargeGroupChange");
		header.add("HydroGroupChange");
		header.add("VanDerWaalsVolumeGroupChange");
		header.add("PloarityGroupChange");
		header.add("PolarizabilityGroupChange");
		header.add("SSWTGroupChange");
		header.add("AsaGroupChange");

		header.add("AsaChange");

		header.add("mutationPhi");
		header.add("mutationPsi");

		List<String> aaPropsHeader = AAprops.getAApropsHeader();

		for (int i = 0; i < aaPropsHeader.size(); i++) {
			header.add(aaPropsHeader.get(i));
		}

		for (int i = 0; i < aaPropsHeader.size(); i++) {
			header.add(aaPropsHeader.get(i) + "Surrounding");
		}

		Map<String, Integer> countMap = new HashMap<String, Integer>();

		boolean foldxHeaderAdded = false;

		for (String[] row : rows) {

			if (singlePDB.equals("all") || (!singlePDB.equals("all") && singlePDB.equals(row[4]))) {

				Integer count = countMap.get(row[4]);
				if (count == null)
					count = 0;

				countMap.put(row[4], count + 1);

				String pdb = row[4];

				if (!foldxHeaderAdded) {
					header = Featurizer.addFoldxHeader(foldxPath, header, pdb);
					foldxHeaderAdded = true;
				}

				List<String> wtFeatures = Featurizer.getMutantFeatures(row, countMap, foldxPath);
				List<String> mutatedFeatures = Featurizer.getWildTypeFeatures(row, countMap);

				annotatedSnps.add(wtFeatures.toArray(new String[wtFeatures.size()]));
				annotatedSnps.add(mutatedFeatures.toArray(new String[mutatedFeatures.size()]));

			} // if all PDBs
		}

		annotatedSnps.add(0, header.toArray(new String[header.size()]));

		return annotatedSnps;
	}

	public static List<String> getMutantFeatures(String[] row, Map<String, Integer> countMap, String foldxPath)
			throws IOException, StructureException {

		String pdb = row[4];

		List<String> annotatedSnp = new ArrayList<String>();
		Collections.addAll(annotatedSnp, row);

		System.out.println(pdb);

		annotatedSnp.add(pdb + "_protein_Repair_" + countMap.get(pdb));
		annotatedSnp.add(row[7] + row[10] + row[12] + row[8] + ";");

		List<SecStrucState> dssp = PdbTools.getDsspForPDB(Config.getProperty("DSSP_PATH") + pdb + ".dssp.gz", pdb);

		String ss = PdbTools.getSnpSecStruc(dssp, Integer.parseInt(row[12]));
		annotatedSnp.add(ss);

		String ssSimple = PdbTools.getSnpHelixOrStrand(dssp, Integer.parseInt(row[12]));
		annotatedSnp.add(ssSimple);

		annotatedSnp.add(row[7].equals("C") || row[8].equals("C") ? "YES" : "NO");
		annotatedSnp.add(row[7].equals("G") || row[8].equals("G") ? "YES" : "NO");
		annotatedSnp.add(row[7].equals("P") || row[8].equals("P") ? "YES" : "NO");

		Convert2Charge con1 = new Convert2Charge();
		Convert2Hydrophobicity con2 = new Convert2Hydrophobicity();
		Convert2NormalizedVanDerWaalsVolume con3 = new Convert2NormalizedVanDerWaalsVolume();
		Convert2Polarity con4 = new Convert2Polarity();
		Convert2Polarizability con5 = new Convert2Polarizability();
		Convert2SecondaryStructure con6 = new Convert2SecondaryStructure();
		Convert2SolventAccessibility con7 = new Convert2SolventAccessibility();

		// ------------ Get AA group change----------------

		annotatedSnp.add(String.valueOf("group" + con1.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con1.convert(row[8].trim().charAt(0))));

		annotatedSnp.add(String.valueOf("group" + con2.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con2.convert(row[8].trim().charAt(0))));

		annotatedSnp.add(String.valueOf("group" + con3.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con3.convert(row[8].trim().charAt(0))));

		annotatedSnp.add(String.valueOf("group" + con4.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con4.convert(row[8].trim().charAt(0))));

		annotatedSnp.add(String.valueOf("group" + con5.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con5.convert(row[8].trim().charAt(0))));

		annotatedSnp.add(String.valueOf("group" + con6.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con6.convert(row[8].trim().charAt(0))));

		annotatedSnp.add(String.valueOf("group" + con7.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con7.convert(row[8].trim().charAt(0))));

		// ------------ Get AA ASA----------------

		String wtProtein = Config.getProperty("VINA_DOCKING_DIR") + pdb + "/proteins/" + pdb + "_protein_Repair_WT/"
				+ pdb + "_protein_Repair_WT_final.pdb";

		String mutatedProtein = Config.getProperty("VINA_DOCKING_DIR") + pdb + "/proteins/" + pdb + "_protein_Repair_"
				+ countMap.get(pdb) + "/" + pdb + "_protein_Repair_" + countMap.get(pdb) + "_final.pdb";

		if (new File(wtProtein).exists()) {

			double asaWT = PdbTools.getResidueASA(wtProtein, row[12]);

			double asaMutation = PdbTools.getResidueASA(mutatedProtein, row[12]);

			if (asaWT == 0.0)
				asaWT = 0.0001;

			annotatedSnp.add(String.valueOf(r(r(asaMutation) / r(asaWT))));

			Pair<Double, Double> torsionAngles = PdbTools.getResiduePhiPsi(mutatedProtein, row[12]);

			annotatedSnp.add(String.valueOf(r(torsionAngles.getLeft())));
			annotatedSnp.add(String.valueOf(r(torsionAngles.getRight())));

		} else {
			annotatedSnp.add(String.valueOf(""));
			annotatedSnp.add(String.valueOf(""));

			annotatedSnp.add(String.valueOf(""));
			annotatedSnp.add(String.valueOf(""));
		}

		// -------- Get AAprops----------------

		List<Double> mutationProps = AAprops.getResidueAndSurroundingProps(mutatedProtein, row[12], row[7]);

		for (Double prop : mutationProps) {
			annotatedSnp.add(String.valueOf(r(prop)));
		}

		// -------- Get Foldx energy terms----------------

		String buildModelPath = foldxPath + "/" + pdb + "/output" + "/Average_" + pdb + "_protein_Repair.fxout";

		BufferedReader reader;
		String line;

		try {

			if (new File(buildModelPath).exists()) {
				reader = new BufferedReader(new FileReader(buildModelPath));
				line = "";

				while ((line = reader.readLine()) != null) {

					if (line.startsWith(pdb + "_protein_Repair_" + countMap.get(pdb))) {
						String[] lineArr = line.split("\t");

						for (int i = 2; i < 24; i++) {
							annotatedSnp.add(String.valueOf(r(Double.parseDouble(lineArr[i]))));
						}
						break;
					}
				}
				reader.close();
			}

		} catch (IOException e) {
			e.printStackTrace();
		}

		return annotatedSnp;
	}

	public static List<String> getWildTypeFeatures(String[] row, Map<String, Integer> countMap)
			throws IOException, StructureException {

		String pdb = row[4];

		List<String> annotatedSnp = new ArrayList<String>();
		Collections.addAll(annotatedSnp, row);

		System.out.println(pdb);

		annotatedSnp.add(pdb + "_protein_Repair_WT");
		annotatedSnp.add(row[7] + row[10] + row[12] + row[8] + ";");

		List<SecStrucState> dssp = PdbTools.getDsspForPDB(Config.getProperty("DSSP_PATH") + pdb + ".dssp.gz", pdb);

		String ss = PdbTools.getSnpSecStruc(dssp, Integer.parseInt(row[12]));
		annotatedSnp.add(ss);

		String ssSimple = PdbTools.getSnpHelixOrStrand(dssp, Integer.parseInt(row[12]));
		annotatedSnp.add(ssSimple);

		annotatedSnp.add(row[7].equals("C") || row[8].equals("C") ? "YES" : "NO");
		annotatedSnp.add(row[7].equals("G") || row[8].equals("G") ? "YES" : "NO");
		annotatedSnp.add(row[7].equals("P") || row[8].equals("P") ? "YES" : "NO");

		Convert2Charge con1 = new Convert2Charge();
		Convert2Hydrophobicity con2 = new Convert2Hydrophobicity();
		Convert2NormalizedVanDerWaalsVolume con3 = new Convert2NormalizedVanDerWaalsVolume();
		Convert2Polarity con4 = new Convert2Polarity();
		Convert2Polarizability con5 = new Convert2Polarizability();
		Convert2SecondaryStructure con6 = new Convert2SecondaryStructure();
		Convert2SolventAccessibility con7 = new Convert2SolventAccessibility();

		// ------------ Get AA group change----------------

		annotatedSnp.add(String.valueOf("group" + con1.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con1.convert(row[7].trim().charAt(0))));

		annotatedSnp.add(String.valueOf("group" + con2.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con2.convert(row[7].trim().charAt(0))));

		annotatedSnp.add(String.valueOf("group" + con3.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con3.convert(row[7].trim().charAt(0))));

		annotatedSnp.add(String.valueOf("group" + con4.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con4.convert(row[7].trim().charAt(0))));

		annotatedSnp.add(String.valueOf("group" + con5.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con5.convert(row[7].trim().charAt(0))));

		annotatedSnp.add(String.valueOf("group" + con6.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con6.convert(row[7].trim().charAt(0))));

		annotatedSnp.add(String.valueOf("group" + con7.convert(row[7].trim().charAt(0)) + "-" + "group"
				+ con7.convert(row[7].trim().charAt(0))));

		// ------------ Get AA ASA----------------

		String wtProtein = Config.getProperty("VINA_DOCKING_DIR") + pdb + "/proteins/" + pdb + "_protein_Repair_WT/"
				+ pdb + "_protein_Repair_WT_final.pdb";

		if (new File(wtProtein).exists()) {

			annotatedSnp.add("1.0");

			Pair<Double, Double> torsionAngles = PdbTools.getResiduePhiPsi(wtProtein, row[12]);

			annotatedSnp.add(String.valueOf(r(torsionAngles.getLeft())));
			annotatedSnp.add(String.valueOf(r(torsionAngles.getRight())));

		} else {
			annotatedSnp.add(String.valueOf(""));
			annotatedSnp.add(String.valueOf(""));

			annotatedSnp.add(String.valueOf(""));
			annotatedSnp.add(String.valueOf(""));
		}

		// -------- Get AAprops----------------

		List<Double> mutationProps = AAprops.getResidueAndSurroundingProps(wtProtein, row[12], row[7]);

		for (Double prop : mutationProps) {
			annotatedSnp.add(String.valueOf(r(prop)));
		}

		// -------- Get Foldx energy terms----------------

		for (int i = 0; i < 22; i++) {
			annotatedSnp.add("0.0");
		}

		return annotatedSnp;

	}

	public static List<String> addFoldxHeader(String foldxPath, List<String> header, String pdb) {

		String buildModelPath = foldxPath + "/" + pdb + "/output" + "/Average_" + pdb + "_protein_Repair.fxout";

		BufferedReader reader;
		String line;

		try {

			if (new File(buildModelPath).exists()) {

				reader = new BufferedReader(new FileReader(buildModelPath));
				line = "";

				while ((line = reader.readLine()) != null) {

					if (line.startsWith("Pdb") && line.split("\t")[0].trim().equals("Pdb")) {

						String[] lineArr = line.split("\t");

						for (int i = 2; i < 24; i++) {
							header.add(lineArr[i].replace(" ", "_").trim());
						}
					}
				}
				reader.close();
			}

		} catch (IOException e) {
			e.printStackTrace();
		}

		return header;
	}

	public static List<String[]> getLigandsFeatures(String ligandsPath, String filePath, String singlePDB)
			throws IOException, StructureException, ClassNotFoundException, CDKException, CloneNotSupportedException {

		TsvParserSettings settings = new TsvParserSettings();
		settings.getFormat().setLineSeparator("\n");

		TsvParser parser = new TsvParser(settings);

		List<String[]> rows = parser.parseAll(new FileReader(filePath));

		List<String[]> annotatedLigands = new ArrayList<String[]>(rows.size());

		List<String> header = new ArrayList<String>();

		header.add("pdb");
		header.add("ligand_file");
		header.add("chembl_id");
		header.add("tanimoto_index");

		header.add("WeightDescriptor");
		header.add("RuleOfFiveDescriptor");
		header.add("RotatableBondsCountDescriptor");
		header.add("HBondDonorCountDescriptor");
		header.add("HBondAcceptorCountDescriptor");
		header.add("APolDescriptor");

		for (int i = 1; i <= 9; i++) {
			header.add("SmallRingDescriptor" + i);
		}

		header.add("VABCDescriptor");
		header.add("TPSADescriptor");
		header.add("FractionalPSADescriptor");

		for (int i = 1; i <= 7; i++) {
			header.add("MomentOfInertiaDescriptor" + i);
		}

		header.add("AromaticAtomsCountDescriptor");
		header.add("AromaticBondsCountDescriptor");
		header.add("AtomCountDescriptor");
		header.add("XLogPDescriptor");
		header.add("VAdjMaDescriptor");
		header.add("BasicGroupCountDescriptor");
		header.add("TotalSurfaceArea");

		header.add("AlogPDescriptor1");
		header.add("AlogPDescriptor2");
		header.add("AlogPDescriptor3");

		header.add("JPlogPDescriptor");

		for (int i = 1; i <= 6; i++) {
			header.add("BCUTDescriptor" + i);
		}

		String[] weights = { "unity", "mass", "volume", "eneg", "polar" };

		for (String weight : weights) {
			for (int i = 1; i <= 17; i++) {
				header.add("WHIMDescriptor" + "_" + weight + i);
			}
		}

		header.add("AtomCountDescriptorN");
		header.add("AtomCountDescriptorO");
		header.add("Electronegativity");

		for (int i = 1; i <= 9; i++) {
			header.add("CarbonTypesDescriptor" + i);
		}

		for (int i = 1; i <= 1024; i++) {
			header.add("CircularFingerprinter" + i);
		}

		List<String> annotatedLigand = null;

		for (String[] row : rows) {

			if (singlePDB.equals("all") || (!singlePDB.equals("all") && singlePDB.equals(row[0]))) {

				annotatedLigand = new ArrayList<String>();
				Collections.addAll(annotatedLigand, row);

				File ligandFile = new File(
						Config.getProperty("LIGANDS_PATH") + row[0] + "/splitted/" + row[1] + "_min.mol2");

				File ligandSmilesFile = new File(
						Config.getProperty("LIGANDS_PATH") + row[0] + "/splitted-smi/" + row[1] + ".smi");

				if (ligandFile.exists() && ligandSmilesFile.exists()) {

					if (ligandFile.length() == 0) {

						for (int i = 4; i < header.size(); i++) {
							annotatedLigand.add("");
						}
						annotatedLigands.add(annotatedLigand.toArray(new String[annotatedLigand.size()]));

						continue;
					}

					IAtomContainer ac = LigandTools.readSmilesFileandAddHydrogens(ligandSmilesFile, true);
					IAtomContainer acCleaned = LigandTools.readMol2andAddHydrogens(ligandFile, false);

					if (ac == null || acCleaned == null) {

						System.out.println(row[0] + " - " + row[1] + " ERROR reading file!!");

						for (int i = 4; i < header.size(); i++) {
							annotatedLigand.add("");
						}
						annotatedLigands.add(annotatedLigand.toArray(new String[annotatedLigand.size()]));

						continue;
					}

					annotatedLigand = Featurizer.getLigandFeatures(ligandFile, ligandSmilesFile, annotatedLigand);

				} else {

					for (int i = 4; i < header.size(); i++) {
						annotatedLigand.add("");
					}
					annotatedLigands.add(annotatedLigand.toArray(new String[annotatedLigand.size()]));

					continue;
				}

				if (annotatedLigand != null)
					annotatedLigands.add(annotatedLigand.toArray(new String[annotatedLigand.size()]));

			} // if all pdbs
		}

		annotatedLigands.add(0, header.toArray(new String[header.size()]));

		return annotatedLigands;
	}

	public static List<String> getLigandFeatures(File ligandFile, File ligandSmilesFile, List<String> annotatedLigand)
			throws CDKException, ClassNotFoundException, IOException {

		// RULES OF FIVE
		// Molecular Weight, it is affected with adding hydrogens
		WeightDescriptor descW = new WeightDescriptor();
		RuleOfFiveDescriptor descRof = new RuleOfFiveDescriptor(); // not affected
		RotatableBondsCountDescriptor descRot = new RotatableBondsCountDescriptor();// not affected
		HBondDonorCountDescriptor descHD = new HBondDonorCountDescriptor();// not affected
		HBondAcceptorCountDescriptor descHA = new HBondAcceptorCountDescriptor();// not affected

		// Sum of atomic polarizability (including implcit hydrogens) so it is affected
		// with adding hydrogens
		APolDescriptor descApol = new APolDescriptor();

		/*
		 * Small ring count descriptor
		 * 
		 * IntegerArrayResult
		 * 
		 * result.add(nSmallRings); result.add(nAromRings); result.add(nRingBlocks);
		 * result.add(nAromBlocks); result.add(nRings3); result.add(nRings4);
		 * result.add(nRings5); result.add(nRings6); result.add(nRings7);
		 * result.add(nRings8); result.add(nRings9);
		 */
		SmallRingDescriptor descSR = new SmallRingDescriptor();

		// Volume descriptor
		VABCDescriptor descVol = new VABCDescriptor(); // not affected

		// total polar surface area
		TPSADescriptor descTpsa = new TPSADescriptor(); // not affected

		// Polar surface adjusted by molecular weight
		FractionalPSADescriptor descFpsa = new FractionalPSADescriptor(); // not affected

		// Gyration radius descriptor, needs 3D
		MomentOfInertiaDescriptor descRoJ = new MomentOfInertiaDescriptor();

		// Atom counts descriptor
		AromaticAtomsCountDescriptor descArCount = new AromaticAtomsCountDescriptor(); // not affected
		AromaticBondsCountDescriptor descArBoCount = new AromaticBondsCountDescriptor(); // not affected
		AtomCountDescriptor descAtCount = new AtomCountDescriptor();

		// XlogP
		XLogPDescriptor descXlogp = new XLogPDescriptor();
		ALOGPDescriptor descAlogp = new ALOGPDescriptor();
		JPlogPDescriptor descJPlogp = new JPlogPDescriptor();

		// Vertex adjacency information
		VAdjMaDescriptor descVadj = new VAdjMaDescriptor(); // not affected

		// Basic group count
		BasicGroupCountDescriptor descBG = new BasicGroupCountDescriptor();// not affected
		descBG.initialise(DefaultChemObjectBuilder.getInstance());// not affected

		BCUTDescriptor descBcut = new BCUTDescriptor();
		CarbonTypesDescriptor descCTD = new CarbonTypesDescriptor();

		WHIMDescriptor descWhim = new WHIMDescriptor();

		// MCFP fingerprint
		CircularFingerprinter circularFingerprinter = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4,
				1024); // not
						// affected

		// circularFingerprinter.setPerceiveStereo(true);

		// -------------------------------------------------------

		if (ligandFile.exists() && ligandSmilesFile.exists()) {

			if (ligandFile.length() == 0 || ligandSmilesFile.length() == 0)
				return null;

			IAtomContainer ac = LigandTools.readSmilesFileandAddHydrogens(ligandSmilesFile, true);

			// read mol2 without adding hydrogen, boolean parameter is false
			IAtomContainer acCleaned = LigandTools.readMol2andAddHydrogens(ligandFile, false);

			if (ac == null || acCleaned == null)
				return null;

			annotatedLigand.add(String.valueOf(r(((DoubleResult) descW.calculate(ac).getValue()).doubleValue())));
			annotatedLigand.add("V" + descRof.calculate(ac).getValue().toString());
			annotatedLigand.add(descRot.calculate(ac).getValue().toString());
			annotatedLigand.add(descHD.calculate(ac).getValue().toString());
			annotatedLigand.add(descHA.calculate(ac).getValue().toString());
			annotatedLigand.add(String.valueOf(r(((DoubleResult) descApol.calculate(ac).getValue()).doubleValue())));

			IntegerArrayResult iarSR = (IntegerArrayResult) descSR.calculate(ac).getValue();

			for (int i = 0; i < 9; i++) {
				annotatedLigand.add(String.valueOf(iarSR.get(i)));
			}

			annotatedLigand.add(String.valueOf(r(((DoubleResult) descVol.calculate(ac).getValue()).doubleValue())));
			annotatedLigand.add(descTpsa.calculate(ac).getValue().toString());
			annotatedLigand.add(String.valueOf(r(((DoubleResult) descFpsa.calculate(ac).getValue()).doubleValue())));

			DoubleArrayResult darRoJ = (DoubleArrayResult) descRoJ.calculate(acCleaned).getValue();

			for (int i = 0; i < 7; i++) {
				annotatedLigand.add(String.valueOf(r(darRoJ.get(i))));
			}

			annotatedLigand.add(descArCount.calculate(ac).getValue().toString());
			annotatedLigand.add(descArBoCount.calculate(ac).getValue().toString());
			descAtCount.setParameters(new Object[] { "*" });
			annotatedLigand.add(descAtCount.calculate(ac).getValue().toString());
			annotatedLigand.add(descXlogp.calculate(ac).getValue().toString());
			annotatedLigand.add(String.valueOf(r(((DoubleResult) descVadj.calculate(ac).getValue()).doubleValue())));
			annotatedLigand.add(descBG.calculate(ac).getValue().toString());
			annotatedLigand.add(String.valueOf(r(new NumericalSurface(acCleaned).getTotalSurfaceArea())));

			DoubleArrayResult darAlogp = (DoubleArrayResult) descAlogp.calculate(ac).getValue();

			annotatedLigand.add(String.valueOf(r(darAlogp.get(0))));
			annotatedLigand.add(String.valueOf(r(darAlogp.get(1))));
			annotatedLigand.add(String.valueOf(r(darAlogp.get(2))));

			annotatedLigand.add(String.valueOf(r(((DoubleResult) descJPlogp.calculate(ac).getValue()).doubleValue())));

			try {

				DoubleArrayResult darBcut = (DoubleArrayResult) descBcut.calculate(ac).getValue();

				for (int i = 0; i < 6; i++) {
					annotatedLigand.add(String.valueOf(r(darBcut.get(i))));
				}
				
			}catch(ArrayIndexOutOfBoundsException ex) {
				
				for (int i = 0; i < 6; i++) {
					annotatedLigand.add(String.valueOf(r(Double.NaN)));
				}
			}

			String[] weights = { "unity", "mass", "volume", "eneg", "polar" };

			DoubleArrayResult darWhim = null;

			for (String weight : weights) {

				descWhim.setParameters(new Object[] { weight });

				darWhim = (DoubleArrayResult) descWhim.calculate(acCleaned).getValue();

				for (int i = 0; i < 17; i++) {
					annotatedLigand.add(String.valueOf(r(darWhim.get(i))));
				}
			}

			descAtCount.setParameters(new Object[] { "N" });
			annotatedLigand.add(descAtCount.calculate(ac).getValue().toString());

			descAtCount.setParameters(new Object[] { "O" });
			annotatedLigand.add(descAtCount.calculate(ac).getValue().toString());

			double totalE = 0.0;

			for (IAtom atA : ac.atoms()) {
				Elements e = Elements.ofString(atA.getSymbol());
				totalE += e.electronegativity();
			}

			annotatedLigand.add(String.valueOf(totalE));

			IntegerArrayResult iarCTD = (IntegerArrayResult) descCTD.calculate(ac).getValue();

			for (int i = 0; i < 9; i++) {
				annotatedLigand.add(String.valueOf(iarCTD.get(i)));
			}

			int[] fparr = circularFingerprinter.getBitFingerprint(ac).getSetbits();

			int[] fp = new int[1024];

			for (int i = 0; i < fparr.length; i++) {
				fp[fparr[i]] = 1;
			}

			for (int i = 0; i < 1024; i++) {
				annotatedLigand.add(String.valueOf(fp[i]));
			}

		} else {
			return null;
		}

		return annotatedLigand;

	}

	public static List<String[]> getPocketsFeatures(String singlePDB) throws IOException, StructureException {
		
		File casf = new File(Config.getProperty("FOLDX_PDB_DIR"));
		File[] mols = casf.listFiles();
		
		List<String[]> annotatedPockets = new ArrayList<String[]>(mols.length);

		List<String> header = new ArrayList<String>();

		header.add("pdb");

		header.add("HelixSS");
		header.add("StrandSS");
		header.add("OtherSS");
		header.add("DominantSS");
		header.add("BuriedASA");
		header.add("ExposedASA");
		header.add("RatioASA");
		header.add("PocketASA");
		
		for(File molFolder: mols) {
			if(molFolder.isDirectory()) {
				
				if (singlePDB.equals("all") || (!singlePDB.equals("all") && singlePDB.equals(molFolder.getName()))) {
				
					List<String> annotatedPocket = Featurizer.getPocketFeatures(molFolder.getName());
					
					annotatedPockets.add(annotatedPocket.toArray(new String[annotatedPocket.size()]));
				}
			}
		}
		
		annotatedPockets.add(0, header.toArray(new String[header.size()]));

		return annotatedPockets;
	}

	public static List<String> getPocketFeatures(String pdb) throws IOException, StructureException {
		
		List<String> annotatedPocket = new ArrayList<String>();

		System.out.println("Pocket "+pdb);
		
		annotatedPocket.add(pdb);

		List<SecStrucState> dssp = PdbTools.getDsspForPDB(Config.getProperty("DSSP_PATH") + pdb + ".dssp.gz", pdb);

		Map<String, Double> frequencyMapSS = PdbTools.getPocketHelixStrandPercentageFromDSSP(dssp, pdb);
		
		annotatedPocket.add(String.valueOf(r(frequencyMapSS.get("Helix"))));
		annotatedPocket.add(String.valueOf(r(frequencyMapSS.get("Strand"))));
		annotatedPocket.add(String.valueOf(r(frequencyMapSS.get("Other"))));
		
		if(frequencyMapSS.get("Dominant") == 1000.0) {
			
			annotatedPocket.add("Helix");
			
		}else if(frequencyMapSS.get("Dominant") == 2000.0) {

			annotatedPocket.add("Strand");
			
		}else if(frequencyMapSS.get("Dominant") == 3000.0) {

			annotatedPocket.add("Other");

		}
		
    	Map<String, Double> frequencyMapASA = PdbTools.getPocketBuriedExposedASA(pdb);

		annotatedPocket.add(String.valueOf(r(frequencyMapASA.get("Buried"))));
		annotatedPocket.add(String.valueOf(r(frequencyMapASA.get("Exposed"))));
		annotatedPocket.add(String.valueOf(r(frequencyMapASA.get("Ratio"))));
		annotatedPocket.add(String.valueOf(r(frequencyMapASA.get("PocketASA"))));

    	
		return annotatedPocket;
	}

}