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

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.Mol2Reader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

/**
 * A utility class to manipulate ligand files mainly using CDK (Mol2 and SMILES formats)
 * 
 * @author Ammar Ammar
 *
 */
public class LigandTools {

	/**
	 * A method to read a molecule from mol2 file and parse it
	 * as a CDK IAtomContainer
	 * @param file the molecule file object
	 * @param addHydrogens boolean if hydrogen should be added to the molecule IAtomContainer
	 * @return a CDK IAtomContainer
	 */
	public static IAtomContainer readMol2andAddHydrogens(File file, boolean addHydrogens)
			throws IOException, ClassNotFoundException, CDKException {

		InputStream ins = new FileInputStream(file);

		Mol2Reader reader = new Mol2Reader(ins);

		IAtomContainer ac = null;

		try {
			ac = reader.read(DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainer.class));
		} catch (CDKException ex) {
			System.out.println("ERROR: reading file failed!!");
			reader.close();
			return null;
		}

		reader.close();

		for (IAtom atA : ac.atoms()) {

			if (atA.getBondCount() == 0) {
				System.out.println("Structure problem:" + file.getName());
				return null;
			}
		}

		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);

		for (IAtom atom : ac.atoms())
			atom.setImplicitHydrogenCount(0);

		if (addHydrogens) {

			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
			adder.addImplicitHydrogens(ac);

			AtomContainerManipulator.convertImplicitToExplicitHydrogens(ac);
		}

		try {

			Kekulization.kekulize(ac);

		} catch (CDKException ex) {
			System.out.println("kekulize failed!!");
		}

		return ac;
	}
	

	/**
	 * A method to read a molecule from mol2 String path and parse it
	 * as a CDK IAtomContainer
	 * @param file a String for the molecule file path
	 * @param addHydrogens boolean if hydrogen should be added to the molecule IAtomContainer
	 * @return a CDK IAtomContainer
	 */
	public static IAtomContainer readMol2andAddHydrogens(String file, boolean addHydrogens) {
		
		IAtomContainer ac = null;
		try {
			ac = readMol2andAddHydrogens(new File(file), addHydrogens);
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (CDKException e) {
			e.printStackTrace();
		}
		
		return ac;
	}
	
	/**
	 * A method to read a molecule from SMILES string file and parse it
	 * as a CDK IAtomContainer
	 * @param file the molecule SMILES file object
	 * @param addHydrogens boolean if hydrogen should be added to the molecule IAtomContainer
	 * @return a CDK IAtomContainer
	 */
	public static IAtomContainer readSMILESandAddHydrogens(File file, boolean addHydrogens)
			throws IOException, ClassNotFoundException, CDKException {

		List<String> lines = Files.readAllLines(Paths.get(file.getAbsolutePath()), StandardCharsets.UTF_8);

		String smile = lines.get(0).substring(0, lines.get(0).indexOf("CHEMBL")).trim();

		System.out.println(smile);
		
		SmilesParser sp = new SmilesParser(SilentChemObjectBuilder.getInstance());

		IAtomContainer ac;
		try {
			ac = sp.parseSmiles(smile);
		} catch (InvalidSmilesException e) {
			System.out.println("ERROR: reading file failed!!");
			return null;
		}

		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(ac);

		for (IAtom atom : ac.atoms())
			atom.setImplicitHydrogenCount(0);

		if (addHydrogens) {

			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(SilentChemObjectBuilder.getInstance());
			adder.addImplicitHydrogens(ac);

			AtomContainerManipulator.convertImplicitToExplicitHydrogens(ac);
		}

		ElectronDonation model = ElectronDonation.daylight();
		CycleFinder cycles = Cycles.or(Cycles.all(), Cycles.all(6));
		Aromaticity aromaticity = new Aromaticity(model, cycles);

		aromaticity.apply(ac);

		try {

			Kekulization.kekulize(ac);

		} catch (CDKException ex) {
			System.out.println("kekulize failed!!");
		}

		return ac;
	}

}
