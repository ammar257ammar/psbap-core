package io.github.ammar257ammar.psnpbind.core;

import static org.junit.Assert.*;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Before;
import org.junit.Test;

public class FoldXTest {

	private List<String[]> pdbbindPocketVariants;

	@Before
	public void setUp() throws Exception {

		pdbbindPocketVariants = new ArrayList<String[]>();
		pdbbindPocketVariants.add(new String[] { "O14757", "p.Leu92Phe", "rs775964760", "missense variant", "3jvr", "1.76",
				"3jvr:AGX;3jvs:AGY;1nvq:UCN;2br1:PFP;2brb:PFQ;", "L", "F", "92", "A", "L", "92", "L", "92", "L", "92", "O14757",
				"3jvr", "LEU", "91" });
	}

	@Test
	public void createMutationConfigFilesTest() throws IOException {

		FoldX.createMutationConfigFiles(pdbbindPocketVariants, Config.getProperty("PDBBIND_ENTRIES_PATH"),
				Config.getProperty("FOLDX_PDB_DIR"));

		File file = new File(Config.getProperty("FOLDX_PDB_DIR") + "3jvr/input/individual_list.txt");

		assertTrue(file.exists());

	}

}
