package io.github.ammar257ammar.psnpbind.core;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

import io.github.ammar257ammar.psnpbind.core.model.PdbBindDataset;
import io.github.ammar257ammar.psnpbind.core.model.PdbBindDataset.PdbbindAttribute;

public class PdbBindDatasetTest {

	private PdbBindDataset pdbbindDataset;

	@Before
	public void setUp() throws Exception {

		pdbbindDataset = PdbBindDataset.create().loadData();
	}

	@Test
	public void getDataTest() {
		assertEquals(16139, pdbbindDataset.getData().size());
	}

	@Test
	public void keepAsFolderMatchTest() {
		assertEquals(285, pdbbindDataset.keepAsFolderMatch().getData().size());
	}

	@Test
	public void filterStringNotEqualUniprotTest() {
		assertEquals(15898, pdbbindDataset.filterStringNotEqual(PdbbindAttribute.UNIPROT, "------").getData().size());
	}

	@Test
	public void filterStringNotEqualNmrTest() {
		assertEquals(15862, pdbbindDataset.filterStringNotEqual(PdbbindAttribute.RESOLUTION, "NMR").getData().size());
	}

}
