package io.github.ammar257ammar.psnpbind.core;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

import io.github.ammar257ammar.psnpbind.core.model.PDBbindEntry;

public class PDBbindEntryTest {

	private PDBbindEntry pdbbindEntry;

	@Before
	public void setUp() throws Exception {
		pdbbindEntry = new PDBbindEntry("2hb1", false, false);
	}

	@Test
	public void hasProteinStructureTest() {
		assertEquals(true, pdbbindEntry.hasProteinStructure());
	}

	@Test
	public void hasPocketStructureTest() {
		assertEquals(true, pdbbindEntry.hasPocketStructure());
	}

	@Test
	public void hasSiftsMappingTest() {
		assertEquals(true, pdbbindEntry.hasSiftsMapping());
	}

	@Test
	public void hasLigandStructureTest() {
		assertNull(pdbbindEntry.getLigandStructure());
	}

	@Test
	public void getSiftResiduesTest() {
		assertEquals(299, pdbbindEntry.getSiftResidues().size());
	}

	@Test
	public void getPocketAminoAcidsTest() {
		assertEquals(39, pdbbindEntry.getPocketAminoAcids().size());
	}

}
