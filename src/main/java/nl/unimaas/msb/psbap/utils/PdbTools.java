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

import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.FileParsingParameters;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.io.LocalPDBDirectory.FetchBehavior;

/**
 * A utility class to perform operations on PDB and SIFTS files related to PdbBind Datasets
 *  
 * @author Ammar Ammar
 * 
 */

public class PdbTools {
	

	/**
	 * A method to configure a PDB file reader 
	 * 
	 * @param alignSeqAndPraseSS boolean value to specify if secondary structure info need to be parsed
	 *        and if SEQRES sequence and ATOM sequence need to be aligned when reading a PDB file
	 * @return a PDB file reader
	 */
	public static PDBFileReader configureReader(boolean alignSeqAndPraseSS){
		
		AtomCache cache = new AtomCache();

    	FileParsingParameters params = cache.getFileParsingParams();
    	
    	params.setAlignSeqRes(alignSeqAndPraseSS);
    	params.setParseSecStruc(alignSeqAndPraseSS);
    	
    	cache.setFileParsingParams(params);
    	
    	PDBFileReader reader = new PDBFileReader();
    	
    	reader.setFetchBehavior(FetchBehavior.LOCAL_ONLY);
    	
    	reader.setFileParsingParameters(params);
    	
    	return reader;
	}
	
	
	
	

}
