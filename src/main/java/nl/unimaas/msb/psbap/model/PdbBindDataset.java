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

package nl.unimaas.msb.psbap.model;

import java.util.ArrayList;
import java.util.List;

import nl.unimaas.msb.psbap.Config;

public class PdbBindDataset {
	
	private String pathGeneralFile;
	private String pathGeneralNamesFile;
	private String entriesPath;
	
	private List<String[]> pdbbindData = new ArrayList<String[]>();
	
	private PdbBindDataset() {
		this.pathGeneralFile = Config.getProperty("PDBBIND_DATA_PATH_1");
		this.pathGeneralNamesFile = Config.getProperty("PDBBIND_DATA_PATH_2");
		this.entriesPath = Config.getProperty("PDBBIND_ENTRIES_PATH");
	}
	
	
	private PdbBindDataset(String pathGeneralFile, String pathGeneralNamesFile, String entriesPath) {
	    this.pathGeneralFile = pathGeneralFile;
	    this.pathGeneralNamesFile = pathGeneralNamesFile;
	    this.entriesPath = entriesPath;
	}
	
	

}
