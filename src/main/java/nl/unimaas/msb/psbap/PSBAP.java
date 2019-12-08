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

import java.util.List;

import nl.unimaas.msb.psbap.model.PdbBindDataset;
import nl.unimaas.msb.psbap.utils.DataHandler;

/**
 * Hello world!
 *
 */
public class PSBAP 
{
    public static void main( String[] args )
    {

    	CliOptions cli = new CliOptions(args);
    	
    	switch(cli.operation){
    	
    	case "print-pdbbind-head":
    		List<String[]> dataset = PdbBindDataset.create().loadData().getData();
			DataHandler.printDatasetHead(dataset);
    	}
    	
    	
    }
}
