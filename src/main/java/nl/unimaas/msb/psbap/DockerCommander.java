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

package nl.unimaas.msb.psbap;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * 
 * @author Ammar Ammar
 * 
 * A basic class to call docker containers to perform specific tasks
 *
 */

public class DockerCommander {
	

	/**
	 * A method to dowload SIFTS files from a TSV file of URLs one in each line
	 * @param path a string of the path of the downloaded SIFTS files
	 * @param downloadUrlsPath a string of the path of the TSV file of SIFTS URLs
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public static void downloadSifts(String path, String downloadUrlsPath) throws IOException, InterruptedException {
				
		ProcessBuilder builder = new ProcessBuilder();
		
		builder.command("docker","run",
				"-v", path+":/sifts",
				"-v", "/home/apc/workspace/PockerSnpBindingAffinity/"+downloadUrlsPath+":/url",
				"sifts", "--quiet", "--no-clobber","--retry-connrefused","--waitretry=1","--read-timeout=60","--timeout=60","-t","0","-i","/url/pdbbind_sifts_urls.tsv");
		
		System.out.println("Process started!!");

		Process process = builder.start();
		
		BufferedReader error = new BufferedReader(new InputStreamReader(process.getErrorStream()));

	    BufferedReader output = new BufferedReader(new InputStreamReader(process.getInputStream()));

		StringBuilder sbuilder = new StringBuilder();
		
		String line = null;
		
		while ( (line = output.readLine()) != null) {
			System.out.println("in the output!!");
			sbuilder.append(line);
		   sbuilder.append(System.getProperty("line.separator"));
		}	    
		
		System.out.println("Wating for download to finish!!");
		
	    int exitVal = process.waitFor();
		if (exitVal == 0) {
			System.out.println("Success!");
			System.out.println(sbuilder.toString());
		} else {
			System.out.println(error.readLine());
		}
	    
		System.out.println("Download finished!!");

		process.destroy();
		if (process.isAlive()) {
		    process.destroyForcibly();
		}	
		
		System.out.println("Process destroyed!!");

	}

}
