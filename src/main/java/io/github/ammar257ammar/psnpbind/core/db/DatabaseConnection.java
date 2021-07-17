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

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

import io.github.ammar257ammar.psnpbind.core.Config;

/**
 * A utility class to create a connection to MySQL database
 * 
 * @author Ammar Ammar
 * 
 */
public class DatabaseConnection {

    /**
     * Create a connectin to PSnpBind MySQL database and return the connection to be used by other methods.
     * @return a database connection object
     * 
     */
	public static Connection getConnection() {

		Connection con = null;
		
		try {

			//Class.forName("com.mysql.cj.jdbc.Driver");

			con = DriverManager.getConnection("jdbc:mysql://"+Config.getProperty("DB_HOST") + ":" + 
					Config.getProperty("DB_PORT") + 
					"/"+Config.getProperty("DB_DATABASE"),
					Config.getProperty("DB_USERNAME"),
					Config.getProperty("DB_PASSWORD"));
			
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return con;
	}
}
