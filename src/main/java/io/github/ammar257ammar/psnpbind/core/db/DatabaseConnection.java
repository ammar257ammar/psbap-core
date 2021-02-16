package io.github.ammar257ammar.psnpbind.core.db;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

import io.github.ammar257ammar.psnpbind.core.Config;


public class DatabaseConnection {

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
