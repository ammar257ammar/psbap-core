package nl.unimaas.msb.psbap;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

public class Config {
	
	private static Properties properties = new Properties();
	private static InputStream input = null;

	static {

		try {

			input = new FileInputStream("config.properties");

			properties.load(input);

		} catch (IOException ex) {

			ex.printStackTrace();

		} finally {

			if (input != null) {
				try {
					input.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}


}
