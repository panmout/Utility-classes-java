package gr.uth.ece.dsel.aknn_hadoop;

import java.io.IOException;
import java.util.Formatter;
import java.util.FormatterClosedException;

public final class WriteLocalFiles
{	
	public static void writeFile(String file, String content)
	{
		try
		{
			Formatter outputTextFile = new Formatter(file);
			outputTextFile.format(content);
			outputTextFile.close();
		}
		catch (FormatterClosedException formatterException)
		{
			System.err.println("Error writing to file, exiting");
			System.exit(2);
		}
		catch (IOException ioException)
		{
			System.err.println("Error writing to file, exiting");
			System.exit(3);
		}
	}
}
