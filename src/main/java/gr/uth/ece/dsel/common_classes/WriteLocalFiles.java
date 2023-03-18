package gr.uth.ece.dsel.common_classes;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Formatter;
import java.util.FormatterClosedException;

public final class WriteLocalFiles
{
    public static void writeFile(String file, String content, boolean appendable)
    {
        try
        {
            Formatter outputTextFile = new Formatter(new FileWriter(file, appendable));
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