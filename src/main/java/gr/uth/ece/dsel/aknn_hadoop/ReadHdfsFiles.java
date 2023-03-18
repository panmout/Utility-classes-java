package gr.uth.ece.dsel.aknn_hadoop;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.util.HashMap;

import gr.uth.ece.dsel.common_classes.Node;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;

public final class ReadHdfsFiles
{
	// read treefile from hdfs
	public static Node getTree(String treeFile, FileSystem fs)
	{
		Node root = null;
		
		try
		{
			Path pt = new Path(treeFile); // create path object from path string
			ObjectInputStream input = new ObjectInputStream(fs.open(pt)); // open HDFS tree file
			root = (Node) input.readObject(); // assign quad tree binary form to root node
		}
		catch (ClassNotFoundException classNotFoundException)
		{
			System.err.println("Invalid object type");
		}
		catch (IOException e)
		{
			System.err.println("hdfs file does not exist");
		}
		return root;
	}
	
	// read mapreduce1 output from hdfs as hashmap
	public static HashMap<String, Integer> getMR1output(String mr1OutFull, FileSystem fs)
	{
		HashMap<String, Integer> cell_tpoints = new HashMap<>();
		
		try // open files
		{
			FileStatus[] status = fs.listStatus(new Path(mr1OutFull)); // FileStatus iterates through contents of hdfs dir
			
			BufferedReader reader;
			
			String line;
			
			// read mapreduce1 output from hdfs
			for (int i = 1; i < status.length; i++) // skipping status[0] = "_SUCCESS" file
			{
				reader = new BufferedReader(new InputStreamReader(fs.open(status[i].getPath()))); // create reader object from java data stream object
				
				while ((line = reader.readLine())!= null) // while input has more lines
				{
					String[] data = line.trim().split("\t");
					String cell = data[0]; // 1st element is point cell
					Integer num = Integer.parseInt(data[1]); // 2nd element is number of training points in cell
					cell_tpoints.put(cell, num); // add to hashmap
				}
				reader.close(); // close file
			}
		}
		catch (IOException e)
		{
			System.err.println("hdfs file does not exist");
		}
		
		return cell_tpoints;
	}
}
