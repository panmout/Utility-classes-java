/*
package gr.uth.ece.dsel.common_classes;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Formatter;
import java.util.FormatterClosedException;
import java.util.HashSet;
import java.util.Random;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;

public class CreateQTreeArray
{
	private static Formatter outputTextFile; // local output tree text file
	private static FileOutputStream treeFileout;
	private static ObjectOutputStream outputTreeFile; // local output tree object file
    private static Point[] sampleArray;
	private static int capacity;
	private static int numCells = 0;
	private static final StringBuilder x = new StringBuilder();
	private static String treeFilePath;
	private static String treeFileName;
	private static String trainingDatasetPath;
	private static int samplerate;

    public CreateQTreeArray(int newCapacity, String newTreeFilePath, String newTreeFileName, String newArrayFilePath, String newArrayFileName, String newTrainingDatasetPath, int newSamplerate)
	{
		capacity = newCapacity;
		treeFilePath = newTreeFilePath;
		treeFileName = newTreeFileName;
		trainingDatasetPath = newTrainingDatasetPath;
		samplerate = newSamplerate;
	}
	
	private static void readSample()
	{
		try // open files
		{
			treeFileout = new FileOutputStream(treeFileName);
			outputTreeFile = new ObjectOutputStream(treeFileout); // open local output tree object file
			outputTextFile = new Formatter("qtree.txt"); // open local output text file
			
			FileSystem fs = FileSystem.get(new Configuration());
			Path trainingPath = new Path(trainingDatasetPath);
			BufferedReader trainingBr = new BufferedReader(new InputStreamReader(fs.open(trainingPath))); // open HDFS training dataset file

            ArrayList<Point> sampleArrayList = new ArrayList<Point>();
			
			HashSet<Integer> randomNumbers = new HashSet<Integer>(samplerate); // [percentSample] size set for random integers
			
			Random random = new Random();
			
			while (randomNumbers.size() < samplerate) // fill list
				randomNumbers.add(random.nextInt(100)); // add a random integer 0 - 99
						
			String line;
			// read training dataset and get sample points
			while ((line = trainingBr.readLine()) != null)
			{
				if (randomNumbers.contains(random.nextInt(100))) // [percentSample]% probability
				{
					String[] data = line.trim().split("\t");
					int pid = Integer.parseInt(data[0]); // tpoint id
					double x = Double.parseDouble(data[1]); // get x
					double y = Double.parseDouble(data[2]); // get y
					sampleArrayList.add(new Point(pid, x, y)); // add {pid, x, y} to ArrayList
				}
			}
			
			// convert ArrayList to Point Array
			sampleArray = sampleArrayList.toArray(new Point[sampleArrayList.size()]);
		}
		catch (IOException ioException)
		{
			System.err.println("Could not open file, exiting");
			System.err.println(ioException);
			System.exit(1);
		}
	}
	
	private static void writeFiles(Node node)
	{		
		// write to files
		try
		{
			// local
			outputTextFile = new Formatter("qtree.txt");
			outputTextFile.format("%s", x);
			outputTreeFile.writeObject(node);
			
			outputTreeFile.close();
			outputTextFile.close();
			treeFileout.close();
			
			// write to hdfs
			FileSystem fs = FileSystem.get(new Configuration());
			
			Path path = new Path(treeFilePath);
			ObjectOutputStream outputStream = new ObjectOutputStream(fs.create(path));
			outputStream.writeObject(node);
			outputStream.close();
			
			
		}
		catch (FormatterClosedException formatterException)
		{
			System.err.println("Error writing to file, exiting");
			System.err.println(formatterException);
			System.exit(2);
		}
		catch (IOException ioException)
		{
			System.err.println("Error writing to file, exiting");
			System.err.println(ioException);
			System.exit(3);
		}
	}
	
	public static void makeQuadArray(Point[] samplePoints, Node node)
	{
		double xmin = node.getXmin();
		double xmax = node.getXmax();
		double ymin = node.getYmin();
		double ymax = node.getYmax();
		
		int i1 = node.getLow();
		int i2 = node.getHigh();
		
		Arrays.sort(samplePoints, i1, i2 + 1, new PointXYComparator("min", 'x')); // sort subarray [i1, i2] by x-ascending
		
		int j2 = xBinarySearchArray(samplePoints, i1, i2, (xmin + xmax) / 2); // get index of middle x
		
		Arrays.sort(samplePoints, i1, j2, new PointXYComparator("min", 'y')); // sort subarray [i1, j2] by y-ascending
		
		int j1 = yBinarySearchArray(samplePoints, i1, j2, (ymin + ymax) / 2); // get index of middle y

		Arrays.sort(samplePoints, j2, i2 + 1, new PointXYComparator("min", 'y')); // sort subarray [j2, i2] by y-ascending
		
		int j3 = yBinarySearchArray(samplePoints, j2, i2, (ymin + ymax) / 2); // get index of middle y
		
		// check capacities and create children
		if (i2 - i1 + 1 > capacity) // root node capacity
		{
			node.setSW(new Node(xmin, ymin, (xmin + xmax) / 2, (ymin + ymax) / 2));
			node.setNW(new Node(xmin, (ymin + ymax) / 2, (xmin + xmax) / 2, ymax));
			node.setSE(new Node((xmin + xmax) / 2, ymin, xmax, (ymin + ymax) / 2));
			node.setNE(new Node((xmin + xmax) / 2, (ymin + ymax) / 2, xmax, ymax));
			
			node.getSW().setLow(i1);
			node.getSW().setHigh(j1 - 1);
			
			node.getNW().setLow(j1);
			node.getNW().setHigh(j2 - 1);
			
			node.getSE().setLow(j2);
			node.getSE().setHigh(j3 - 1);
			
			node.getNE().setLow(j3);
			node.getNE().setHigh(i2);
		}
		
		// children capacities
		if (j1 - i1 > capacity)
			makeQuadArray(samplePoints, node.getSW());
		
		if (j2 - j1 > capacity)
			makeQuadArray(samplePoints, node.getNW());
		
		if (j3 - j2 > capacity)
			makeQuadArray(samplePoints, node.getSE());
		
		if (i2 - j3 + 1 > capacity)
			makeQuadArray(samplePoints, node.getNE());
	}
	
	public static int xBinarySearchArray(Point[] array, int low, int high, double key)
	{
		int l = low;
		int r = high;
		
		int mid = -1;
		
		while (l <= r)
		{
			mid = (l + r) / 2;
			
			if (key > array[mid].getX())
				l = mid + 1;
			else if (key < array[mid].getX())
				r = mid - 1;
			else
				return mid;
		}
		
		return mid;
	}
	
	public static int yBinarySearchArray(Point[] array, int low, int high, double key)
	{
		int l = low;
		int r = high;
		
		int mid = -1;
		
		while (l <= r)
		{
			mid = (l + r) / 2;
			
			if (key > array[mid].getY())
				l = mid + 1;
			else if (key < array[mid].getY())
				r = mid - 1;
			else
				return mid;
		}
		
		return mid;
	}
	
	// create quadtree in string form
	private static void df_repr(Node node)
	{
		if (node.getNW() == null)
		{
			x.append("0");
			numCells++;
		}
		else
		{
			x.append("1");
			df_repr(node.getNW());
			df_repr(node.getNE());
			df_repr(node.getSW());
			df_repr(node.getSE());
		}
	}
	
	// create quadtree (capacity based only)
	public static void createQTree()
	{
		readSample();
		
		// create quad tree from sample dataset
        Node root = new Node(0.0, 0.0, 1.0, 1.0); // create root node
        root.setLow(0);
		root.setHigh(sampleArray.length - 1);
		
		// create quadtree
		makeQuadArray(sampleArray, root);
		
		// create quadtree in string form
		df_repr(root);
		
		System.out.printf("number of cells: %d\n", numCells);
		
		writeFiles(root);
	}
}
*/