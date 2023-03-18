package gr.uth.ece.dsel.common_classes;

import java.io.IOException;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.ObjectOutputStream;
import java.io.FileOutputStream;
import java.util.Formatter;
import java.util.FormatterClosedException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;

import gr.uth.ece.dsel.common_classes.Node;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;

public final class CreateQTree
{
	private Formatter outputTextFile; // local output text file
	private FileOutputStream fileout;
	private ObjectOutputStream outputObjectFile; // local output object file
	private final int capacity;
	private HashMap<Integer, Double[]> sample_dataset;
	private int numCells = 0;
	private String x = "";
	private final String treeFilePath;
	private final String treeFileName;
	private final String trainingDatasetPath;
	private final int samplerate;
	private final HashSet<Double> widths = new HashSet<>(); // store cell widths
	private boolean is3d = false; // 2d or 3d quad tree (set in readSample() method)
	
	public CreateQTree(int newCapacity, String newTreeFilePath, String newTreeFileName, String newTrainingDatasetPath, int newSamplerate)
	{
		this.capacity = newCapacity;
		this.treeFilePath = newTreeFilePath;
		this.treeFileName = newTreeFileName;
		this.trainingDatasetPath = newTrainingDatasetPath;
		this.samplerate = newSamplerate;
	}
	
	// create root node (maximum capacity method)
	private Node createQT(Node node)
	{
		if (node.getContPoints().size() > this.capacity)
		{
			node = createChildren(node);
			
			if (!this.is3d) // 2d case
			{
				node.setNW(createQT(node.getNW()));
				node.setNE(createQT(node.getNE()));
				node.setSW(createQT(node.getSW()));
				node.setSE(createQT(node.getSE()));
			}
			else if (this.is3d) // 3d case
			{
				node.setCNW(createQT(node.getCNW()));
				node.setCNE(createQT(node.getCNE()));
				node.setCSW(createQT(node.getCSW()));
				node.setCSE(createQT(node.getCSE()));
				
				node.setFNW(createQT(node.getFNW()));
				node.setFNE(createQT(node.getFNE()));
				node.setFSW(createQT(node.getFSW()));
				node.setFSE(createQT(node.getFSE()));
			}
		}
		return node;
	}
	
	// create root node (all children split method)
	private Node createQT(Node node, boolean force)
	{
		if (node.getContPoints().size() > this.capacity || (force))
		{
			int nodeNumPoints = node.getContPoints().size();
			
			node = createChildren(node);
			
			if (nodeNumPoints <= this.capacity)
			{
				if (!this.is3d) // 2d case
				{
					node.setNW(createQT(node.getNW(), false));
					node.setNE(createQT(node.getNE(), false));
					node.setSW(createQT(node.getSW(), false));
					node.setSE(createQT(node.getSE(), false));
				}
				else if (this.is3d) // 3d case
				{
					node.setCNW(createQT(node.getCNW(), false));
					node.setCNE(createQT(node.getCNE(), false));
					node.setCSW(createQT(node.getCSW(), false));
					node.setCSE(createQT(node.getCSE(), false));
					
					node.setFNW(createQT(node.getFNW(), false));
					node.setFNE(createQT(node.getFNE(), false));
					node.setFSW(createQT(node.getFSW(), false));
					node.setFSE(createQT(node.getFSE(), false));
				}
			}
			else
			{
				force = false;
				
				if (!this.is3d) // 2d case
				{
					if (node.getNW().getContPoints().size() > this.capacity)
						force = true;
					if (node.getNE().getContPoints().size() > this.capacity)
						force = true;
					if (node.getSW().getContPoints().size() > this.capacity)
						force = true;
					if (node.getSE().getContPoints().size() > this.capacity)
						force = true;
					
					node.setNW(createQT(node.getNW(), force));
					node.setNE(createQT(node.getNE(), force));
					node.setSW(createQT(node.getSW(), force));
					node.setSE(createQT(node.getSE(), force));
				}
				else if (this.is3d) // 3d case
				{
					if (node.getCNW().getContPoints().size() > this.capacity)
						force = true;
					if (node.getCNE().getContPoints().size() > this.capacity)
						force = true;
					if (node.getCSW().getContPoints().size() > this.capacity)
						force = true;
					if (node.getCSE().getContPoints().size() > this.capacity)
						force = true;
					
					if (node.getFNW().getContPoints().size() > this.capacity)
						force = true;
					if (node.getFNE().getContPoints().size() > this.capacity)
						force = true;
					if (node.getFSW().getContPoints().size() > this.capacity)
						force = true;
					if (node.getFSE().getContPoints().size() > this.capacity)
						force = true;
					
					node.setCNW(createQT(node.getCNW(), force));
					node.setCNE(createQT(node.getCNE(), force));
					node.setCSW(createQT(node.getCSW(), force));
					node.setCSE(createQT(node.getCSE(), force));
					
					node.setFNW(createQT(node.getFNW(), force));
					node.setFNE(createQT(node.getFNE(), force));
					node.setFSW(createQT(node.getFSW(), force));
					node.setFSE(createQT(node.getFSE(), force));
				}
			}
		}
		return node;
	}
	
	// create root node (average width split method)
	private Node createQT(Node node, double avgWidth)
	{
		if ((node.getContPoints().size() > this.capacity) || (node.getXmax() - node.getXmin() > avgWidth)) // divide node only if it has many points or is bigger than average size
		{
			node = createChildren(node);
			
			if (!this.is3d) // 2d case
			{
				node.setNW(createQT(node.getNW(), avgWidth));
				node.setNE(createQT(node.getNE(), avgWidth));
				node.setSW(createQT(node.getSW(), avgWidth));
				node.setSE(createQT(node.getSE(), avgWidth));
			}
			else if (this.is3d) // 3d case
			{
				node.setCNW(createQT(node.getCNW(), avgWidth));
				node.setCNE(createQT(node.getCNE(), avgWidth));
				node.setCSW(createQT(node.getCSW(), avgWidth));
				node.setCSE(createQT(node.getCSE(), avgWidth));
				
				node.setFNW(createQT(node.getFNW(), avgWidth));
				node.setFNE(createQT(node.getFNE(), avgWidth));
				node.setFSW(createQT(node.getFSW(), avgWidth));
				node.setFSE(createQT(node.getFSE(), avgWidth));
			}
		}
		return node;
	}
	
	// create children and split sample training points
	private Node createChildren(Node node)
	{
		// define x, y
		double xmin = node.getXmin();
		double xmax = node.getXmax();
		double xmid = (xmin + xmax) / 2;
		
		double ymin = node.getYmin();
		double ymax = node.getYmax();
		double ymid = (ymin + ymax) / 2;
		
		if (!this.is3d) // 2d case
		{
			// create child nodes
			node.setNW(new Node(xmin, ymid, xmid, ymax));
			node.setNE(new Node(xmid, ymid, xmax, ymax));
			node.setSW(new Node(xmin, ymin, xmid, ymid));
			node.setSE(new Node(xmid, ymin, xmax, ymid));
			
			// partition dataset to child nodes
			Iterator<Integer> iterator = node.getContPoints().iterator(); // create iterator
			while (iterator.hasNext()) // while set has elements
			{
				int pid = iterator.next();
				Double[] coords = this.sample_dataset.get(pid);
				double x = coords[0];
				double y = coords[1];
				if (x >= xmin && x < xmid) // point inside SW or NW
				{
					if (y >= ymin && y < ymid) // point inside SW
						node.getSW().addPoints(pid);
					else if (y >= ymid && y < ymax) // point inside NW
						node.getNW().addPoints(pid);
				}
				else if (x >= xmid && x < xmax) // point inside SE or NE
				{
					if (y >= ymin && y < ymid) // point inside SE
						node.getSE().addPoints(pid);
					else if (y >= ymid && y < ymax) // point inside NE
						node.getNE().addPoints(pid);
				}
				iterator.remove();
				node.removePoint(pid); // remove point from parent node
			}
		}
		else if (this.is3d) // 3d case
		{
			// define z
			double zmin = node.getZmin();
			double zmax = node.getZmax();
			double zmid = (zmin + zmax) / 2;
			
			// create child nodes
			// CNW (xmin, ymid, zmid, xmid, ymax, zmax)
			node.setCNW(new Node(xmin, ymid, zmid, xmid, ymax, zmax));
			// CNE (xmid, ymid, zmid, xmax, ymax, zmax)
			node.setCNE(new Node(xmid, ymid, zmid, xmax, ymax, zmax));
			// CSW (xmin, ymin, zmid, xmid, ymid, zmax)
			node.setCSW(new Node(xmin, ymin, zmid, xmid, ymid, zmax));
			// CSE (xmid, ymin, zmid, xmax, ymid, zmax)
			node.setCSE(new Node(xmid, ymin, zmid, xmax, ymid, zmax));
			// FNW (xmin, ymid, zmin, xmid, ymax, zmid)
			node.setFNW(new Node(xmin, ymid, zmin, xmid, ymax, zmid));
			// FNE (xmid, ymid, zmin, xmax, ymax, zmid)
			node.setFNE(new Node(xmid, ymid, zmin, xmax, ymax, zmid));
			// FSW (xmin, ymin, zmin, xmid, ymid, zmid)
			node.setFSW(new Node(xmin, ymin, zmin, xmid, ymid, zmid));
			// FSE (xmid, ymin, zmin, xmax, ymid, zmid)
			node.setFSE(new Node(xmid, ymin, zmin, xmax, ymid, zmid));
			
			// partition dataset to child nodes
			Iterator<Integer> iterator = node.getContPoints().iterator(); // create iterator
			while (iterator.hasNext()) // while set has elements
			{
				int pid = iterator.next();
				Double[] coords = this.sample_dataset.get(pid);
				double x = coords[0];
				double y = coords[1];
				double z = coords[2];
				
				if (x >= xmin && x < xmid) // point inside SW or NW (Floor or Ceiling)
				{
					if (y >= ymin && y < ymid) // point inside SW (Floor or Ceiling)
					{
						if (z >= zmin && z < zmid) // point inside FSW
							node.getFSW().addPoints(pid);
						else if (z >= zmid && z < zmax) // point inside CSW
							node.getCSW().addPoints(pid);
					}
					else if (y >= ymid && y < ymax) // point inside NW (Floor or Ceiling)
					{
						if (z >= zmin && z < zmid) // point inside FNW
							node.getFNW().addPoints(pid);
						else if (z >= zmid && z < zmax) // point inside CNW
							node.getCNW().addPoints(pid);
					}
				}
				else if (x >= xmid && x < xmax) // point inside SE or NE (Floor or Ceiling)
				{
					if (y >= ymin && y < ymid) // point inside SE (Floor or Ceiling)
					{
						if (z >= zmin && z < zmid) // point inside FSE
							node.getFSE().addPoints(pid);
						else if (z >= zmid && z < zmax) // point inside CSE
							node.getCSE().addPoints(pid);
					}
					else if (y >= ymid && y < ymax) // point inside NE (Floor or Ceiling)
					{
						if (z >= zmin && z < zmid) // point inside FNE
							node.getFNE().addPoints(pid);
						else if (z >= zmid && z < zmax) // point inside CNE
							node.getCNE().addPoints(pid);
					}
				}
				iterator.remove();
				node.removePoint(pid); // remove point from parent node
			}
		}
		return node;
	}
	
	private void df_repr(Node node) // create qtree in string form
	{
		if (!this.is3d) // 2d case
		{
			if (node.getNW() == null)
			{
				this.x = this.x.concat("0");
				this.numCells++;
			}
			else
			{
				this.x = this.x.concat("1");
				
				df_repr(node.getNW());
				df_repr(node.getNE());
				df_repr(node.getSW());
				df_repr(node.getSE());
			}
		}
		else if (this.is3d) // 3d case
		{
			if (node.getCNW() == null)
			{
				this.x = this.x.concat("0");
				this.numCells++;
			}
			else
			{
				this.x = this.x.concat("1");
				
				df_repr(node.getFNW());
				df_repr(node.getFNE());
				df_repr(node.getFSW());
				df_repr(node.getFSE());
				
				df_repr(node.getCNW());
				df_repr(node.getCNE());
				df_repr(node.getCSW());
				df_repr(node.getCSE());
			}
		}
	}
	
	// get leaves widths
	private void getWidths(Node node)
	{		
		if (node.getNW() == null || node.getCNW() == null)
			this.widths.add(node.getXmax() - node.getXmin());
		else
		{
			if (!this.is3d) // 2d case
			{
				getWidths(node.getNW());
				getWidths(node.getNE());
				getWidths(node.getSW());
				getWidths(node.getSE());
			}
			else if (this.is3d) // 3d case
			{
				getWidths(node.getFNW());
				getWidths(node.getFNE());
				getWidths(node.getFSW());
				getWidths(node.getFSE());
				
				getWidths(node.getCNW());
				getWidths(node.getCNE());
				getWidths(node.getCSW());
				getWidths(node.getCSE());
			}
		}
	}
	
	private void readSample()
	{
		try // open files
		{
			this.fileout = new FileOutputStream(this.treeFileName);
			this.outputObjectFile = new ObjectOutputStream(this.fileout); // open local output object file
			this.outputTextFile = new Formatter("qtree.txt"); // open local output text file
			
			FileSystem fs = FileSystem.get(new Configuration());
			Path trainingPath = new Path(this.trainingDatasetPath);
			BufferedReader trainingBr = new BufferedReader(new InputStreamReader(fs.open(trainingPath))); // open HDFS training dataset file
			
			this.sample_dataset = new HashMap<>();
			
			HashSet<Integer> randomNumbers = new HashSet<>(this.samplerate); // [percentSample] size set for random integers
			
			Random random = new Random();
			
			while (randomNumbers.size() < this.samplerate) // fill list
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
					
					Double[] tpoint = null;
					
					if (data.length == 3) // 2d dataset
						tpoint = new Double[]{x, y};
					else if (data.length == 4) // 3d dataset
					{
						double z = Double.parseDouble(data[3]); // get z
						tpoint = new Double[]{x, y, z};
						this.is3d = true;
					}
					this.sample_dataset.put(pid, tpoint); // add {pid, x, y, z} to hashmap
				}
			}
		}
		catch (IOException ioException)
		{
			System.err.println("Could not open file, exiting");
			System.exit(1);
		}
	}
	
	private void writeFiles(Node node)
	{		
		// write to files
		try
		{
			// local
			this.outputTextFile = new Formatter("qtree.txt");
			this.outputTextFile.format("%s", this.x);
			this.outputObjectFile.writeObject(node);
			
			this.outputObjectFile.close();
			this.outputTextFile.close();
			this.fileout.close();
			
			// write to hdfs
			FileSystem fs = FileSystem.get(new Configuration());
			Path path = new Path(this.treeFilePath);
			ObjectOutputStream outputStream = new ObjectOutputStream(fs.create(path));
			outputStream.writeObject(node);
			outputStream.close();
		}
		catch (FormatterClosedException formatterException)
		{
			System.err.println("Error writing to file, exiting");
			System.exit(2);
		}
		catch (IOException ioException)
		{
			System.err.println("Error writing to file, exiting");
			System.exit(2);
		}
	}
	
	// create quadtree (capacity based only)
	public void createQTree()
	{
		readSample();
		
		Node root = null;
		
		// create 2d or 3d quad tree from sample dataset
		if (!this.is3d) // 2d
			root = new Node(0.0, 0.0, 1.0, 1.0); // create root node
		else if (this.is3d) // 3d
			root = new Node(0.0, 0.0, 0.0, 1.0, 1.0, 1.0); // create root node
		
		for (int i : this.sample_dataset.keySet()) // add all sample points to root node
			root.addPoints(i);
		
		root = createQT(root); // create tree from root
		
		root.removePoints(); // remove all tpoints from tree
		
		df_repr(root); // create tree string
		
		System.out.printf("number of cells: %d\n", this.numCells);
		
		writeFiles(root);
	}
	
	// create quadtree (all children split method = if one child splits, all will)
	public void createAllChldSplitQTree()
	{
		readSample();
		
		Node root = null;
		
		// create 2d or 3d quad tree from sample dataset
		if (!this.is3d) // 2d
			root = new Node(0.0, 0.0, 1.0, 1.0); // create root node
		else if (this.is3d) // 3d
			root = new Node(0.0, 0.0, 0.0, 1.0, 1.0, 1.0); // create root node
		
		for (int i : this.sample_dataset.keySet()) // add all sample points to root node
			root.addPoints(i);
		
		root = createQT(root, false); // create tree from root
		
		root.removePoints(); // remove all tpoints from tree
		
		df_repr(root); // create tree string
		
		System.out.printf("number of cells: %d\n", this.numCells);
		
		writeFiles(root);
	}
	
	// create quadtree (capacity based only)
	public void createAvgWidthQTree()
	{
		readSample();
		
		Node root1 = null;
		
		// create 2d or 3d quad tree from sample dataset
		if (!this.is3d) // 2d
			root1 = new Node(0.0, 0.0, 1.0, 1.0); // create root node
		else if (this.is3d) // 3d
			root1 = new Node(0.0, 0.0, 0.0, 1.0, 1.0, 1.0); // create root node
		
		for (int i : this.sample_dataset.keySet()) // add all sample points to root node
			root1.addPoints(i);
		
		// create initial tree from root, 'POSITIVE_INFINITY' is used to make the right OR statement 'false'
		root1 = createQT(root1, Double.POSITIVE_INFINITY);
		
		getWidths(root1); // get all widths of leaves
		
		double averageWidth = 0;
		
		for (double i : this.widths)
			averageWidth += i;
		
		averageWidth = averageWidth / this.widths.size();
		
		System.out.printf("average width: %9.8f\n", averageWidth);
		
		Node root2 = null;
		
		// create 2d or 3d quad tree from sample dataset
		if (!this.is3d) // 2d
			root2 = new Node(0.0, 0.0, 1.0, 1.0); // create root node
		else if (this.is3d) // 3d
			root2 = new Node(0.0, 0.0, 0.0, 1.0, 1.0, 1.0); // create root node
		
		for (int i : this.sample_dataset.keySet()) // add all sample points to root node
			root2.addPoints(i);
		
		root2 = createQT(root2, averageWidth); // create final tree from root and average width
		
		root2.removePoints(); // remove all tpoints from tree
		
		df_repr(root2); // create tree string
		
		System.out.printf("number of cells: %d\n", this.numCells);
		
		writeFiles(root2);
	}
}