package gr.uth.ece.dsel.common_classes;

public final class Qtree
{
	private static String trainingDataset; // training dataset name in HDFS
	private static String nameNode; // hostname
	private static String trainingDir; // HDFS dir containing training dataset
	private static int samplerate; // percent sample of dataset
	private static int capacity; // quad tree node maximum capacity
	private static String treeDir; // HDFS dir containing sample trees
	private static int type; // 1 for simple capacity based quadtree, 2 for all children split method, 3 for average width method
	
	public static void main(String[] args)
	{
		final long t0 = System.currentTimeMillis();
		
		for (String arg: args)
		{
			String[] newarg;
			if (arg.contains("="))
			{
				newarg = arg.split("=");
				
				if (newarg[0].equals("nameNode"))
					nameNode = newarg[1];
				if (newarg[0].equals("trainingDir"))
					trainingDir = newarg[1];
				if (newarg[0].equals("treeDir"))
					treeDir = newarg[1];
				if (newarg[0].equals("trainingDataset"))
					trainingDataset = newarg[1];
				if (newarg[0].equals("samplerate"))
					samplerate = Integer.parseInt(newarg[1]);
				if (newarg[0].equals("capacity"))
					capacity = Integer.parseInt(newarg[1]);
				if (newarg[0].equals("type"))
					type = Integer.parseInt(newarg[1]);
			}
			else
				throw new IllegalArgumentException("not a valid argument, must be \"name=arg\", : " + arg);
		}

		// username
		String username = System.getProperty("user.name");
		// full HDFS path+name of training dataset
		String trainingDatasetPath = String.format("hdfs://%s:9000/user/%s/%s/%s", nameNode, username, trainingDir, trainingDataset);
		// tree file name
		String treeFileName = "qtree.ser";
		// full hdfs path name
		String treeFilePath = String.format("hdfs://%s:9000/user/%s/%s/%s", nameNode, username, treeDir, treeFileName);
		
		CreateQTree qtree = new CreateQTree(capacity, treeFilePath, treeFileName, trainingDatasetPath, samplerate);
		
		String qtreeType = "";
		
		switch(type)
		{
			case 1:
				qtree.createQTree();
				qtreeType = qtreeType.concat("maximum capacity method");
				break;
			case 2:
				qtree.createAllChldSplitQTree();
				qtreeType = qtreeType.concat("all children split method");
				break;
			case 3:
				qtree.createAvgWidthQTree();
				qtreeType = qtreeType.concat("average width method");
				break;
		}
		
		System.out.printf("Quadtree {%s, capacity: %d, samplerate: %d} creation time: %d millis\n", qtreeType, capacity, samplerate, System.currentTimeMillis() - t0);
	}
}