/*
package gr.uth.ece.dsel.common_classes;

public final class QuadtreeArray
{
	private static String trainingDataset; // training dataset name in HDFS
	private static String nameNode; // hostname
    private static String trainingDir; // HDFS dir containing training dataset
    private static int samplerate; // percent sample of dataset
	private static int capacity; // quad tree node maximum capacity
    private static String treeDir; // HDFS dir containing sample trees

    public static void main(String[] args)
	{
		long t0 = System.currentTimeMillis();
		
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
        // full hdfs path name for qtree object file
        String treeFilePath = String.format("hdfs://%s:9000/user/%s/%s/%s", nameNode, username, treeDir, treeFileName);

        // array object file name
        String arrayFileName = "qtreeArray.ser";
        // full hdfs path name for qtree array file
        String arrayFilePath = String.format("hdfs://%s:9000/user/%s/%s/%s", nameNode, username, treeDir, arrayFileName);
		
		
		new CreateQTreeArray(capacity, treeFilePath, treeFileName, arrayFilePath, arrayFileName, trainingDatasetPath, samplerate);
		
		CreateQTreeArray.createQTree();
		
		long treetime = System.currentTimeMillis() - t0;
		
		System.out.printf("Quadtree {capacity: %d, samplerate: %d} creation time: %d millis\n", capacity, samplerate, treetime);
	}
}
*/