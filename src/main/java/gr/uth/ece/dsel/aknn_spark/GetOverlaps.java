package gr.uth.ece.dsel.aknn_spark;

import org.apache.spark.api.java.function.Function;
import scala.Tuple2;
import gr.uth.ece.dsel.common_classes.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;

public final class GetOverlaps implements Function<Tuple2<String, ArrayList<Tuple2<Point, PriorityQueue<IdDist>>>>, ArrayList<Tuple2<String, Point>>>
{
	private final HashMap<String, Integer> cell_tpoints; // hashmap of training points per cell list from Phase 1 <cell_id, number of training points>
	private final String partitioning; // gd or qt
	private int N; // N*N or N*N*N cells
	private final int K; // K neighbors
	private Node root; // create root node
	
	public GetOverlaps (HashMap<String, Integer> cell_tpoints, int k, String partitioning)
	{
		this.K = k;
		this.partitioning = partitioning;
		this.cell_tpoints = cell_tpoints;
	}
	
	public void setN (int n)
	{
		this.N = n;
	}
	
	public void setRoot (Node root)
	{
		this.root = root;
	}
	
	@Override
	public ArrayList<Tuple2<String, Point>> call (Tuple2<String, ArrayList<Tuple2<Point, PriorityQueue<IdDist>>>> input)
	{
		final String qcell = input._1; // query points' cell

		HashSet<String> overlaps = new HashSet<>(); // list of overlapped cells

		PriorityQueue<IdDist> neighbors ; // neighbors list

		final ArrayList<Tuple2<String, Point>> outValue = new ArrayList<>(); // output value: <overlapped cell, qpoint>

		final GetOverlapsFunctions ovl = new GetOverlapsFunctions (this.K, this.cell_tpoints);

		if (this.partitioning.equals("gd"))
			ovl.setN(this.N);
		else if (this.partitioning.equals("qt"))
			ovl.setRoot(this.root);
		
		for (Tuple2<Point, PriorityQueue<IdDist>> tuple: input._2) // read next <query point, neighbors list> tuple
		{
			Point qpoint = tuple._1; // get qpoint

			neighbors = tuple._2; // get its neighbor list
			
			// call appropriate overlaps discovery function according to partitioning
			if (this.partitioning.equals("gd"))
				overlaps = ovl.getOverlapsGD(qcell, qpoint, neighbors); // find overlaps
			else if (this.partitioning.equals("qt"))
				overlaps = ovl.getOverlapsQT(qcell, qpoint, neighbors); // find overlaps
			
			// if overlaps contains only qcell, qpoint gets 'true' status and is not going to output
			boolean listComplete = overlaps.size() == 1 && overlaps.contains(qcell);

			// fill output list with tuples <cell, qpoint> only for "false" qpoints
			if (!listComplete)
				for (String cell: overlaps)
					outValue.add(new Tuple2<>(cell, qpoint));
		}
		
		return outValue;
	}
}