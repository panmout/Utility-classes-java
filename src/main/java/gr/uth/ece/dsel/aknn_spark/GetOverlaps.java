package gr.uth.ece.dsel.aknn_spark;

import org.apache.spark.api.java.function.Function;
import scala.Tuple2;
import gr.uth.ece.dsel.common_classes.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;

public final class GetOverlaps implements Function<Tuple2<String, ArrayList<Tuple2<Point, PriorityQueue<IdDist>>>>, ArrayList<Tuple2<String, Tuple2<Point, Boolean>>>>
{
	private final HashMap<String, Integer> cell_tpoints; // hashmap of training points per cell list from Phase 1 <cell_id, number of training points>
	private final String partitioning; // gd or qt
	private int N; // N*N or N*N*N cells
	private final int K; // K neighbors
	private Node root; // create root node
	private final HashSet<String> overlaps; // list of overlapped cells
	private final PriorityQueue<IdDist> neighbors; // neighbors list
	private final ArrayList<Tuple2<String, Tuple2<Point, Boolean>>> outValue; // output value: <overlapped cell, qpoint, true/false>
	
	public GetOverlaps (HashMap<String, Integer> cell_tpoints, int k, String partitioning)
	{
		this.K = k;
		this.partitioning = partitioning;
		this.overlaps = new HashSet<>();
		this.neighbors = new PriorityQueue<>(this.K, new IdDistComparator("max"));
		this.outValue = new ArrayList<>();
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
	public ArrayList<Tuple2<String, Tuple2<Point, Boolean>>> call (Tuple2<String, ArrayList<Tuple2<Point, PriorityQueue<IdDist>>>> input)
	{
		final String qcell = input._1; // query points' cell
		
		// clear list
		this.outValue.clear();
		
		for (Tuple2<Point, PriorityQueue<IdDist>> tuple: input._2) // read next <query point, neighbors list> tuple
		{
			// clear lists
			this.overlaps.clear();
			this.neighbors.clear();
			
			Point qpoint = tuple._1; // get qpoint
			this.neighbors.addAll(tuple._2); // get its neighbor list
			
			// call appropriate overlaps discovery function according to partitioning
			if (this.partitioning.equals("gd"))
				this.overlaps.addAll(GetOverlapsFunctions.getOverlapsGD(qcell, qpoint, this.K, this.N, this.cell_tpoints, this.neighbors)); // find overlaps
			else if (this.partitioning.equals("qt"))
				this.overlaps.addAll(GetOverlapsFunctions.getOverlapsQT(qcell, qpoint, this.K, this.root, this.cell_tpoints, this.neighbors)); // find overlaps
			
			boolean listComplete = this.overlaps.size() == 1 && this.overlaps.contains(qcell); // if overlaps contains only qcell, qpoint gets 'true' status

			// fill output list with tuples <cell, <qpoint, true/false>>
			for (String cell: this.overlaps)
				this.outValue.add(new Tuple2<>(cell, new Tuple2<>(qpoint, listComplete)));
		}
		
		return new ArrayList<>(this.outValue);
	}
}