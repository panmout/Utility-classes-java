package gr.uth.ece.dsel.aknn_hadoop;

import java.util.ArrayList;
import java.util.PriorityQueue;
import gr.uth.ece.dsel.common_classes.*;
import org.apache.hadoop.mapreduce.Reducer.Context;

public final class FindNeighbors
{
	private final ArrayList<Point> tpoints;
	private final int k;

	public FindNeighbors(ArrayList<Point> tp, int K, Context con)
	{
		this.tpoints = new ArrayList<>(tp);
		this.k = K;
	}
	
	public PriorityQueue<IdDist> getBfNeighbors(Point qpoint)
	{
	    return FindNeighborsFunctions.getBfNeighbors(qpoint, this.tpoints, this.k);
	}

	public PriorityQueue<IdDist> getPsNeighbors(Point qpoint)
	{
	    return FindNeighborsFunctions.getPsNeighbors(qpoint, this.tpoints, this.k);
	}
}