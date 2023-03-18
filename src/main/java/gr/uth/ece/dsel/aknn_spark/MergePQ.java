package gr.uth.ece.dsel.aknn_spark;

import java.util.PriorityQueue;
import gr.uth.ece.dsel.common_classes.*;
import gr.uth.ece.dsel.UtilityFunctions;
import org.apache.spark.api.java.function.Function;
import scala.Tuple2;

public final class MergePQ implements Function<Tuple2<Iterable<PriorityQueue<IdDist>>, Iterable<PriorityQueue<IdDist>>>, PriorityQueue<IdDist>>
{
	private final int k;
	private final PriorityQueue<IdDist> neighbors;
	
	public MergePQ (int k)
	{
		this.k = k;
		this.neighbors = new PriorityQueue<>(this.k, new IdDistComparator("max"));
	}
	
	@Override
	public PriorityQueue<IdDist> call (Tuple2<Iterable<PriorityQueue<IdDist>>, Iterable<PriorityQueue<IdDist>>> tuple)
	{
		this.neighbors.clear();
		
		// first iterable of priority queues are neighbor lists from phase 2
		for (PriorityQueue<IdDist> pq1 : tuple._1)
			joinPQ(pq1);
			
		// second iterable of priority queues are neighbor lists from phase 3
		for (PriorityQueue<IdDist> pq2 : tuple._2)
			joinPQ(pq2);
		
		while (this.neighbors.size() > this.k)
			this.neighbors.poll();
		
		return new PriorityQueue<>(this.neighbors);
	}
	
	private void joinPQ (PriorityQueue<IdDist> pq)
	{
		while (!pq.isEmpty())
		{
			IdDist neighbor = pq.poll();
			if (!UtilityFunctions.isDuplicate(this.neighbors, neighbor))
				this.neighbors.add(neighbor);
		}
	}
}
