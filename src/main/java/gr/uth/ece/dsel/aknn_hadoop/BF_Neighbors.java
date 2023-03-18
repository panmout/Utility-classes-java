package gr.uth.ece.dsel.aknn_hadoop;

import java.util.ArrayList;
import java.util.PriorityQueue;

import gr.uth.ece.dsel.UtilityFunctions;
import gr.uth.ece.dsel.common_classes.IdDist;
import gr.uth.ece.dsel.common_classes.IdDistComparator;
import gr.uth.ece.dsel.common_classes.Point;
import org.apache.hadoop.mapreduce.Reducer.Context;

public final class BF_Neighbors
{
	private final ArrayList<Point> tpoints;
	private final int k;
	private final PriorityQueue<IdDist> neighbors;

	public BF_Neighbors(ArrayList<Point> tp, int K, Context con)
	{
		this.tpoints = new ArrayList<>(tp);
		this.k = K;
		this.neighbors = new PriorityQueue<>(this.k, new IdDistComparator("max")); // max heap of K neighbors
	}
	
	public PriorityQueue<IdDist> getNeighbors(Point qpoint)
	{
		this.neighbors.clear();
		
	    for (Point tpoint : this.tpoints)
	    {	    	
	    	final double dist = UtilityFunctions.distance(qpoint, tpoint); // distance calculation
	    	final IdDist neighbor = new IdDist(tpoint.getId(), dist); // create neighbor
	    	
	    	// if PriorityQueue not full, add new tpoint (IdDist)
	    	if (this.neighbors.size() < this.k)
	    	{
		    	if (!UtilityFunctions.isDuplicate(this.neighbors, neighbor))
		    		this.neighbors.offer(neighbor); // insert to queue
	    	}
	    	else // if queue is full, run some checks and replace elements
	    	{
	    		final double dm = this.neighbors.peek().getDist(); // get (not remove) distance of neighbor with maximum distance
	    		
  				if (dist < dm) // compare distance
  				{  					
  					if (!UtilityFunctions.isDuplicate(this.neighbors, neighbor))
  					{
  						this.neighbors.poll(); // remove top element
  			    		this.neighbors.offer(neighbor); // insert to queue
  					}
  				} // end if
	    	} // end else
		} // end tpoints for
	    
	    return this.neighbors;
	} // end gdBfNeighbors
}