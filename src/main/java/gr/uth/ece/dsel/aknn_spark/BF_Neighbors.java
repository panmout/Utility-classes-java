package gr.uth.ece.dsel.aknn_spark;

import java.util.ArrayList;
import java.util.PriorityQueue;
import org.apache.spark.api.java.function.Function;
import gr.uth.ece.dsel.common_classes.*;
import gr.uth.ece.dsel.UtilityFunctions;
import scala.Tuple2;

public final class BF_Neighbors implements Function<Tuple2<Iterable<Point>, Iterable<Point>>, ArrayList<Tuple2<Point, PriorityQueue<IdDist>>>>
{
	private final int k;
	private final PriorityQueue<IdDist> neighbors;
	private final ArrayList<Tuple2<Point, PriorityQueue<IdDist>>> qpoint_neighbors;
	
	public BF_Neighbors(int k)
	{
		this.k = k;
		this.neighbors = new PriorityQueue<>(this.k, new IdDistComparator("max"));
		this.qpoint_neighbors = new ArrayList<>();
	}
	
	@Override
	public ArrayList<Tuple2<Point, PriorityQueue<IdDist>>> call (Tuple2<Iterable<Point>, Iterable<Point>> tuple)
	{
		this.qpoint_neighbors.clear();
		
		// traverse <query points, neighbors> tuples
	    for (Point qpoint: tuple._1)
	    {
    		// traverse training points
	    	for (Point tpoint: tuple._2)
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
			} // end training points traverse
	    	//System.out.printf("qpoint: %d, neighbors: %s\n", qpoint.getId(), AknnFunctions.pqToString(this.neighbors, this.k, "min"));
	    	this.qpoint_neighbors.add(new Tuple2<>(qpoint, new PriorityQueue<>(this.neighbors)));
	    	this.neighbors.clear();
	    } // end query points traverse
	    // return <list of qpoints with their neighbor list>
	    return new ArrayList<>(this.qpoint_neighbors);
	} // end gdBfNeighbors
}
