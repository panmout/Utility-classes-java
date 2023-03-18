package gr.uth.ece.dsel.aknn_spark;

import org.apache.spark.api.java.function.Function;
import scala.Tuple2;
import gr.uth.ece.dsel.common_classes.*;
import gr.uth.ece.dsel.UtilityFunctions;
import java.util.ArrayList;
import java.util.PriorityQueue;

public final class PS_Neighbors implements Function<Tuple2<Iterable<Point>, Iterable<Point>>, ArrayList<Tuple2<Point, PriorityQueue<IdDist>>>>
{
	private final int k;
	private final ArrayList<Point> qpoints;
	private final ArrayList<Point> tpoints;
	private final PriorityQueue<IdDist> neighbors;
	private final ArrayList<Tuple2<Point, PriorityQueue<IdDist>>> qpoint_neighbors;
	
	public PS_Neighbors(int k)
	{
		this.k = k;
		this.neighbors = new PriorityQueue<>(this.k, new IdDistComparator("max"));
		this.qpoint_neighbors = new ArrayList<>();
		this.qpoints = new ArrayList<>();
		this.tpoints = new ArrayList<>();
	}
	
	@Override
	public ArrayList<Tuple2<Point, PriorityQueue<IdDist>>> call (Tuple2<Iterable<Point>, Iterable<Point>> tuple)
	{
		this.qpoint_neighbors.clear();
		
		this.qpoints.clear();
		this.tpoints.clear();
	    
		// fill points lists
	    for (Point qpoint: tuple._1)
	    	this.qpoints.add(qpoint);
	    
	    for (Point tpoint: tuple._2)
	    	this.tpoints.add(tpoint);
	    
	    // sort lists by x ascending
	    this.qpoints.sort(new PointXYComparator("min", 'x'));
	 	this.tpoints.sort(new PointXYComparator("min", 'x'));
	 	
		 // traverse query points
		for (Point qpoint: this.qpoints)
	 	{
  			if (!this.tpoints.isEmpty()) // if this cell has any tpoints
	 		{
	 			// get pqpoint's x
	 			final double xq = qpoint.getX();
	 			
	 			final double x_left = this.tpoints.get(0).getX(); // leftmost tpoint's x
				final double x_right = this.tpoints.get(this.tpoints.size() - 1).getX(); // rightmost tpoint's x
				
				boolean check_right = false;
				boolean check_left = false;
				
				int low = 0;
				int high = 0;

				// end if
				if (xq < x_left) // pqpoint is at left of all tpoints
					check_right = true;
				else if (x_right < xq) // pqpoint is at right of all tpoints
				{
					check_left = true;
					high = this.tpoints.size() - 1;
				} // end else if
				else // pqpoint is among tpoints
				{
					check_left = true;
					check_right = true;
					
					int tindex = UtilityFunctions.binarySearchTpoints(xq, this.tpoints); // get tpoints array index for qpoint interpolation
					low = tindex + 1;
					high = tindex;
				} // end else
				
				boolean cont_search = true; // set flag to true
				
				if (check_right)
					while (low < this.tpoints.size() && cont_search) // scanning for neighbors to the right of tindex
						cont_search = getPsNeighbors(qpoint, low++);
				
				cont_search = true; // reset flag to true
				
				if (check_left)
					while (high >= 0 && cont_search) // scanning for neighbors to the left of tindex
						cont_search = getPsNeighbors(qpoint, high--);
	 		} // end empty tpoints check
  			
 			this.qpoint_neighbors.add(new Tuple2<>(qpoint, new PriorityQueue<>(this.neighbors)));
			this.neighbors.clear();
	 	} // end query points traverse
 		
 		// return <list of qpoints with their neighbor list>
	    return new ArrayList<>(this.qpoint_neighbors);
	} // end call
	
	// calculate neighbors
	private boolean getPsNeighbors(Point qpoint, int index)
	{
		final Point tpoint = this.tpoints.get(index); // get tpoint
		
		if (this.neighbors.size() < this.k) // if queue is not full, add new tpoints 
		{
			final double dist = UtilityFunctions.distance(qpoint, tpoint); // distance calculation
	    	final IdDist neighbor = new IdDist(tpoint.getId(), dist); // create neighbor
	    	
	    	if (!UtilityFunctions.isDuplicate(this.neighbors, neighbor))
	    		this.neighbors.offer(neighbor); // insert to queue
		} // end if
		else // if queue is full, run some checks and replace elements
		{
			final double dm = this.neighbors.peek().getDist(); // get (not remove) neighbor with maximum distance
			
			if (UtilityFunctions.xDistance(qpoint, tpoint) > dm)  // if tpoint's x distance is greater than neighbors max distance
				return false; // end for loop, no other points to the right will have smaller distance
			else
			{
				final double dist = UtilityFunctions.distance(qpoint, tpoint); // distance calculation
				
				if (dist < dm) // compare distance
				{
					final IdDist neighbor = new IdDist(tpoint.getId(), dist); // create neighbor
					
					if (!UtilityFunctions.isDuplicate(this.neighbors, neighbor))
					{
						this.neighbors.poll(); // remove top element
						this.neighbors.offer(neighbor); // insert to queue
					}
				} // end if
			} // end else
		} // end else
		return true; // neighbors updated, proceed to next ipoint
	}
}
