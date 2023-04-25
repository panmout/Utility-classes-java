package gr.uth.ece.dsel.common_classes;

import java.util.ArrayList;
import java.util.PriorityQueue;

public final class FindNeighborsFunctions
{
    private PriorityQueue<IdDist> neighbors;

	// BF method
    public PriorityQueue<IdDist> getBfNeighbors(Point qpoint, ArrayList<Point> tpoints, int K)
	{
		this.neighbors = new PriorityQueue<>(K, new IdDistComparator("max")); // max heap of K neighbors;

	    for (Point tpoint : tpoints)
	    {	    	
	    	final double dist = UtilityFunctions.distance(qpoint, tpoint); // distance calculation
	    	final IdDist neighbor = new IdDist(tpoint.getId(), dist); // create neighbor
	    	
	    	// if PriorityQueue not full, add new tpoint (IdDist)
	    	if (this.neighbors.size() < K)
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
	} // end BF_Neighbors

	// PS method
	public PriorityQueue<IdDist> getPsNeighbors(Point qpoint, ArrayList<Point> tpoints, int K)
	{
		this.neighbors = new PriorityQueue<>(K, new IdDistComparator("max")); // max heap of K neighbors;
		
	 	if (!tpoints.isEmpty()) // if this cell has any tpoints
		{	
			final double x_left = tpoints.get(0).getX(); // leftmost tpoint's x
			final double x_right = tpoints.get(tpoints.size() - 1).getX(); // rightmost tpoint's x
			
			// check for neighbors to the left or to the right of the query point
			boolean check_right = false;
			boolean check_left = false;
			
			// tpoints list indexes
			int low = 0;
			int high = 0;
			
			if (qpoint.getX() < x_left) // qpoint is at left of all tpoints
				check_right = true;
			else if (x_right < qpoint.getX()) // qpoint is at right of all tpoints
			{
				check_left = true;
				high = tpoints.size() - 1;
			}
			else // qpoint is among tpoints
			{
				check_left = true;
				check_right = true;
				
				int tindex = UtilityFunctions.binarySearchTpoints(qpoint.getX(), tpoints); // get tpoints array index for qpoint interpolation
				low = tindex + 1;
				high = tindex;
			}
			
			boolean cont_search = true; // set flag to true
			
			if (check_right)
				while (low < tpoints.size() && cont_search) // scanning for neighbors to the right of tindex
					cont_search = psNeighbors(qpoint, tpoints, K, low++);
			
			cont_search = true; // reset flag to true
			
			if (check_left)
				while (high >= 0 && cont_search) // scanning for neighbors to the left of tindex
					cont_search = psNeighbors(qpoint, tpoints, K, high--);
		} // end if
	 	
	    return this.neighbors;
	 	// end PS_Neighbors
	}
	
	// calculates PS neighbors
	private boolean psNeighbors(Point qpoint, ArrayList<Point> tpoints, int K, int i)
	{
		final Point tpoint = tpoints.get(i); // get tpoint
		
		if (this.neighbors.size() < K) // if queue is not full, add new tpoints 
		{
			final double dist = UtilityFunctions.distance(qpoint, tpoint); // distance calculation
			final IdDist neighbor = new IdDist(tpoint.getId(), dist); // create neighbor
			
	    	if (!UtilityFunctions.isDuplicate(this.neighbors, neighbor))
				this.neighbors.offer(neighbor); // insert to queue
		}
		else  // if queue is full, run some checks and replace elements
		{
			final double dm = this.neighbors.peek().getDist(); // get (not remove) distance of neighbor with maximum distance
			
			if (UtilityFunctions.xDistance(qpoint, tpoint) > dm) // if tpoint's x distance is greater than neighbors max distance
				return false; // end for loop, no other points at this side will have smaller distance
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
		return true; // neighbors updated, proceed to next qpoint
	} // end psNeighbors
}
