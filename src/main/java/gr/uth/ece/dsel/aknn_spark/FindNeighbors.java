package gr.uth.ece.dsel.aknn_spark;

import org.apache.spark.api.java.function.Function;
import scala.Tuple2;
import gr.uth.ece.dsel.common_classes.*;
import java.util.ArrayList;
import java.util.PriorityQueue;

public final class FindNeighbors implements Function<Tuple2<Iterable<Point>, Iterable<Point>>, ArrayList<Tuple2<Point, PriorityQueue<IdDist>>>>
{
    private final int k;
	private final ArrayList<Point> qpoints;
	private final ArrayList<Point> tpoints;
	private final ArrayList<Tuple2<Point, PriorityQueue<IdDist>>> qpoint_neighbors;
    private final String method;
	
	public FindNeighbors(int k, String method)
	{
		this.k = k;
		this.qpoint_neighbors = new ArrayList<>();
		this.qpoints = new ArrayList<>();
		this.tpoints = new ArrayList<>();
        this.method = method;
	}

    @Override
	public ArrayList<Tuple2<Point, PriorityQueue<IdDist>>> call (Tuple2<Iterable<Point>, Iterable<Point>> tuple)
	{
        this.qpoint_neighbors.clear();
        this.qpoints.clear();
		this.tpoints.clear();

        // convert iterable to ArrayList
        tuple._1.forEach(this.qpoints::add);
		tuple._2.forEach(this.tpoints::add);

		PriorityQueue<IdDist> neighbors;

        if (method.equals("ps"))
        {	    
            // sort lists by x ascending
            this.qpoints.sort(new PointXYComparator("min", 'x'));
            this.tpoints.sort(new PointXYComparator("min", 'x'));
        }

        // traverse query points
		for (Point qpoint: this.qpoints)
        {
            // get neighbors
            neighbors = FindNeighborsFunctions.getPsNeighbors(qpoint, tpoints, this.k);
             
            this.qpoint_neighbors.add(new Tuple2<>(qpoint, new PriorityQueue<>(neighbors)));

            neighbors.clear();
        } // end query points traverse
        
        // return <list of qpoints with their neighbor list>
       return new ArrayList<>(this.qpoint_neighbors);
    }
}
