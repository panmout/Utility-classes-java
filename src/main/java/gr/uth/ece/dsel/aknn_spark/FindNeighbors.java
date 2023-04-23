package gr.uth.ece.dsel.aknn_spark;

import org.apache.spark.api.java.function.Function;
import scala.Tuple2;
import gr.uth.ece.dsel.common_classes.*;
import java.util.ArrayList;
import java.util.PriorityQueue;

public final class FindNeighbors implements Function<Tuple2<Iterable<Point>, Iterable<Point>>, ArrayList<Tuple2<Point, PriorityQueue<IdDist>>>>
{
    private final int k;
    private final String method;
	
	public FindNeighbors(int k, String method)
	{
		this.k = k;
        this.method = method;
	}

    @Override
	public ArrayList<Tuple2<Point, PriorityQueue<IdDist>>> call (Tuple2<Iterable<Point>, Iterable<Point>> tuple)
	{
        final ArrayList<Point> qpoints = new ArrayList<>();
        final ArrayList<Point> tpoints = new ArrayList<>();
        final ArrayList<Tuple2<Point, PriorityQueue<IdDist>>> qpoint_neighbors = new ArrayList<>();

        // convert iterable to ArrayList
        tuple._1.forEach(qpoints::add);
		tuple._2.forEach(tpoints::add);

        if (this.method.equals("ps"))
        {
            // sort lists by x ascending
            qpoints.sort(new PointXYComparator("min", 'x'));
            tpoints.sort(new PointXYComparator("min", 'x'));
        }

        final FindNeighborsFunctions fnf = new FindNeighborsFunctions();

        // traverse query points
		for (Point qpoint: qpoints)
        {
            // get neighbors
            if (this.method.equals("bf"))
                qpoint_neighbors.add(new Tuple2<>(qpoint, fnf.getBfNeighbors(qpoint, tpoints, this.k)));
            else if (this.method.equals("ps"))
                qpoint_neighbors.add(new Tuple2<>(qpoint, fnf.getPsNeighbors(qpoint, tpoints, this.k)));
        }
        
        // return <list of qpoints with their neighbor list>
       return qpoint_neighbors;
    }
}
