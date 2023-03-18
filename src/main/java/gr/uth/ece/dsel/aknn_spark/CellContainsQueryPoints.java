package gr.uth.ece.dsel.aknn_spark;

import java.util.Iterator;

import gr.uth.ece.dsel.common_classes.Point;
import org.apache.spark.api.java.function.Function;
import scala.Tuple2;

public final class CellContainsQueryPoints implements Function<Tuple2<String, Tuple2<Iterable<Point>, Iterable<Point>>>, Boolean>
{
	@Override
	public Boolean call(Tuple2<String, Tuple2<Iterable<Point>, Iterable<Point>>> qtIterable)
	{
		Iterator<Point> qpoint = qtIterable._2._1.iterator(); // iterator on query points tuples

		return qpoint.hasNext();
	}
}
