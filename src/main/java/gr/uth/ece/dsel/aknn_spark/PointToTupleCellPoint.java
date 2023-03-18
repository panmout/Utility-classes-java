package gr.uth.ece.dsel.aknn_spark;

import org.apache.spark.api.java.function.PairFunction;
import gr.uth.ece.dsel.common_classes.*;
import gr.uth.ece.dsel.UtilityFunctions;
import scala.Tuple2;

public final class PointToTupleCellPoint implements PairFunction<Point, String, Point>
{
	private int N = 0;
	private Node node = null;
	
	public PointToTupleCellPoint (int n)
	{
		this.N = n;
	}
	
	public PointToTupleCellPoint (Node node)
	{
		this.node = node;
	}
	
	@Override
	public Tuple2<String, Point> call(Point p)
	{
		String cell = null;
		if (this.N != 0)
			cell = UtilityFunctions.pointToCellGD(p, this.N);
		else if (this.node != null)
		{
			if (this.node.getCNE() == null) // 2d
				cell = UtilityFunctions.pointToCellQT(p.getX(), p.getY(), this.node);
			else // 3d
				cell = UtilityFunctions.pointToCellQT(p.getX(), p.getY(), p.getZ(), this.node);
		}
		return new Tuple2<>(cell, p);
	}
}
