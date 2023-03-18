package gr.uth.ece.dsel.common_classes;

import java.util.Comparator;

public final class PointXYComparator implements Comparator<Point>
{
	private final int type;
	private final char xy;
	
	public PointXYComparator(String type, char xy)
	{
		if (type.equals("min")) // x or y ascending comparator for Point objects
			this.type = -1;
		else if (type.equals("max")) // x or y descending comparator for Point objects
			this.type = 1;
		else throw new IllegalArgumentException("first argument must be 'min' or 'max'");
		
		if (xy == 'x')
			this.xy = 'x';
		else if (xy == 'y')
			this.xy = 'y';
		else throw new IllegalArgumentException("second argument must be 'x' or 'y'");
	}
	
	@Override
	public int compare(Point element1, Point element2)
	{
		switch (this.xy)
		{
			case 'x':
				if (element1.getX() < element2.getX())
					return this.type;
				else if (element1.getX() == element2.getX())
					return 0;
				else
					return -this.type;
			case 'y':
				if (element1.getY() < element2.getY())
					return this.type;
				else if (element1.getY() == element2.getY())
					return 0;
				else
					return -this.type;
			default:
				return 1000;
		}
	}
}
