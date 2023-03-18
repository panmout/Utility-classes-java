package gr.uth.ece.dsel.common_classes;

import java.io.Serializable;
import java.util.Comparator;

public final class IdDistComparator implements Comparator<IdDist>, Serializable
{
	private final int x;
	
	public IdDistComparator(String type)
	{
		if (type.equals("min")) // ascending comparator for IdDist objects (point_id, dist)
			this.x = -1;
		else if (type.equals("max")) // descending comparator for IdDist objects (point_id, dist)
			this.x = 1;
		else throw new IllegalArgumentException("argument must be 'min' or 'max'");
	}
	
	@Override
	public int compare(IdDist element1, IdDist element2)
	{
		if (element1.getDist() < element2.getDist())
			return this.x;
		else if (element1.getDist() == element2.getDist())
			return 0;
		else
			return -this.x;
	}
}

