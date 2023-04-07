package gr.uth.ece.dsel;

import gr.uth.ece.dsel.common_classes.*;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.PriorityQueue;

public final class UtilityFunctions
{
	// return euclidean distance between two 2d points (x1, y1) and (x2, y2)
	public static double distance (double x1, double y1, double x2, double y2)
	{
		return Math.sqrt(square_distance (x1, y1, x2, y2));
	}// end euclidean distance
	
	// return euclidean distance between two 3d points (x1, y1, z1) and (x2, y2, z2)
	public static double distance (double x1, double y1, double z1, double x2, double y2, double z2)
	{
		return Math.sqrt(square_distance (x1, y1, z1, x2, y2, z2));
	}// end euclidean distance
	
	// return euclidean distance between two points
	public static double distance (Point ipoint, Point tpoint)
	{
		// 3d points
		if (ipoint.getZ() != Double.NEGATIVE_INFINITY && tpoint.getZ() != Double.NEGATIVE_INFINITY)
			return distance(ipoint.getX(), ipoint.getY(), ipoint.getZ(), tpoint.getX(), tpoint.getY(), tpoint.getZ());
		else // 2d points
			return distance(ipoint.getX(), ipoint.getY(), tpoint.getX(), tpoint.getY());
	}// end euclidean distance
	
	// return square of euclidean distance between two 2d points (x1, y1) and (x2, y2)
	public static double square_distance (double x1, double y1, double x2, double y2)
	{
		return Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2);
	}// end square of euclidean distance
	
	// return square of euclidean distance between two 3d points (x1, y1, z1) and (x2, y2, z2)
	public static double square_distance (double x1, double y1, double z1, double x2, double y2, double z2)
	{
		return Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2) + Math.pow(z1 - z2,  2);
	}// end square of euclidean distance
	
	// return x-distance between two points
	public static double xDistance (Point qpoint, Point tpoint)
	{
		return Math.abs(qpoint.getX() - tpoint.getX());
	}
		
	// String to point
	public static Point stringToPoint(String line, String sep)
	{
		final String[] data = line.trim().split(sep);
		final int id = Integer.parseInt(data[0]);
		final double x = Double.parseDouble(data[1]);
		final double y = Double.parseDouble(data[2]);
		
		if (data.length == 3) // 2d point
			return new Point(id, x, y);
		else if (data.length == 4) // 3d point
		{
			final double z = Double.parseDouble(data[3]);
			return new Point(id, x, y, z);
		}
		else
			throw new IllegalArgumentException();
	}
	
	/*
	Cell array (numbers inside cells are cell_id)

	    n*ds |---------|----------|----------|----------|----------|----------|--------------|
	         | (n-1)n  | (n-1)n+1 | (n-1)n+2 |          | (n-1)n+i |          | (n-1)n+(n-1) |
	(n-1)*ds |---------|----------|----------|----------|----------|----------|--------------|
	         |         |          |          |          |          |          |              |
	         |---------|----------|----------|----------|----------|----------|--------------|
	         |   j*n   |  j*n+1   |  j*n+2   |          |  j*n+i   |          |   j*n+(n-1)  |
	    j*ds |---------|----------|----------|----------|----------|----------|--------------|
	         |         |          |          |          |          |          |              |
	         |---------|----------|----------|----------|----------|----------|--------------|
	         |   2n    |   2n+1   |   2n+2   |          |   2n+i   |          |     3n-1     |
	    2*ds |---------|----------|----------|----------|----------|----------|--------------|
	         |    n    |    n+1   |    n+2   |          |    n+i   |          |     2n-1     |
	      ds |---------|----------|----------|----------|----------|----------|--------------|
	         |    0    |     1    |     2    |          |     i    |          |      n-1     |
	         |---------|----------|----------|----------|----------|----------|--------------|
	       0          ds         2*ds                  i*ds               (n-1)*ds          n*ds


	point to 2d GD cell
	-------------------
	cell_id(i,j) = j*n+i
	
	ds = 1.0/N;
	i = (int) (x/ds)
	j = (int) (y/ds)
	cell = j*N+i
	
	point to 3d GD cell
	-------------------
	z dimension extends vertically towards the user

	we continue the numbering from the top right cell by adding its next, which is (n-1)n+(n-1)+1 = n*n-n+n-1+1 = n*n,
	to the bottom left 0 and so on, every cell is just increased by n*n

	the first level cell above j*n+i will be i+j*n+n*n

	the second level cell above j*n+i will be i+j*n+2*n*n

	the k-th level cell above j*n+i will be i+j*n+k*n*n

	So, cell_id(i,j,k) = i+j*n+k*n*n

	How to get i,j,k from cell_id:

	cell = i+(j+k*n)*n --> i = mod(cell, n)

	cell-i = j*n+k*n*n --> (cell-i)/n = j+k*n --> j = mod{(cell-i)/n, n}

	k = {(cell-i)/n-j}/n
	*/
	
	public static String pointToCellGD(Point p, int n)
	{
		final double ds = 1.0 / n; // interval ds (cell width)
		final double x = p.getX();  // p.x
		final double y = p.getY();  // p.y
		final int i = (int) (x / ds); // i = (int) x/ds
		final int j = (int) (y / ds); // j = (int) y/ds
		int cellId;
		
		if (p.getZ() == Double.NEGATIVE_INFINITY) // 2d point
			cellId = j * n + i;
		else // 3d point
		{
			final double z = p.getZ(); // p.z
			final int k = (int) (z / ds); // k = (int) z/ds
			cellId = i + j * n + k * n * n;
		}
		return String.valueOf(cellId); // return cellId
	}
	
//	// node to cell
//	public static final String nodeToCell(Node node)
//	{
//		return pointToCellQT((node.getXmin() + node.getXmax()) / 2, (node.getYmin() + node.getYmax()) / 2, node);
//	}
	
	// point to QT cell
	public static String pointToCellQT(Point p, Node node)
	{
		if (node.getCNE() == null) // 2d
			return pointToCellQT(p.getX(), p.getY(), node);
		else // 3d
			return pointToCellQT(p.getX(), p.getY(), p.getZ(), node);
	}
	
	// point to QT cell 2d
	public static String pointToCellQT (double x, double y, Node node)
	{
		// define x, y
		double xmin = node.getXmin();
		double xmax = node.getXmax();
		double xmid = (xmin + xmax) / 2;
		
		double ymin = node.getYmin();
		double ymax = node.getYmax();
		double ymid = (ymin + ymax) / 2;
				
		if (node.getNW() != null)
		{
			if (x >= xmin && x < xmid) // point inside SW or NW
			{
				if (y >= ymin && y < ymid) // point inside SW
					return "2" + pointToCellQT(x, y, node.getSW());
				else if (y >= ymid && y < ymax) // point inside NW
					return "0" + pointToCellQT(x, y, node.getNW());
			}
			else if (x >= xmid && x < xmax) // point inside SE or NE
			{
				if (y >= ymin && y < ymid) // point inside SE
					return "3" + pointToCellQT(x, y, node.getSE());
				else if (y >= ymid && y < ymax) // point inside NE
					return "1" + pointToCellQT(x, y, node.getNE());
			}
		}
		return "";
	}
	
	// point to QT cell 3d
	public static String pointToCellQT (double x, double y, double z, Node node)
	{
		// define x, y, z
		double xmin = node.getXmin();
		double xmax = node.getXmax();
		double xmid = (xmin + xmax) / 2;
		
		double ymin = node.getYmin();
		double ymax = node.getYmax();
		double ymid = (ymin + ymax) / 2;
		
		double zmin = node.getZmin();
		double zmax = node.getZmax();
		double zmid = (zmin + zmax) / 2;
		
		if (node.getFNW() != null)
		{
			if (x >= xmin && x < xmid) // point inside SW or NW (Floor or Ceiling)
			{
				if (y >= ymin && y < ymid) // point inside SW (Floor or Ceiling)
				{
					if (z >= zmin && z < zmid) // point inside FSW
						return "2" + pointToCellQT(x, y, z, node.getFSW());
					else if (z >= zmid && z < zmax) // point inside CSW
						return "6" + pointToCellQT(x, y, z, node.getCSW());
				}
				else if (y >= ymid && y < ymax) // point inside NW (Floor or Ceiling)
				{
					if (z >= zmin && z < zmid) // point inside FNW
						return "0" + pointToCellQT(x, y, z, node.getFNW());
					else if (z >= zmid && z < zmax) // point inside CNW
						return "4" + pointToCellQT(x, y, z, node.getCNW());
				}
			}
			else if (x >= xmid && x < xmax) // point inside SE or NE (Floor or Ceiling)
			{
				if (y >= ymin && y < ymid) // point inside SE (Floor or Ceiling)
				{
					if (z >= zmin && z < zmid) // point inside FSE
						return "3" + pointToCellQT(x, y, z, node.getFSE());
					else if (z >= zmid && z < zmax) // point inside CSE
						return "7" + pointToCellQT(x, y, z, node.getCSE());
				}
				else if (y >= ymid && y < ymax) // point inside NE (Floor or Ceiling)
				{
					if (z >= zmin && z < zmid) // point inside FNE
						return "1" + pointToCellQT(x, y, z, node.getFNE());
					else if (z >= zmid && z < zmax) // point inside CNE
						return "5" + pointToCellQT(x, y, z, node.getCNE());
				}
			}
		}
		return "";
	}

	// check for duplicates in PriorityQueue
	public static boolean isDuplicate(PriorityQueue<IdDist> pq, IdDist neighbor)
	{
		for (IdDist elem : pq)
			if (elem.getId() == neighbor.getId())
				return true;

		return false;
	}
	
	// PriorityQueue<IdDist> to String
	public static String pqToString(PriorityQueue<IdDist> pq, int k, String comp)
	{
		// if we use pq directly, it will modify the original PQ, so we make a copy
		PriorityQueue<IdDist> newPQ = new PriorityQueue<>(k, new IdDistComparator(comp));

		newPQ.addAll(pq);

		StringBuilder output = new StringBuilder();

		int counter = 0;

		while (!newPQ.isEmpty() && counter < k) // add neighbors to output
		{
			IdDist elem = newPQ.poll();
			output.append(String.format("%s\t", elem));
			counter++;
		}

		return output.toString();
	}

	// merge two priority queues
	public static PriorityQueue<IdDist> joinPQ (PriorityQueue<IdDist> pq1, PriorityQueue<IdDist> pq2, int k)
	{
		PriorityQueue<IdDist> pq = new PriorityQueue<>(k, new IdDistComparator("max"));

		while (!pq1.isEmpty())
		{
			IdDist n1 = pq1.poll();
			if (!isDuplicate(pq, n1))
				pq.offer(n1);
		}

		while (!pq2.isEmpty())
		{
			IdDist n2 = pq2.poll();
			if (!isDuplicate(pq, n2))
				pq.offer(n2);
		}

		while (pq.size() > k)
			pq.poll();

		return pq;
	}

	// return point array index for point interpolation
	public static int binarySearchTpoints(double x, ArrayList<Point> points)
	{
		int low = 0;
		int high = points.size() - 1;
		int middle = (low + high + 1) / 2;
		int location = -1;
		
		do
		{
			if (x >= points.get(middle).getX())
			{
				if (middle == points.size() - 1) // middle = array length
					location = middle;
				else if (x < points.get(middle + 1).getX()) // x between middle and high
					location = middle;
				else // x greater than middle but not smaller than middle+1
					low = middle + 1;
			}
			else // x smaller than middle
				high = middle - 1;
			
			middle = (low + high + 1) / 2; // recalculate middle
			
		} while ((low < high) && (location == -1));
		
		return location;
	}
}
