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
	}
	
	// return euclidean distance between two 3d points (x1, y1, z1) and (x2, y2, z2)
	public static double distance (double x1, double y1, double z1, double x2, double y2, double z2)
	{
		return Math.sqrt(square_distance (x1, y1, z1, x2, y2, z2));
	}
	
	// return Euclidean distance between two points
	public static double distance (Point ipoint, Point tpoint)
	{
		// 3d points
		if (ipoint.getZ() != Double.NEGATIVE_INFINITY && tpoint.getZ() != Double.NEGATIVE_INFINITY)
			return distance(ipoint.getX(), ipoint.getY(), ipoint.getZ(), tpoint.getX(), tpoint.getY(), tpoint.getZ());
		else // 2d points
			return distance(ipoint.getX(), ipoint.getY(), tpoint.getX(), tpoint.getY());
	}
	
	// return square of euclidean distance between two 2d points (x1, y1) and (x2, y2)
	public static double square_distance (double x1, double y1, double x2, double y2)
	{
		return Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2);
	}
	
	// return square of Euclidean distance between two 3d points (x1, y1, z1) and (x2, y2, z2)
	public static double square_distance (double x1, double y1, double z1, double x2, double y2, double z2)
	{
		return Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2) + Math.pow(z1 - z2,  2);
	}
	
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
	
	// point to GD cell
	public static String pointToCell(Point p, int n)
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
//		return pointToCell((node.getXmin() + node.getXmax()) / 2, (node.getYmin() + node.getYmax()) / 2, node);
//	}
	
	// point to QT cell
	public static String pointToCell(Point p, Node node)
	{
		if (node.getCNE() == null) // 2d
			return pointToCell(p.getX(), p.getY(), node);
		else // 3d
			return pointToCell(p.getX(), p.getY(), p.getZ(), node);
	}
	
	// point to QT cell 2d
	public static String pointToCell (double x, double y, Node node)
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
					return "2" + pointToCell(x, y, node.getSW());
				else if (y >= ymid && y < ymax) // point inside NW
					return "0" + pointToCell(x, y, node.getNW());
			}
			else if (x >= xmid && x < xmax) // point inside SE or NE
			{
				if (y >= ymin && y < ymid) // point inside SE
					return "3" + pointToCell(x, y, node.getSE());
				else if (y >= ymid && y < ymax) // point inside NE
					return "1" + pointToCell(x, y, node.getNE());
			}
		}
		return "";
	}
	
	// point to QT cell 3d
	public static String pointToCell (double x, double y, double z, Node node)
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
						return "2" + pointToCell(x, y, z, node.getFSW());
					else if (z >= zmid && z < zmax) // point inside CSW
						return "6" + pointToCell(x, y, z, node.getCSW());
				}
				else if (y >= ymid && y < ymax) // point inside NW (Floor or Ceiling)
				{
					if (z >= zmin && z < zmid) // point inside FNW
						return "0" + pointToCell(x, y, z, node.getFNW());
					else if (z >= zmid && z < zmax) // point inside CNW
						return "4" + pointToCell(x, y, z, node.getCNW());
				}
			}
			else if (x >= xmid && x < xmax) // point inside SE or NE (Floor or Ceiling)
			{
				if (y >= ymin && y < ymid) // point inside SE (Floor or Ceiling)
				{
					if (z >= zmin && z < zmid) // point inside FSE
						return "3" + pointToCell(x, y, z, node.getFSE());
					else if (z >= zmid && z < zmax) // point inside CSE
						return "7" + pointToCell(x, y, z, node.getCSE());
				}
				else if (y >= ymid && y < ymax) // point inside NE (Floor or Ceiling)
				{
					if (z >= zmin && z < zmid) // point inside FNE
						return "1" + pointToCell(x, y, z, node.getFNE());
					else if (z >= zmid && z < zmax) // point inside CNE
						return "5" + pointToCell(x, y, z, node.getCNE());
				}
			}
		}
		return "";
	}

	// 2d quadtree intersect
	public static boolean intersect (double x, double y, double r, Node node)
	{
		return circleSquareIntersect(x, y, r, node.getXmin(), node.getXmax(), node.getYmin(), node.getYmax());
	}

	// 3d quadtree intersect
	public static boolean intersect (double x, double y, double z, double r, Node node)
	{
		return sphereCubeIntersect(x, y, z, r, node.getXmin(), node.getXmax(), node.getYmin(), node.getYmax(), node.getZmin(), node.getZmax());
	}

	// 2d circle - square intersection check
	public static boolean circleSquareIntersect (double x, double y, double r, double xmin, double xmax, double ymin, double ymax)
	{
		// if point is inside cell return true
		if (x >= xmin && x <= xmax && y >= ymin && y <= ymax)
			return true;

		// check circle - cell collision
		final double ds = xmax - xmin; // cell's width

		// get cell center coordinates
		final double xc = (xmin + xmax) / 2;
		final double yc = (ymin + ymax) / 2;

		// circle center to cell center distance
		final double centers_dist_x = Math.abs(x - xc);
		final double centers_dist_y = Math.abs(y - yc);

		// if centers are far in either direction, return false
		if (centers_dist_x > r + ds / 2)
			return false;
		if (centers_dist_y > r + ds / 2)
			return false;

		// if control reaches here, centers are close enough

		// the next two cases mean that circle center is within a stripe of width r around the square
		if (centers_dist_x < ds / 2)
			return true;
		if (centers_dist_y < ds / 2)
			return true;

		// else check the corner distance
		final double corner_dist_sq = UtilityFunctions.square_distance(centers_dist_x, centers_dist_y, ds / 2, ds / 2);

		return corner_dist_sq <= r * r;
	}

	// 3d sphere - cube intersection check
	public static boolean sphereCubeIntersect (double x, double y, double z, double r, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
	{
		// if point is inside cell return true
		if (x >= xmin && x <= xmax && y >= ymin && y <= ymax && z >= zmin && z <= zmax)
			return true;

		// check sphere - cell collision
		final double ds = xmax - xmin; // cell's width
		// get cell center coordinates
		final double xc = (xmin + xmax) / 2;
		final double yc = (ymin + ymax) / 2;
		final double zc = (zmin + zmax) / 2;
		// sphere center to cell center distance
		final double centers_dist_x = Math.abs(x - xc);
		final double centers_dist_y = Math.abs(y - yc);
		final double centers_dist_z = Math.abs(z - zc);

		// if centers are far in either direction, return false
		if (centers_dist_x > r + ds / 2)
			return false;
		if (centers_dist_y > r + ds / 2)
			return false;
		if (centers_dist_z > r + ds / 2)
			return false;

		// if control reaches here, centers are close enough

		// the next three cases mean that sphere center is within a stripe of width r around the cell
		if (centers_dist_x < ds / 2)
			return true;
		if (centers_dist_y < ds / 2)
			return true;
		if (centers_dist_z < ds / 2)
			return true;

		// else check the corner distance
		final double corner_dist_sq = UtilityFunctions.square_distance(centers_dist_x, centers_dist_y, centers_dist_z, ds / 2, ds / 2, ds / 2);

		return corner_dist_sq <= r * r;
	}

	// get GD cell's i, j, k
	public static int[] cellIJK (int cell, int N, boolean is3d)
	{
		final int i = cell % N; // 2d/3d
		final int j = (cell - i) / N; // 2d
		
		if (is3d)
		{
			final int j3d = j % N;
			final int k = ((cell - i) / N - j) / N;

			return new int[]{i, j3d, k};
		}

		return new int[]{i, j};
	}

	// get min & max x, y, z of GD cell in integer form
	public static double[] cellBorders (int cell, int N, boolean is3d)
	{
		// get cell's i, j, k
		final int[] ijk = cellIJK(cell, N, is3d);
		final int i = ijk[0];
		final int j = ijk[1];

		final double ds = 1.0 / N;

		final double xmin = i * ds;
		final double xmax = (i + 1) * ds;
		final double ymin = j * ds;
		final double ymax = (j + 1) * ds;

		if (is3d)
		{
			final int k = ijk[2];

			final double zmin = k * ds;
			final double zmax = (k + 1) * ds;

			return new double[]{xmin, xmax, ymin, ymax, zmin, zmax};
		}

		return new double[]{xmin, xmax, ymin, ymax};
	}

	// get min & max x, y, z of QT cell in string form
	public static double[] cellBorders (String cell)
	{

		double xmin = 0; // cell's floor-south-west corner coords initialization
		double ymin = 0;
		double zmin = 0;

		for (int i = 0; i < cell.length(); i++) // check cellname's digits
		{
			switch(cell.charAt(i))
			{
				case '0': // 2d / 3d
					ymin += 1.0 / Math.pow(2, i + 1); // if digit = 0 increase y0
					break;
				case '1': // 2d / 3d
					xmin += 1.0 / Math.pow(2, i + 1); // if digit = 1 increase x0
					ymin += 1.0 / Math.pow(2, i + 1); // and y0
					break;
				case '3': // 2d / 3d
					xmin += 1.0 / Math.pow(2, i + 1); // if digit = 3 increase x0
					break;
				case '4': // 3d only
					ymin += 1.0 / Math.pow(2, i + 1); // if digit = 4 increase y0
					zmin += 1.0 / Math.pow(2, i + 1); // and z0
					break;
				case '5': // 3d only
					xmin += 1.0 / Math.pow(2, i + 1); // if digit = 5 increase x0
					ymin += 1.0 / Math.pow(2, i + 1); // and y0
					zmin += 1.0 / Math.pow(2, i + 1); // and z0
					break;
				case '6':
					zmin += 1.0 / Math.pow(2, i + 1); // if digit = 6 increase z0
					break;
				case '7': // 3d only
					xmin += 1.0 / Math.pow(2, i + 1); // if digit = 7 increase x0
					zmin += 1.0 / Math.pow(2, i + 1); // and z0
					break;
			}
		}

		final double ds = 1.0 / Math.pow(2, cell.length()); // cell side length

		final double xmax = xmin + ds;
		final double ymax = ymin + ds;
		final double zmax = zmin + ds;

		/*
		 * cell's lower left corner: xmin, ymin
		 *        upper left corner: xmin, ymax
		 *        upper right corner: xmax, ymax
		 *        lower right corner: xmax, ymin
		 *
		 */

		return new double[]{xmin, xmax, ymin, ymax, zmin, zmax};
	}

	// return true if circle (x, y, r) is completely inside 2d cell (xmin, xmax, ymin, ymax)
	public static boolean circleInsideCell (double x, double y, double r, double xmin, double xmax, double ymin, double ymax)
	{
		// if r > distance of center to cell borders, circle is not completely inside cell
		if (r > Math.abs(x - xmin))
			return false;
		if (r > Math.abs(x - xmax))
			return false;
		if (r > Math.abs(y - ymin))
			return false;
		if (r > Math.abs(y - ymax))
			return false;

		return true;
	}

	// return true if sphere (x, y, z, r) is completely inside 3d cell (xmin, xmax, ymin, ymax, zmin, zmax)
	public static boolean circleInsideCell (double x, double y, double z, double r, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
	{
		// if r > distance of center to cell borders, sphere is not completely inside cell
		if (r > Math.abs(x - xmin))
			return false;
		if (r > Math.abs(x - xmax))
			return false;
		if (r > Math.abs(y - ymin))
			return false;
		if (r > Math.abs(y - ymax))
			return false;
		if (r > Math.abs(z - zmin))
			return false;
		if (r > Math.abs(z - zmax))
			return false;

		return true;
	}

	// return surrounding cells of GD cell in integer form
	public static HashSet<Integer> surroundingCells (int cell, int N, boolean is3d, HashSet<Integer> southRow, HashSet<Integer> northRow, HashSet<Integer> westColumn, HashSet<Integer> eastColumn, HashSet<Integer> bottomLevel, HashSet<Integer> topLevel)
	{
		final HashSet<Integer> surCells = new HashSet<>();

		if (!westColumn.contains(cell)) // excluding west column
			surCells.add(cell - 1); // W
		if (!eastColumn.contains(cell)) // excluding east column
			surCells.add(cell + 1); // E
		if (!southRow.contains(cell)) // excluding southRow
			surCells.add(cell - N); // S
		if (!northRow.contains(cell)) // excluding northRow
			surCells.add(cell + N); // N
		if (!southRow.contains(cell) && !westColumn.contains(cell)) // excluding south row and west column
			surCells.add(cell - N - 1); // SW
		if (!southRow.contains(cell) && !eastColumn.contains(cell)) // excluding south row and east column
			surCells.add(cell - N + 1); // SE
		if (!northRow.contains(cell) && !westColumn.contains(cell)) // excluding north row and west column
			surCells.add(cell + N - 1); // NW
		if (!northRow.contains(cell) && !eastColumn.contains(cell)) // excluding north row and east column
			surCells.add(cell + N + 1); // NE

		if (is3d) // 3d
		{
			final int zfloor = N * N;

			if (!topLevel.contains(cell)) // excluding top level
				surCells.add(cell + zfloor); // above
			if (!topLevel.contains(cell) && !westColumn.contains(cell)) // excluding top level & west wall
				surCells.add(cell + zfloor - 1); // above-west
			if (!topLevel.contains(cell) && !eastColumn.contains(cell)) // excluding top level & east wall
				surCells.add(cell + zfloor + 1);// above-east
			if (!topLevel.contains(cell) && !northRow.contains(cell)) // excluding top level & north wall
				surCells.add(cell + zfloor + N); // above-N
			if (!topLevel.contains(cell) && !southRow.contains(cell)) // excluding top level & south wall
				surCells.add(cell + zfloor - N); // above-S
			if (!topLevel.contains(cell) && !southRow.contains(cell) && !westColumn.contains(cell)) // excluding top level & south wall & west wall
				surCells.add(cell + zfloor - N - 1); // above-SW
			if (!topLevel.contains(cell) && !southRow.contains(cell) && !eastColumn.contains(cell)) // excluding top level & south wall & east wall
				surCells.add(cell + zfloor - N + 1); // above-SE
			if (!topLevel.contains(cell) && !northRow.contains(cell) && !westColumn.contains(cell)) // excluding top level & north wall & west wall
				surCells.add(cell + zfloor + N - 1); // above-NW
			if (!topLevel.contains(cell) && !northRow.contains(cell) && !eastColumn.contains(cell)) // excluding top level & north wall & east wall
				surCells.add(cell + zfloor + N + 1); // above-NE
			if (!bottomLevel.contains(cell)) // excluding bottom level
				surCells.add(cell - zfloor); // below
			if (!bottomLevel.contains(cell) && !westColumn.contains(cell)) // excluding bottom level & west wall
				surCells.add(cell - zfloor - 1); // below-W
			if (!bottomLevel.contains(cell) && !eastColumn.contains(cell)) // excluding bottom level & east wall
				surCells.add(cell - zfloor + 1); // below-E
			if (!bottomLevel.contains(cell) && !northRow.contains(cell)) // excluding bottom level & north wall
				surCells.add(cell - zfloor + N); // below-N
			if (!bottomLevel.contains(cell) && !southRow.contains(cell)) // excluding bottom level & south wall
				surCells.add(cell - zfloor - N); // below-S
			if (!bottomLevel.contains(cell) && !southRow.contains(cell) && !westColumn.contains(cell)) // excluding bottom level & south wall & west wall
				surCells.add(cell - zfloor - N - 1); // below-SW
			if (!bottomLevel.contains(cell) && !southRow.contains(cell) && !eastColumn.contains(cell)) // excluding bottom level & south wall & east wall
				surCells.add(cell - zfloor - N + 1); // below-SE
			if (!bottomLevel.contains(cell) && !northRow.contains(cell) && !westColumn.contains(cell)) // excluding bottom level & north wall & west wall
				surCells.add(cell - zfloor + N - 1); // below-NW
			if (!bottomLevel.contains(cell) && !northRow.contains(cell) && !eastColumn.contains(cell)) // excluding bottom level & north wall & east wall
				surCells.add(cell - zfloor + N + 1); // below-NE
		}

		return surCells;
	}

	// check for duplicates in PriorityQueue
	public static boolean isDuplicate (PriorityQueue<IdDist> pq, IdDist neighbor)
	{		
		return pq.stream().anyMatch(el -> el.getId() == neighbor.getId());
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
