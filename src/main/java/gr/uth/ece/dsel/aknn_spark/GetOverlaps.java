package gr.uth.ece.dsel.aknn_spark;

import org.apache.spark.api.java.function.Function;
import scala.Tuple2;
import gr.uth.ece.dsel.common_classes.*;
import gr.uth.ece.dsel.UtilityFunctions;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;

public final class GetOverlaps implements Function<Tuple2<String, ArrayList<Tuple2<Point, PriorityQueue<IdDist>>>>, ArrayList<Tuple2<String, Tuple2<Point, Boolean>>>>
{
	private final HashMap<String, Long> cell_tpoints; // hashmap of training points per cell list from Phase 1 <cell_id, number of training points>
	private final String partitioning; // gd or qt
	private int N; // N*N or N*N*N cells
	private final int K; // K neighbors
	private Node root; // create root node
	private final HashSet<String> overlaps; // list of overlapped cells
	private final PriorityQueue<IdDist> neighbors; // neighbors list
	private final ArrayList<Tuple2<String, Tuple2<Point, Boolean>>> outValue; // output value: <overlapped cell, qpoint, true/false>
	private boolean listComplete; // true/false if neighbors list is complete or not
	private boolean is3d = false; // false = 2d, true = 3d
	
	public GetOverlaps (HashMap<String, Long> cell_tpoints, int k, String partitioning)
	{
		this.cell_tpoints = cell_tpoints;
		this.K = k;
		this.partitioning = partitioning;
		this.overlaps = new HashSet<>();
		this.neighbors = new PriorityQueue<>(this.K, new IdDistComparator("max"));
		this.outValue = new ArrayList<>();
	}
	
	public void setN (int n)
	{
		this.N = n;
	}
	
	public void setRoot (Node root)
	{
		this.root = root;
	}
	
	@Override
	public ArrayList<Tuple2<String, Tuple2<Point, Boolean>>> call (Tuple2<String, ArrayList<Tuple2<Point, PriorityQueue<IdDist>>>> input)
	{
		final String qcell = input._1; // query points' cell
		
		// clear list
		this.outValue.clear();
		
		for (Tuple2<Point, PriorityQueue<IdDist>> tuple: input._2) // read next <query point, neighbors list> tuple
		{
			// clear lists
			this.overlaps.clear();
			this.neighbors.clear();
			//reset list status
			this.listComplete = false;
			// set 2d as default
			this.is3d = false;
			
			Point qpoint = tuple._1; // get qpoint
			this.neighbors.addAll(tuple._2); // get its neighbor list
			
			// call appropriate overlaps discovery function according to partitioning
			if (this.partitioning.equals("gd"))
				getOverlapsGD(qcell, qpoint); // find overlaps
			else if (this.partitioning.equals("qt"))
				getOverlapsQT(qcell, qpoint); // find overlaps
			
			// fill output list with tuples <cell, <qpoint, true/false>>
			for (String cell: this.overlaps)
				this.outValue.add(new Tuple2<>(cell, new Tuple2<>(qpoint, this.listComplete)));
		}
		
		return new ArrayList<>(this.outValue);
	}
	
	// find grid query overlaps
	private void getOverlapsGD (String qcell, Point qpoint)
	{
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
		
		
		How to find neighboring cells:

		If current is cell_id (x-y plane):
		 W is (cell_id - 1)
		 E is (cell_id + 1)
		 N is (cell_id + n)
		 S is (cell_id - n)
		NE is (cell_id + n + 1)
		NW is (cell_id + n - 1)
		SE is (cell_id - n + 1)
		SW is (cell_id - n - 1)
		
		On z-axis:
		above is (cell_id + n*n)
		below is (cell_id - n*n)
		
		2d:
		south row cells ( 0, 1,..., n-1 ) don't have S, SE, SW neighbors
		north row cells ( (n-1)n, (n-1)n+1,..., (n-1)n+(n-1) ) don't have N, NE, NW neighbors
		west column cells ( 0, n,..., (n-1)n ) don't have W, NW, SW neighbors
		east column cells ( n-1, 2n-1,..., (n-1)n+(n-1) ) don't have E, NE, SE neighbors
		
		
		                           xi mod ds (part of xi inside cell)
		     |----------xi------------->.(xi, yi)
		     |---------|----------|-----^--
	         |         |          |     |   yi mod ds (part of yi inside cell)
	    2*ds |---------|----------|-----|--
	         |         |          |     yi
	      ds |---------|----------|-----|--
	         |         |          |     |
	         |---------|----------|--------
	       0          ds         2*ds                  
	       

		3d:
		south wall cells ( 0 + k*N*N, 1 + k*N*N,..., n-1 + k*N*N) don't have S, SE, SW neighbors
		north wall cell's ( (n-1)n + k*N*N, (n-1)n+1 + k*N*N,..., (n-1)n+(n-1) + k*N*N ) don't have N, NE, NW neighbors
		west wall cells ( 0 + k*N*N, n + k*N*N,..., (n-1)n + k*N*N ) don't have W, NW, SW neighbors
		east wall cells ( n-1 + k*N*N, 2n-1 + k*N*N,..., (n-1)n+(n-1) + k*N*N ) don't have E, NE, SE neighbors
		k = 0, 1,..., n - 1
		
		On z-axis:
		bottom level cells ( 0, 1,..., (n-1)n+(n-1) ) don't have neighbors below
		top level cells ( (n-1)*n*n, (n-1)*n*n+1,..., (n-1)*n*n+(n-1)n+(n-1) ) don't have neighbors above
		*/
		
		// read query point coordinates and neighbors list
    	final double xq = qpoint.getX();
    	final double yq = qpoint.getY();
		final double zq = qpoint.getZ();

    	// check if 3d
    	this.is3d = (zq != Double.NEGATIVE_INFINITY);
    	
    	// find query point cell
    	final double ds = 1.0 / this.N; // interval ds (cell width)
    	final int intQCell = Integer.parseInt(qcell); // get int value of query point cell
    	
    	//final int iq = (int) (xq / ds); // get i
    	//final int jq = (int) (yq / ds); // get j
    	//final int intQCell = jq * this.N + iq; // calculate int cell_id
    	
    	// if neighbors list not empty, set circle radius R to the distance of farthest neighbor
    	// else half the cell width
		double R = !this.neighbors.isEmpty() ? this.neighbors.peek().getDist() : 0.5 * ds;
		
		// add number of training points in this cell, 0 if null
		final int tpointsInQcell = this.cell_tpoints.getOrDefault(qcell, 0L).intValue();
		
		// top-bottom rows, far left-right columns (rows and columns become walls in 3d)
		final HashSet<Integer> south_row = new HashSet<>(); // no S, SE, SW for cells in this set
		final HashSet<Integer> north_row = new HashSet<>(); // no N, NE, NW for cells in this set
		final HashSet<Integer> west_column = new HashSet<>(); // no W, NW, SW for cells in this set
		final HashSet<Integer> east_column = new HashSet<>(); // no E, NE, SE for cells in this set
		
		// top-bottom z-level cells
		final HashSet<Integer> bottom_level = new HashSet<>(); // no lower level (-n*n) for cells in this set
		final HashSet<Integer> top_level = new HashSet<>(); // no upper level (+n*n) for cells in this set
		
		final int zfloor = this.N * this.N; // N*N: when added/subtracted goes one floor up/down
		
		for (int i = 0; i < this.N; i++) // fill sets
		{
			for (int j = 0; j < this.N; j++) {
				south_row.add(i + (j * zfloor));
				north_row.add((this.N - 1) * this.N + i + (j * zfloor));
				west_column.add(i * this.N + (j * zfloor));
				east_column.add(i * this.N + this.N - 1 + (j * zfloor));

				if (this.is3d) // 3d only
				{
					bottom_level.add(j * this.N + i);
					top_level.add(j * this.N + i + (this.N - 1) * zfloor);
				}
			}
		}
		
		// case 1: there are at least knn in this cell
		boolean case1 = false;

		if (tpointsInQcell >= this.K)
		{
			double[] borders = UtilityFunctions.cellBorders(intQCell, tpointsInQcell, is3d);

			final double xmin = borders[0];
			final double xmax = borders[1];
			final double ymin = borders[2];
			final double ymax = borders[3];

			if (borders.length == 4) // 2d
			{
				if (UtilityFunctions.circleInsideCell(xq, yq, R, xmin, xmax, ymin, ymax))
					case1 = true;
			}
			else if (borders.length == 6) // 3d
			{
				final double zmin = borders[4];
				final double zmax = borders[5];

				if (UtilityFunctions.circleInsideCell(xq, yq, zq, R, xmin, xmax, ymin, ymax, zmin, zmax))
					case1 = true;
			}
			if (case1) // point goes straight to next phase
			{
				this.overlaps.add(qcell); // add qpoint cell
				this.listComplete = true; // set status to true
			}
		}
		// case 2: not enough neighbors in query point cell or circle overlaps other cells
		else
		{
			int sentinel = 0; // loop control variable

			// list of possible overlaps with circle
			final HashSet<Integer> candidateOverlaps = new HashSet<>();

			candidateOverlaps.add(intQCell);

			// dummy set of cells to be added (throws ConcurrentModificationException if trying to modify set while traversing it)
			final HashSet<Integer> tempOverlaps = new HashSet<>();

			// runs until it finds >=k tpoints, then once more
			while (sentinel < 2)
			{
				int overlaps_points = 0; // total number of training points in overlaps

				// get new layer of surrounding cells
				for (int cell : candidateOverlaps)
					tempOverlaps.addAll(UtilityFunctions.surroundingCells(cell, this.N, this.is3d, south_row, north_row, west_column, east_column, bottom_level, top_level));

				candidateOverlaps.addAll(tempOverlaps);

				tempOverlaps.clear();

				// check each cell in list if overlaps with circle/sphere
				for (int cell : candidateOverlaps)
				{
					final String strCell = String.valueOf(cell);

					// proceed only if this cell contains any training points
					// and skip query point cell
					if (cell != intQCell && this.cell_tpoints.containsKey(strCell))
					{
						// get cell's borders (xmin, xmax, ymin, ymax, zmin, zmax)
						final double[] borders = UtilityFunctions.cellBorders(cell, this.N, this.is3d);
						final double xmin = borders[0];
						final double xmax = borders[1];
						final double ymin = borders[2];
						final double ymax = borders[3];

						if (!this.is3d) // 2d
						{
							if (UtilityFunctions.circleSquareIntersect(xq, yq, R, xmin, xmax, ymin, ymax))
								this.overlaps.add(strCell);
						}
						else // 3d
						{
							final double zmin = borders[4];
							final double zmax = borders[5];

							if (UtilityFunctions.sphereCubeIntersect(xq, yq, zq, R, xmin, xmax, ymin, ymax, zmin, zmax))
								this.overlaps.add(strCell);
						}
					}
				}

				// remove query point cell (because its training points have been already counted)
				// also remove empty cells
				this.overlaps.removeIf(cell -> (cell.equals(qcell) || !this.cell_tpoints.containsKey(cell)));

				// count total training points from overlaps
				if (!this.overlaps.isEmpty())
					for (String cell : this.overlaps)
						overlaps_points += this.cell_tpoints.get(cell);

				R += 0.5 * ds; // increase radius by half ds

				// if k neighbors found, run loop one more time and increase radius by the diagonal of a cell
				if (overlaps_points + tpointsInQcell >= this.K)
				{
					sentinel++;

					if (!this.is3d) // 2d pythagorean
						R += Math.sqrt(2) * ds;
					else // 3d pythagorean
						R += Math.sqrt(3) * ds;
				}
			}
		}
	} // end getOverlapsGD
	
	// find quadtree query overlaps
	private void getOverlapsQT (String qcell, Point qpoint)
	{
		// read query point coordinates and neighbors list
    	final double xq = qpoint.getX();
    	final double yq = qpoint.getY();
    	final double zq = qpoint.getZ();

		// check if 3d
		this.is3d = (zq != Double.NEGATIVE_INFINITY);
    	
		/* If
		 * root cell side length = L
		 * and for example
		 * cell id = 3012 (4 digits)
		 * then cell's length = L / (2 ^ 4)
		 */
    	
    	// total number of training points in this cell, 0 if null
    	final int tpointsInQcell = this.cell_tpoints.getOrDefault(qcell, 0L).intValue();
    	
    	final double ds = 1.0 / Math.pow(2, qcell.length()); // ds = query cell width
		
    	// if neighbors list not empty, set circle/sphere radius R to the distance of farthest neighbor
    	// else half the cell width
		double R = !this.neighbors.isEmpty() ? this.neighbors.peek().getDist() : 0.5 * ds;
		
		// case 1: there are at least knn in this cell
		boolean case1 = false;

		if (tpointsInQcell >= this.K)
		{
			double[] borders = UtilityFunctions.cellBorders(qcell);

			final double xmin = borders[0];
			final double xmax = borders[1];
			final double ymin = borders[2];
			final double ymax = borders[3];

			if (borders.length == 4) // 2d
			{
				if (UtilityFunctions.circleInsideCell(xq, yq, R, xmin, xmax, ymin, ymax))
					case1 = true;
			}
			else if (borders.length == 6) // 3d
			{
				final double zmin = borders[4];
				final double zmax = borders[5];

				if (UtilityFunctions.circleInsideCell(xq, yq, zq, R, xmin, xmax, ymin, ymax, zmin, zmax))
					case1 = true;
			}
			if (case1) // point goes straight to next phase
			{
				this.overlaps.add(qcell); // add qpoint cell
				this.listComplete = true; // set status to true
			}
		}
		
		// case 2: there are less than knn in this cell
		else
		{
			/* Define a new increasing radius r1:
			 * - if there are already x < k neighbors in this cell, we suppose a constant density of training points, so
			 *   (2d) in a circle of radius r and area pi*r*r we get x neighbors
			 *        in a circle of radius r1 and area pi*r1*r1 we will get k neighbors
			 *        so r1 = sqrt(k/x)*r
			 *   (3d) in a sphere of radius r and volume 4/3*pi*r*r*r we get x neighbors
			 *        in a sphere of radius r1 and volume 4/3*pi*r1*r1*r1 we will get k neighbors
			 *        so r1 = cubic_root(k/x)*r
			 * - if no neighbors are found in this cell, then we suppose that in radius r = ds/2 we have found only one neighbor (x = 1)
			 *   and so 
			 *   (2d) r1 = 1.2*sqrt(k)*r
			 *   (3d) r1 = 1.2*cubic_root(k)*r
			 *   (we gave it an additional 20% boost)
			 */
			double r1 = 0;

			int n = this.neighbors.size() / 2; // divide by 2 because knnlist also contains distances (size = 2*[number of neighbors])

			// if x > 0 (there are some neighbors in the list) set first value
			// else (no neighbors) set second value
			
			// 2d
			if (!this.is3d)
				r1 = (n > 0) ? Math.sqrt((double) this.K / n) * R : 1.2 * Math.sqrt(this.K) * R;
			// 3d
			else
				r1 = (n > 0) ? Math.cbrt((double) this.K / n) * R : 1.2 * Math.cbrt(this.K) * R;
			
			int sentinel = 0; // loop control variable

			while (sentinel < 2) // trying to find overlaps to fill k-nn
			{
				int overlaps_points = 0; // total number of training points in overlaps

				this.overlaps.clear(); // clear overlaps list
				
				// draw circle/sphere and check for overlaps
				// 2d
				if (!this.is3d)
					rangeQuery(xq, yq, r1, this.root, "");
				// 3d
				else
					rangeQuery(xq, yq, zq, r1, this.root, "");
				
				// remove query point cell (because its training points have been already counted)
				// also remove empty cells
				this.overlaps.removeIf(cell -> (cell.equals(qcell) || !this.cell_tpoints.containsKey(cell)));

				// count total training points from overlaps
				if (!this.overlaps.isEmpty())
					for (String cell : this.overlaps)
						overlaps_points += this.cell_tpoints.get(cell);
				
				r1 += 0.1 * r1; // increase radius by 10%
				
				// if k neighbors found (first time only), run loop one more time and set r1 equal to the maximum distance from ipoint to all overlapped cells
				if ((overlaps_points + tpointsInQcell >= this.K) && (sentinel == 0))
				{
					sentinel++;

					double maxSqrDist = 0; // square of maximum distance found so far (ipoint to cell)

					for (String cell : this.overlaps) // for every cell in overlaps
					{
						// get cell's borders
						final double[] borders = UtilityFunctions.cellBorders(cell);
						final double xmin = borders[0];
						final double xmax = borders[1];
						final double ymin = borders[2];
						final double ymax = borders[3];

						// maximum ipoint distance from cell in x-direction
						final double maxX = Math.max(Math.abs(xq - xmin), Math.abs(xq - xmax));
						// maximum ipoint distance from cell in y-direction
						final double maxY = Math.max(Math.abs(yq - ymin), Math.abs(yq - ymax));

						if (!this.is3d) // 2d
						{
							// replace current maximum squared distance if new is bigger
							maxSqrDist = Math.max(maxSqrDist, maxX * maxX + maxY * maxY);
						}
						else // 3d
						{
							final double zmin = borders[4];
							final double zmax = borders[5];

							// maximum ipoint distance from cell in z-direction
							final double maxZ = Math.max(Math.abs(zq - zmin), Math.abs(zq - zmax));

							// replace current maximum squared distance if new is bigger
							maxSqrDist = Math.max(maxSqrDist, maxX * maxX + maxY * maxY + maxZ * maxZ);
						}
					} // end for
					r1 = Math.sqrt(maxSqrDist); // run overlaps check once more with this radius
				} // end if
				else if (sentinel == 1)
					sentinel++;
			} // end while
		} // end else
	} // end getOverlapsQT
	
	// 2d quadtree range query
	private void rangeQuery (double x, double y, double r, Node node, String address)
	{
		// leaf node
		if (node.getNW() == null)
			this.overlaps.add(address);
		
		// internal node
		else
		{
			if (UtilityFunctions.intersect(x, y, r, node.getNW()))
				rangeQuery(x, y, r, node.getNW(), address + "0");
			
			if (UtilityFunctions.intersect(x, y, r, node.getNE()))
				rangeQuery(x, y, r, node.getNE(), address + "1");
			
			if (UtilityFunctions.intersect(x, y, r, node.getSW()))
				rangeQuery(x, y, r, node.getSW(), address + "2");
			
			if (UtilityFunctions.intersect(x, y, r, node.getSE()))
				rangeQuery(x, y, r, node.getSE(), address + "3");
		}
	}
	
	// 3d quadtree range query
	private void rangeQuery(double x, double y, double z, double r, Node node, String address)
	{
		// leaf node
		if (node.getFNW() == null)
			this.overlaps.add(address);

		// internal node
		else
		{
			if (UtilityFunctions.intersect(x, y, z, r, node.getFNW()))
				rangeQuery(x, y, z, r, node.getFNW(), address + "0");
			
			if (UtilityFunctions.intersect(x, y, z, r, node.getFNE()))
				rangeQuery(x, y, z, r, node.getFNE(), address + "1");
			
			if (UtilityFunctions.intersect(x, y, z, r, node.getFSW()))
				rangeQuery(x, y, z, r, node.getFSW(), address + "2");
			
			if (UtilityFunctions.intersect(x, y, z, r, node.getFSE()))
				rangeQuery(x, y, z, r, node.getFSE(), address + "3");
			
			if (UtilityFunctions.intersect(x, y, z, r, node.getCNW()))
				rangeQuery(x, y, z, r, node.getCNW(), address + "4");
			
			if (UtilityFunctions.intersect(x, y, z, r, node.getCNE()))
				rangeQuery(x, y, z, r, node.getCNE(), address + "5");
			
			if (UtilityFunctions.intersect(x, y, z, r, node.getCSW()))
				rangeQuery(x, y, z, r, node.getCSW(), address + "6");
			
			if (UtilityFunctions.intersect(x, y, z, r, node.getCSE()))
				rangeQuery(x, y, z, r, node.getCSE(), address + "7");
		}
	}
}
