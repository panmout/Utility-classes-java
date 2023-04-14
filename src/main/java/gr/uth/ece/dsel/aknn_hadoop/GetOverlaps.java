package gr.uth.ece.dsel.aknn_hadoop;

import gr.uth.ece.dsel.common_classes.IdDist;
import gr.uth.ece.dsel.common_classes.Node;
import gr.uth.ece.dsel.common_classes.Point;

import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;

public final class GetOverlaps
{
	private final HashMap<String, Integer> cell_tpoints; // hashmap of training points per cell list from Phase 1 <cell_id, number of training points>
	private final String partitioning; // gd or qt
	private int N; // (2d) N*N or (3d) N*N*N cells
	private final int K; // AKNN K
	private Node root; // create root node
	private final HashSet<String> overlaps;
	private String qcell;
	private Point qpoint;
	private PriorityQueue<IdDist> neighbors;
	private boolean is3d = false; // false = 2d, true = 3d

	public GetOverlaps (HashMap<String, Integer> cell_tpoints, int k, String partitioning) {
		this.cell_tpoints = new HashMap<>(cell_tpoints);
		this.K = k;
		this.partitioning = partitioning;
		this.overlaps = new HashSet<>();
	}

	public void initializeFields (Point qp, String qc, PriorityQueue<IdDist> phase2neighbors) {
		this.qpoint = qp;
		this.qcell = qc;
		this.neighbors = new PriorityQueue<>(phase2neighbors);
	}

	public void setN (int n) {
		this.N = n;
	}

	public void setRoot (Node root) {
		this.root = root;
	}

	public HashSet<String> getOverlaps() {
		this.overlaps.clear();

		// set 2d as default
		this.is3d = false;

		// call appropriate function according to partitioning and iterate query points of this cell
		if (this.partitioning.equals("gd"))
			getOverlapsGD(this.qcell, this.qpoint); // find its overlaps
		else if (this.partitioning.equals("qt"))
			getOverlapsQT(this.qcell, this.qpoint); // find its overlaps

		return this.overlaps;
	}

	// find grid overlaps
	private void getOverlapsGD (String qcell, Point qpoint) {

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
		final int tpointsInQcell = this.cell_tpoints.getOrDefault(qcell, 0);

		// top-bottom rows, far left-right columns (rows and columns become walls in 3d)
		final HashSet<Integer> south_row = new HashSet<>(); // no S, SE, SW for cells in this set
		final HashSet<Integer> north_row = new HashSet<>(); // no N, NE, NW for cells in this set
		final HashSet<Integer> west_column = new HashSet<>(); // no W, NW, SW for cells in this set
		final HashSet<Integer> east_column = new HashSet<>(); // no E, NE, SE for cells in this set

		// top-bottom z-level cells
		final HashSet<Integer> bottom_level = new HashSet<>(); // no lower level (-n*n) for cells in this set
		final HashSet<Integer> top_level = new HashSet<>(); // no upper level (+n*n) for cells in this set

		final int zfloor = this.is3d ? this.N * this.N : 0; // 0 for 2d and N*N for 3d: when added/subtracted goes one floor up/down

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

		// case 1: there are at least knn in this cell and circle/sphere with radius R is completely inside the cell
		if (tpointsInQcell >= this.K && circleContainedInCell(xq, yq, zq, R, intQCell)) {
			// point goes straight to next phase
			this.overlaps.add(qcell); // add qpoint cell
		} else // case 2: not enough neighbors in query point cell or circle overlaps other cells
		{
			//int overlaps_points = tpointsInQcell; // total number of training points in overlaps

			int sentinel = 0; // loop control variable

			// list of possible overlaps with circle
			final HashSet<Integer> candidateOverlaps = new HashSet<>();

			candidateOverlaps.add(intQCell);

			// dummy set of cells to be added (throws ConcurrentModificationException if trying to modify set while traversing it)
			final HashSet<Integer> tempOverlaps = new HashSet<>();

			// runs until it finds >=k tpoints, then once more
			while (sentinel < 2) {
				int overlaps_points = 0; // total number of training points in overlaps

				// get new layer of surrounding cells
				for (int cell : candidateOverlaps)
					tempOverlaps.addAll(surroundingCells(cell, south_row, north_row, west_column, east_column, bottom_level, top_level));

				candidateOverlaps.addAll(tempOverlaps);

				tempOverlaps.clear();

				// check each cell in list if overlaps with circle/sphere
				for (int cell : candidateOverlaps) {
					final String strCell = String.valueOf(cell);

					// proceed only if this cell contains any training points
					if (this.cell_tpoints.containsKey(strCell)) {
						// get cell's borders (xmin, xmax, ymin, ymax, zmin, zmax)
						final double xmin = cellBorders(cell)[0];
						final double xmax = cellBorders(cell)[1];
						final double ymin = cellBorders(cell)[2];
						final double ymax = cellBorders(cell)[3];

						if (!this.is3d) // 2d
						{
							if (circleSquareIntersect(xq, yq, R, xmin, xmax, ymin, ymax))
								this.overlaps.add(strCell);
						} else // 3d
						{
							final double zmin = cellBorders(cell)[4];
							final double zmax = cellBorders(cell)[5];

							if (sphereCubeIntersect(xq, yq, zq, R, xmin, xmax, ymin, ymax, zmin, zmax))
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
				if (overlaps_points + tpointsInQcell >= this.K) {
					sentinel++;

					if (!this.is3d) // 2d pythagorean
						R += Math.sqrt(2) * ds;
					else // 3d pythagorean
						R += Math.sqrt(3) * ds;
				}
			}
		}
	} // end getOverlapsGD

	// find quadtree overlaps
	private void getOverlapsQT (String qcell, Point qpoint) {
		// read query point coordinates and neighbors list
		final double xq = qpoint.getX();
		final double yq = qpoint.getY();
		final double zq = qpoint.getZ();

		// check if 3d
		this.is3d = (zq != Double.NEGATIVE_INFINITY);

		/* If
		 * root cell's side length = L
		 * and cell id has n digits
		 * then cell's side length = L / (2 ^ n)
		 */

		// total number of training points in this cell, 0 if null
		final int tpointsInQcell = this.cell_tpoints.getOrDefault(qcell, 0);

		final double ds = 1.0 / Math.pow(2, qcell.length()); // ds = query cell width

		// if neighbors list not empty, set circle/sphere radius R to the distance of farthest neighbor
		// else half the cell width
		double R = !this.neighbors.isEmpty() ? this.neighbors.peek().getDist() : 0.5 * ds;

		// case 1: there are at least knn in this cell and circle/sphere with radius R is completely inside the cell
		if (tpointsInQcell >= this.K && circleContainedInCell(xq, yq, zq, R, qcell)) {
			// point goes straight to next phase
			this.overlaps.add(qcell); // add qpoint cell
		} else // case 2: not enough neighbors in query point cell or circle overlaps other cells
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

			double r1;

			int n = this.neighbors.size() / 2; // divide by 2 because knnlist also contains distances (size = 2*[number of neighbors])

			// if x > 0 (there are some neighbors in the list) set first value
			// else (no neighbors) set second value

			if (!this.is3d) // 2d
				r1 = (n > 0) ? Math.sqrt((double) this.K / n) * R : 1.2 * Math.sqrt(this.K) * R;
			else // 3d
				r1 = (n > 0) ? Math.cbrt((double) this.K / n) * R : 1.2 * Math.cbrt(this.K) * R;

			int sentinel = 0; // loop control variable

			// runs until it finds >=k tpoints, then once more
			while (sentinel < 2)
			{
				int overlaps_points = 0; // total number of training points in overlaps

				// draw circle/sphere and check for overlaps
				if (!this.is3d) // 2d
					rangeQuery(xq, yq, r1, this.root, "");
				else // 3d
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
				if ((overlaps_points + tpointsInQcell >= this.K) && (sentinel == 0)) {
					sentinel++;

					double maxSqrDist = 0; // square of maximum distance found so far (ipoint to cell)

					for (String cell : this.overlaps) // for every cell in overlaps
					{
						// get cell's borders
						final double xmin = cellBorders(cell)[0];
						final double xmax = cellBorders(cell)[1];
						final double ymin = cellBorders(cell)[2];
						final double ymax = cellBorders(cell)[3];
						final double zmin = cellBorders(cell)[4];
						final double zmax = cellBorders(cell)[5];

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
	private void rangeQuery (double x, double y, double r, Node node, String address) {
		if (node.getNW() == null) // leaf node
			this.overlaps.add(address);

			// internal node
		else {
			if (intersect(x, y, r, node.getNW()))
				rangeQuery(x, y, r, node.getNW(), address + "0");

			if (intersect(x, y, r, node.getNE()))
				rangeQuery(x, y, r, node.getNE(), address + "1");

			if (intersect(x, y, r, node.getSW()))
				rangeQuery(x, y, r, node.getSW(), address + "2");

			if (intersect(x, y, r, node.getSE()))
				rangeQuery(x, y, r, node.getSE(), address + "3");
		}
	}

	// 2d quadtree intersect
	private boolean intersect (double x, double y, double r, Node node) {
		return circleSquareIntersect(x, y, r, node.getXmin(), node.getXmax(), node.getYmin(), node.getYmax());
	}

	// 3d quadtree range query
	private void rangeQuery (double x, double y, double z, double r, Node node, String address) {
		if (node.getFNW() == null) // leaf node
			this.overlaps.add(address);
			// internal node
		else {
			if (intersect(x, y, z, r, node.getFNW()))
				rangeQuery(x, y, z, r, node.getFNW(), address + "0");

			if (intersect(x, y, z, r, node.getFNE()))
				rangeQuery(x, y, z, r, node.getFNE(), address + "1");

			if (intersect(x, y, z, r, node.getFSW()))
				rangeQuery(x, y, z, r, node.getFSW(), address + "2");

			if (intersect(x, y, z, r, node.getFSE()))
				rangeQuery(x, y, z, r, node.getFSE(), address + "3");

			if (intersect(x, y, z, r, node.getCNW()))
				rangeQuery(x, y, z, r, node.getCNW(), address + "4");

			if (intersect(x, y, z, r, node.getCNE()))
				rangeQuery(x, y, z, r, node.getCNE(), address + "5");

			if (intersect(x, y, z, r, node.getCSW()))
				rangeQuery(x, y, z, r, node.getCSW(), address + "6");

			if (intersect(x, y, z, r, node.getCSE()))
				rangeQuery(x, y, z, r, node.getCSE(), address + "7");
		}
	}

	// 3d quadtree intersect
	private boolean intersect (double x, double y, double z, double r, Node node) {
		return sphereCubeIntersect(x, y, z, r, node.getXmin(), node.getXmax(), node.getYmin(), node.getYmax(), node.getZmin(), node.getZmax());
	}

	// 2d circle - square intersection check
	private boolean circleSquareIntersect (double x, double y, double r, double xmin, double xmax, double ymin, double ymax) {
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
		final double corner_dist_sq = (centers_dist_x - ds / 2)*(centers_dist_x - ds / 2) + (centers_dist_y - ds / 2)*(centers_dist_y - ds / 2);

		return corner_dist_sq <= r * r;
	}

	// 3d sphere - cube intersection check
	private boolean sphereCubeIntersect (double x, double y, double z, double r, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
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
		final double corner_dist_sq = Math.pow(centers_dist_x - ds / 2, 2) + Math.pow(centers_dist_y - ds / 2, 2) + Math.pow(centers_dist_z - ds / 2, 2);

		return corner_dist_sq <= r * r;
	}

	// return surrounding cells of GD cell in integer form
	private HashSet<Integer> surroundingCells (int cell, HashSet<Integer> southRow, HashSet<Integer> northRow, HashSet<Integer> westColumn, HashSet<Integer> eastColumn, HashSet<Integer> bottomLevel, HashSet<Integer> topLevel) {
		final HashSet<Integer> surCells = new HashSet<>();

		if (!westColumn.contains(cell)) // excluding west column
			surCells.add(cell - 1); // W
		if (!eastColumn.contains(cell)) // excluding east column
			surCells.add(cell + 1); // E
		if (!southRow.contains(cell)) // excluding southRow
			surCells.add(cell - this.N); // S
		if (!northRow.contains(cell)) // excluding northRow
			surCells.add(cell + this.N); // N
		if (!southRow.contains(cell) && !westColumn.contains(cell)) // excluding south row and west column
			surCells.add(cell - this.N - 1); // SW
		if (!southRow.contains(cell) && !eastColumn.contains(cell)) // excluding south row and east column
			surCells.add(cell - this.N + 1); // SE
		if (!northRow.contains(cell) && !westColumn.contains(cell)) // excluding north row and west column
			surCells.add(cell + this.N - 1); // NW
		if (!northRow.contains(cell) && !eastColumn.contains(cell)) // excluding north row and east column
			surCells.add(cell + this.N + 1); // NE

		if (this.is3d) // 3d
		{
			final int zfloor = this.N * this.N;

			if (!topLevel.contains(cell)) // excluding top level
				surCells.add(cell + zfloor); // above
			if (!topLevel.contains(cell) && !westColumn.contains(cell)) // excluding top level & west wall
				surCells.add(cell + zfloor - 1); // above-west
			if (!topLevel.contains(cell) && !eastColumn.contains(cell)) // excluding top level & east wall
				surCells.add(cell + zfloor + 1);// above-east
			if (!topLevel.contains(cell) && !northRow.contains(cell)) // excluding top level & north wall
				surCells.add(cell + zfloor + this.N); // above-N
			if (!topLevel.contains(cell) && !southRow.contains(cell)) // excluding top level & south wall
				surCells.add(cell + zfloor - this.N); // above-S
			if (!topLevel.contains(cell) && !southRow.contains(cell) && !westColumn.contains(cell)) // excluding top level & south wall & west wall
				surCells.add(cell + zfloor - this.N - 1); // above-SW
			if (!topLevel.contains(cell) && !southRow.contains(cell) && !eastColumn.contains(cell)) // excluding top level & south wall & east wall
				surCells.add(cell + zfloor - this.N + 1); // above-SE
			if (!topLevel.contains(cell) && !northRow.contains(cell) && !westColumn.contains(cell)) // excluding top level & north wall & west wall
				surCells.add(cell + zfloor + this.N - 1); // above-NW
			if (!topLevel.contains(cell) && !northRow.contains(cell) && !eastColumn.contains(cell)) // excluding top level & north wall & east wall
				surCells.add(cell + zfloor + this.N + 1); // above-NE
			if (!bottomLevel.contains(cell)) // excluding bottom level
				surCells.add(cell - zfloor); // below
			if (!bottomLevel.contains(cell) && !westColumn.contains(cell)) // excluding bottom level & west wall
				surCells.add(cell - zfloor - 1); // below-W
			if (!bottomLevel.contains(cell) && !eastColumn.contains(cell)) // excluding bottom level & east wall
				surCells.add(cell - zfloor + 1); // below-E
			if (!bottomLevel.contains(cell) && !northRow.contains(cell)) // excluding bottom level & north wall
				surCells.add(cell - zfloor + this.N); // below-N
			if (!bottomLevel.contains(cell) && !southRow.contains(cell)) // excluding bottom level & south wall
				surCells.add(cell - zfloor - this.N); // below-S
			if (!bottomLevel.contains(cell) && !southRow.contains(cell) && !westColumn.contains(cell)) // excluding bottom level & south wall & west wall
				surCells.add(cell - zfloor - this.N - 1); // below-SW
			if (!bottomLevel.contains(cell) && !southRow.contains(cell) && !eastColumn.contains(cell)) // excluding bottom level & south wall & east wall
				surCells.add(cell - zfloor - this.N + 1); // below-SE
			if (!bottomLevel.contains(cell) && !northRow.contains(cell) && !westColumn.contains(cell)) // excluding bottom level & north wall & west wall
				surCells.add(cell - zfloor + this.N - 1); // below-NW
			if (!bottomLevel.contains(cell) && !northRow.contains(cell) && !eastColumn.contains(cell)) // excluding bottom level & north wall & east wall
				surCells.add(cell - zfloor + this.N + 1); // below-NE
		}

		return surCells;
	}

	// get min & max x, y, z of GD cell in integer form
	private double[] cellBorders (int cell) {
		final double[] borders = new double[6];

		final int i = cellIJK(cell)[0];
		final int j = cellIJK(cell)[1];
		final int k = cellIJK(cell)[2];

		final double ds = 1.0 / this.N;

		borders[0] = i * ds; // xmin
		borders[1] = (i + 1) * ds; // xmax
		borders[2] = j * ds; // ymin
		borders[3] = (j + 1) * ds; // ymax
		borders[4] = this.is3d ? k * ds : Double.NEGATIVE_INFINITY; // zmin
		borders[5] = this.is3d ? (k + 1) * ds : Double.NEGATIVE_INFINITY; // zmax

		return borders;
	}

	// get min & max x, y, z of QT cell in string form
	private double[] cellBorders (String cell) {

		double xmin = 0; // cell's floor-south-west corner coords initialization
		double ymin = 0;
		double zmin = this.is3d ? 0 : Double.NEGATIVE_INFINITY;

		for (int i = 0; i < cell.length(); i++) // check cellname's digits
		{
			switch(cell.charAt(i)) {
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

	// get GD cell's i, j, k
	private int[] cellIJK (int cell) {
		final int iq = cell % this.N; // get i (2d/3d)
		final int jq = !this.is3d ? (cell - iq) / this.N : ((cell - iq) / this.N) % this.N; // get j (2d/3d)
		final int kq = !this.is3d ? Integer.MIN_VALUE : ((cell - iq) / this.N - jq) / this.N; // get k (3d only)

		return new int[]{iq, jq, kq};
	}

	// return true if circle/sphere (x, y, z, r) is completely inside GD cell
	private boolean circleContainedInCell (double x, double y, double z, double r, int cell) {
		// get cell's borders (xmin, xmax, ymin, ymax, zmin, zmax)
		final double xmin = cellBorders(cell)[0];
		final double xmax = cellBorders(cell)[1];
		final double ymin = cellBorders(cell)[2];
		final double ymax = cellBorders(cell)[3];
		final double zmin = cellBorders(cell)[4];
		final double zmax = cellBorders(cell)[5];

		// if r > distance of center to cell borders, circle is not completely inside cell
		if (r > Math.abs(x - xmin))
			return false;
		if (r > Math.abs(x - xmax))
			return false;
		if (r > Math.abs(y - ymin))
			return false;
		if (r > Math.abs(y - ymax))
			return false;
		if (this.is3d && r > Math.abs(z - zmin))
			return false;
		if (this.is3d && r > Math.abs(z - zmax))
			return false;

		return true;
	}

	// return true if circle/sphere (x, y, z, r) is completely inside QT cell
	private boolean circleContainedInCell (double x, double y, double z, double r, String cell) {
		// get cell's borders (xmin, xmax, ymin, ymax, zmin, zmax)
		final double xmin = cellBorders(cell)[0];
		final double xmax = cellBorders(cell)[1];
		final double ymin = cellBorders(cell)[2];
		final double ymax = cellBorders(cell)[3];
		final double zmin = cellBorders(cell)[4];
		final double zmax = cellBorders(cell)[5];

		// if r > distance of center to cell borders, circle is not completely inside cell
		if (r > Math.abs(x - xmin))
			return false;
		if (r > Math.abs(x - xmax))
			return false;
		if (r > Math.abs(y - ymin))
			return false;
		if (r > Math.abs(y - ymax))
			return false;
		if (this.is3d && r > Math.abs(z - zmin))
			return false;
		if (this.is3d && r > Math.abs(z - zmax))
			return false;

		return true;
	}
}
