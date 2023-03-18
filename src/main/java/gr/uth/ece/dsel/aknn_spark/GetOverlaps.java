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
    	// check if 3d
    	double zq = Double.NEGATIVE_INFINITY;
    	if (qpoint.getZ() != Double.NEGATIVE_INFINITY)
    	{
    		zq = qpoint.getZ();
    		this.is3d = true;
    	}
    	
    	// find query point cell
    	final double ds = 1.0 / this.N; // interval ds (cell width)
    	final int intQCell = Integer.parseInt(qcell); // get int value of query point cell
    	final int iq = intQCell % this.N; // get i (2d/3d)
    	final int jq = !this.is3d ? (intQCell - iq) / this.N : ((intQCell - iq) / this.N) % this.N; // get j (2d/3d)
    	final int kq = !this.is3d ? Integer.MIN_VALUE : ((intQCell - iq) / this.N - jq) / this.N; // get k (3d only)
    	//final int iq = (int) (xq / ds); // get i
    	//final int jq = (int) (yq / ds); // get j
    	//final int intQCell = jq * this.N + iq; // calculate int cell_id
    	
    	// if neighbors list not empty, set circle radius R to the distance of farthest neighbor
    	// else half the cell width
		double R = !this.neighbors.isEmpty() ? this.neighbors.peek().getDist() : 0.5 * ds;
		
		final double x = xq % ds; // relative x from cell (floor in 3d) SW corner
		final double y = yq % ds; // relative y from cell (floor in 3d) SW corner
		final double z = this.is3d ? zq % ds : Double.NEGATIVE_INFINITY; // relative z from cell floor SW corner (3d only)
		
		// add number of training points in this cell, 0 if null
		final long num = this.cell_tpoints.getOrDefault(qcell, 0L);
		
		// top-bottom rows, far left-right columns (rows and columns become walls in 3d)
		final HashSet<Integer> south_row = new HashSet<>(); // no S, SE, SW for cells in this set
		final HashSet<Integer> north_row = new HashSet<>(); // no N, NE, NW for cells in this set
		final HashSet<Integer> west_column = new HashSet<>(); // no W, NW, SW for cells in this set
		final HashSet<Integer> east_column = new HashSet<>(); // no E, NE, SE for cells in this set
		
		// top-bottom z-level cells
		final HashSet<Integer> bottom_level = new HashSet<>(); // no lower level (-n*n) for cells in this set
		final HashSet<Integer> top_level = new HashSet<>(); // no upper level (+n*n) for cells in this set
		
		final int zfloor = this.N * this.N; // N*N: when added/subtracted goes one floor up/down
		
		for (int i = 0; i < this.N; i++) // filling sets
		{
			if (!this.is3d) // 2d case
			{
				south_row.add(i);
				north_row.add((this.N - 1) * this.N + i);
				west_column.add(i * this.N);
				east_column.add(i * this.N + this.N - 1);
			}
			else // 3d case
			{
				for (int j = 0; j < this.N; j++)
				{
					// if 3d add j*N*N
					south_row.add(i + j * zfloor);
					north_row.add((this.N - 1) * this.N + i + j * zfloor);
					west_column.add(i * this.N + j * zfloor);
					east_column.add(i * this.N + this.N - 1 + j * zfloor);
					
					bottom_level.add(j * this.N + i);
					top_level.add(j * this.N + i + (this.N - 1) * zfloor);
				}
			}
		}
		
		// case 1: there are at least knn in this cell
		if (num >= this.K && R <= ds)
		{
			final HashSet<Integer> int_overlaps = new HashSet<>(); // set of overlapping cells

			// 2d
			// draw circle and check for overlaps
			if (x + R > ds && !east_column.contains(intQCell)) // E (excluding east column)
				int_overlaps.add(intQCell + 1);
			
			if (x < R && !west_column.contains(intQCell)) // W (excluding west column)
				int_overlaps.add(intQCell - 1);
			
			if (y + R > ds && !north_row.contains(intQCell)) // N (excluding north row)
				int_overlaps.add(intQCell + this.N);
			
			if (y < R && !south_row.contains(intQCell)) // S (excluding south row)
				int_overlaps.add(intQCell - this.N);
			
			if (x + R > ds && y + R > ds && !north_row.contains(intQCell) && !east_column.contains(intQCell)) // NE (excluding north row and east column)
			{
				final double xne = (iq + 1) * ds; // NE corner coords
				final double yne = (jq + 1) * ds;
				
				if (UtilityFunctions.square_distance(xq, yq, xne, yne) < R * R) // if NE corner is inside circle, NE cell is overlapped
					int_overlaps.add(intQCell + this.N + 1);
			}
			if (x + R > ds && y < R && !south_row.contains(intQCell) && !east_column.contains(intQCell)) // SE (excluding south row and east column)
			{
				final double xse = (iq + 1) * ds; // SE corner coords
				final double yse = jq * ds;
				
				if (UtilityFunctions.square_distance(xq, yq, xse, yse) < R * R) // if SE corner is inside circle, SE cell is overlapped
					int_overlaps.add(intQCell - this.N + 1);
			}
			if (x < R && y + R > ds && !north_row.contains(intQCell) && !west_column.contains(intQCell)) // NW (excluding north row and west column)
			{
				final double xnw = iq * ds; // NW corner coords
				final double ynw = (jq + 1) * ds;
				
				if (UtilityFunctions.square_distance(xq, yq, xnw, ynw) < R * R) // if NW corner is inside circle, NW cell is overlapped
					int_overlaps.add(intQCell + this.N - 1);
			}
			if (x < R && y < R && !south_row.contains(intQCell) && !west_column.contains(intQCell)) // SW (excluding south row and west column)
			{
				final double xsw = iq * ds; // SW corner coords
				final double ysw = jq * ds;
				
				if (UtilityFunctions.square_distance(xq, yq, xsw, ysw) < R * R) // if SW corner is inside circle, SW cell is overlapped
					int_overlaps.add(intQCell - this.N - 1);
			}
			
			// 3d
			if (this.is3d)
			{
				// z-axis cells above
				if (z + R > ds && !top_level.contains(intQCell)) // above cell, excluding top level
				{
					int_overlaps.add(intQCell + zfloor);
					
					// x-z plane
					if (x + R > ds && !east_column.contains(intQCell)) // E (excluding east wall)
					{
						final double xce = (iq + 1) * ds; // ceiling-east cell acme coords
						final double zce = (kq + 1) * ds;
						
						if (UtilityFunctions.square_distance(xq, zq, xce, zce) < R * R) // if ceiling-east cell acme is inside circle, above-east cell is overlapped
							int_overlaps.add(intQCell + 1 + zfloor);
					}
					if (x < R && !west_column.contains(intQCell)) // W (excluding west wall)
					{
						final double xcw = iq * ds;       // ceiling-west cell acme coords
						final double zcw = (kq + 1) * ds;
						
						if (UtilityFunctions.square_distance(xq, zq, xcw, zcw) < R * R) // if ceiling-west cell acme is inside circle, above-west cell is overlapped
							int_overlaps.add(intQCell - 1 + zfloor);
					}
					// y-z plane
					if (y + R > ds && !north_row.contains(intQCell)) // N (excluding north wall)
					{
						final double ycn = (jq + 1) * ds; // ceiling-north cell acme coords
						final double zcn = (kq + 1) * ds;
						
						if (UtilityFunctions.square_distance(yq, zq, ycn, zcn) < R * R) // if ceiling-north cell acme is inside circle, above-north cell is overlapped
							int_overlaps.add(intQCell + this.N + zfloor);
					}
					if (y < R && !south_row.contains(intQCell)) // S (excluding south wall)
					{
						final double ycs = jq * ds;       // ceiling-south cell acme coords
						final double zcs = (kq + 1) * ds;
						
						if (UtilityFunctions.square_distance(yq, zq, ycs, zcs) < R * R) // if ceiling-south cell acme is inside circle, above-south cell is overlapped
							int_overlaps.add(intQCell - this.N + zfloor);
					}
					// x-y-z space
					if (x + R > ds && y + R > ds && !north_row.contains(intQCell) && !east_column.contains(intQCell)) // NE (excluding north wall and east wall)
					{
						final double xcne = (iq + 1) * ds; // ceiling-NE corner coords
						final double ycne = (jq + 1) * ds;
						final double zcne = (kq + 1) * ds;
						
						if (UtilityFunctions.square_distance(xq, yq, zq, xcne, ycne, zcne) < R * R) // if ceiling-NE corner is inside sphere, above-NE cell is overlapped
							int_overlaps.add(intQCell + this.N + 1 + zfloor);
					}
					if (x + R > ds && y < R && !south_row.contains(intQCell) && !east_column.contains(intQCell)) // SE (excluding south wall and east wall)
					{
						final double xcse = (iq + 1) * ds; // ceiling-SE corner coords
						final double ycse = jq * ds;
						final double zcse = (kq + 1) * ds;
						
						if (UtilityFunctions.square_distance(xq, yq, zq, xcse, ycse, zcse) < R * R) // if ceiling-SE corner is inside circle, above-SE cell is overlapped
							int_overlaps.add(intQCell - this.N + 1 + zfloor);
					}
					if (x < R && y + R > ds && !north_row.contains(intQCell) && !west_column.contains(intQCell)) // NW (excluding north wall and west wall)
					{
						final double xcnw = iq * ds;       // ceiling-NW corner coords
						final double ycnw = (jq + 1) * ds;
						final double zcnw = (kq + 1) * ds;
						
						if (UtilityFunctions.square_distance(xq, yq, zq, xcnw, ycnw, zcnw) < R * R) // if ceiling-NW corner is inside circle, above-NW cell is overlapped
							int_overlaps.add(intQCell + this.N - 1 + zfloor);
					}
					if (x < R && y < R && !south_row.contains(intQCell) && !west_column.contains(intQCell)) // SW (excluding south wall and west wall)
					{
						final double xcsw = iq * ds; // ceiling-SW corner coords
						final double ycsw = jq * ds;
						final double zcsw = (kq + 1) * ds;
						
						if (UtilityFunctions.square_distance(xq, yq, zq, xcsw, ycsw, zcsw) < R * R) // if ceiling-SW corner is inside circle, above-SW cell is overlapped
							int_overlaps.add(intQCell - this.N - 1 + zfloor);
					}
				}
				// z-axis cells below
				if (z < R && !bottom_level.contains(intQCell)) // below cell, excluding bottom level
				{
					int_overlaps.add(intQCell - zfloor);
					
					// x-z plane
					if (x + R > ds && !east_column.contains(intQCell)) // E (excluding east wall)
					{
						final double xfe = (iq + 1) * ds; // floor-east cell acme coords
						final double zfe = kq * ds;
						
						if (UtilityFunctions.square_distance(xq, zq, xfe, zfe) < R * R) // if floor-east cell acme is inside circle, below-east cell is overlapped
							int_overlaps.add(intQCell + 1 - zfloor);
					}
					if (x < R && !west_column.contains(intQCell)) // W (excluding west wall)
					{
						final double xfw = iq * ds;       // floor-west cell acme coords
						final double zfw = kq * ds;
						
						if (UtilityFunctions.square_distance(xq, zq, xfw, zfw) < R * R) // if floor-west cell acme is inside circle, below-west cell is overlapped
							int_overlaps.add(intQCell - 1 - zfloor);
					}
					// y-z plane
					if (y + R > ds && !north_row.contains(intQCell)) // N (excluding north wall)
					{
						final double yfn = (jq + 1) * ds; // floor-north cell acme coords
						final double zfn = kq * ds;
						
						if (UtilityFunctions.square_distance(yq, zq, yfn, zfn) < R * R) // if floor-north cell acme is inside circle, below-north cell is overlapped
							int_overlaps.add(intQCell + this.N - zfloor);
					}
					if (y < R && !south_row.contains(intQCell)) // S (excluding south wall)
					{
						final double yfs = jq * ds;       // floor-south cell acme coords
						final double zfs = kq * ds;
						
						if (UtilityFunctions.square_distance(yq, zq, yfs, zfs) < R * R) // if floor-south cell acme is inside circle, below-south cell is overlapped
							int_overlaps.add(intQCell - this.N - zfloor);
					}
					if (x + R > ds && y + R > ds && !north_row.contains(intQCell) && !east_column.contains(intQCell)) // NE (excluding north wall and east wall)
					{
						final double xfne = (iq + 1) * ds; // floor-NE corner coords
						final double yfne = (jq + 1) * ds;
						final double zfne = kq * ds;
						
						if (UtilityFunctions.square_distance(xq, yq, zq, xfne, yfne, zfne) < R * R) // if floor-NE corner is inside sphere, below-NE cell is overlapped
							int_overlaps.add(intQCell + this.N + 1 - zfloor);
					}
					if (x + R > ds && y < R && !south_row.contains(intQCell) && !east_column.contains(intQCell)) // SE (excluding south wall and east wall)
					{
						final double xfse = (iq + 1) * ds; // floor-SE corner coords
						final double yfse = jq * ds;
						final double zfse = kq * ds;
						
						if (UtilityFunctions.square_distance(xq, yq, zq, xfse, yfse, zfse) < R * R) // if floor-SE corner is inside circle, below-SE cell is overlapped
							int_overlaps.add(intQCell - this.N + 1 - zfloor);
					}
					if (x < R && y + R > ds && !north_row.contains(intQCell) && !west_column.contains(intQCell)) // NW (excluding north wall and west wall)
					{
						final double xfnw = iq * ds;       // floor-NW corner coords
						final double yfnw = (jq + 1) * ds;
						final double zfnw = kq * ds;
						
						if (UtilityFunctions.square_distance(xq, yq, zq, xfnw, yfnw, zfnw) < R * R) // if floor-NW corner is inside circle, below-NW cell is overlapped
							int_overlaps.add(intQCell + this.N - 1 - zfloor);
					}
					if (x < R && y < R && !south_row.contains(intQCell) && !west_column.contains(intQCell)) // SW (excluding south wall and west wall)
					{
						final double xfsw = iq * ds; // floor-SW corner coords
						final double yfsw = jq * ds;
						final double zfsw = kq * ds;
						
						if (UtilityFunctions.square_distance(xq, yq, zq, xfsw, yfsw, zfsw) < R * R) // if floor-SW corner is inside circle, below-SW cell is overlapped
							int_overlaps.add(intQCell - this.N - 1 - zfloor);
					}
				}
			}
			
			// remove overlaps not containing training points
			int_overlaps.removeIf(cell -> !this.cell_tpoints.containsKey(String.valueOf(cell)));
			
			// subcase 1: no overlaps containing any tpoints, point goes straight to next phase
			if (int_overlaps.isEmpty())
			{
				this.overlaps.add(qcell); // add qpoint cell
				this.listComplete = true; // set status to true
			}
			
			// subcase 2: there are overlaps, additional checks must be made in next phase
			else // update overlaps list
				for (int over_cell : int_overlaps)
					this.overlaps.add(String.valueOf(over_cell));
		} // end if case 1
		
		// case 2: there are less than knn in this cell
		else
		{
			// set of surrounding cells
			final HashSet<Integer> surrounding_cells = new HashSet<>();

			// dummy set of cells to be added (throws ConcurrentModificationException if trying to modify set while traversing it)
			final HashSet<Integer> addSquaresList = new HashSet<>();
			
			int overlaps_points = 0; // total number of training points in overlaps
			
			// first element is query cell
			surrounding_cells.add(intQCell);
			
			int loopvar = 0; // loop control variable (runs until it finds >=k tpoints, then once more)
			
			while (loopvar < 2) // trying to find overlaps to fill neighbors list
			{
				overlaps_points = 0; // reset value
				addSquaresList.clear(); // clear list
				
				// getting all surrounding cells of qCell
				
				boolean runAgain = true;
				
				// keep filling set until it contains circle R inside it
				while (runAgain)
				{
					for (int square : surrounding_cells)
					{
						// x-y plane
						if (!west_column.contains(square)) // W (excluding west column)
							addSquaresList.add(square - 1);
						
						if (!east_column.contains(square)) // E (excluding east column)
							addSquaresList.add(square + 1);
						
						if (!north_row.contains(square)) // N (excluding north_row)
							addSquaresList.add(square + this.N);
						
						if (!south_row.contains(square)) // S (excluding south_row)
							addSquaresList.add(square - this.N);
						
						if (!south_row.contains(square) && !west_column.contains(square)) // SW (excluding south row and west column)
							addSquaresList.add(square - this.N - 1);
						
						if (!south_row.contains(square) && !east_column.contains(square)) // SE (excluding south row and east column)
							addSquaresList.add(square - this.N + 1);
						
						if (!north_row.contains(square) && !west_column.contains(square)) // NW (excluding north row and west column)
							addSquaresList.add(square + this.N - 1);
						
						if (!north_row.contains(square) && !east_column.contains(square)) // NE (excluding north row and east column)
							addSquaresList.add(square + this.N + 1);
						
						// 3d
						if (this.is3d)
						{
							// above cell
							if (!top_level.contains(square)) // excluding top level
								addSquaresList.add(square + zfloor);
							
							// above-west cell
							if (!top_level.contains(square) && !west_column.contains(square)) // excluding top level & west wall
								addSquaresList.add(square - 1 + zfloor);
							
							// above-east cell
							if (!top_level.contains(square) && !east_column.contains(square)) // excluding top level & east wall
								addSquaresList.add(square + 1 + zfloor);
							
							// above-north cell
							if (!top_level.contains(square) && !east_column.contains(square)) // excluding top level & north wall
								addSquaresList.add(square + this.N + zfloor);
							
							// above-south cell
							if (!top_level.contains(square) && !east_column.contains(square)) // excluding top level & south wall
								addSquaresList.add(square - this.N + zfloor);
							
							// above-south-west cell
							if (!top_level.contains(square) && !south_row.contains(square) && !west_column.contains(square)) // excluding top level & south wall & west wall
								addSquaresList.add(square - this.N - 1 + zfloor);
							
							// above-south-east cell
							if (!top_level.contains(square) && !south_row.contains(square) && !east_column.contains(square)) // excluding top level & south wall & east wall
								addSquaresList.add(square - this.N + 1 + zfloor);
							
							// above-north-west cell
							if (!top_level.contains(square) && !north_row.contains(square) && !west_column.contains(square)) // excluding top level & north wall & west wall
								addSquaresList.add(square + this.N - 1 + zfloor);
							
							// above-north-east cell
							if (!top_level.contains(square) && !north_row.contains(square) && !east_column.contains(square)) // excluding top level & north wall & east wall
								addSquaresList.add(square + this.N + 1 + zfloor);
							
							// below cell
							if (!bottom_level.contains(square)) // excluding bottom level
								addSquaresList.add(square - zfloor);
							
							// below-west cell
							if (!bottom_level.contains(square) && !west_column.contains(square)) // excluding bottom level & west wall
								addSquaresList.add(square - 1 - zfloor);
							
							// below-east cell
							if (!bottom_level.contains(square) && !east_column.contains(square)) // excluding bottom level & east wall
								addSquaresList.add(square + 1 - zfloor);
							
							// below-north cell
							if (!bottom_level.contains(square) && !east_column.contains(square)) // excluding bottom level & north wall
								addSquaresList.add(square + this.N - zfloor);
							
							// below-south cell
							if (!bottom_level.contains(square) && !east_column.contains(square)) // excluding bottom level & south wall
								addSquaresList.add(square - this.N - zfloor);
							
							// below-south-west cell
							if (!bottom_level.contains(square) && !south_row.contains(square) && !west_column.contains(square)) // excluding bottom level & south wall & west wall
								addSquaresList.add(square - this.N - 1 - zfloor);
							
							// below-south-east cell
							if (!bottom_level.contains(square) && !south_row.contains(square) && !east_column.contains(square)) // excluding bottom level & south wall & east wall
								addSquaresList.add(square - this.N + 1 - zfloor);
							
							// below-north-west cell
							if (!bottom_level.contains(square) && !north_row.contains(square) && !west_column.contains(square)) // excluding bottom level & north wall & west wall
								addSquaresList.add(square + this.N - 1 - zfloor);
							
							// below-north-east cell
							if (!bottom_level.contains(square) && !north_row.contains(square) && !east_column.contains(square)) // excluding bottom level & north wall & east wall
								addSquaresList.add(square + this.N + 1 - zfloor);
						}
					}
					
					surrounding_cells.addAll(addSquaresList); // add new squares to original set
					
					// 2d
					if (!this.is3d)
					{
						// boolean variables to check if surrounding cells include the circle with radius R
						boolean stopRunX = false;
						boolean stopRunY = false;
						
						int maxI = iq; // min & max column index of surrounding cells at query cell row
						int minI = iq;
						
						for (int i = 0; i < this.N; i++) // running through columns 0 to N
						{
							if (surrounding_cells.contains(jq * this.N + i)) // getting cells at query cell row (jq)
							{
								maxI = Math.max(i, maxI);
								minI = Math.min(i, minI);
							}
						}
						
						if ((maxI - minI) * ds > 2 * R) // if surrounding cells width is more than 2*R, set stop var to 'true'
							stopRunX = true;
						
						int maxJ = jq; // min & max row index of surrounding cells at query cell column
						int minJ = jq;
						
						for (int j = 0; j < this.N; j++) // running through columns 0 to N
						{
							if (surrounding_cells.contains(j * this.N + iq)) // getting cells at query cell column (iq)
							{
								maxJ = Math.max(j, maxJ);	
								minJ = Math.min(j, minJ);
							}
						}
						
						if ((maxJ - minJ) * ds > 2 * R) // if surrounding cells width is more than 2*R, set stop var to 'true'
							stopRunY = true;
						
						// if all stop vars are set to 'true', stop loop
						if (stopRunX && stopRunY)
							runAgain = false;
					}
					// 3d case
					else
					{
						// boolean variables to check if surrounding cells include the sphere with radius R
						boolean stopRunX = false;
						boolean stopRunY = false;
						boolean stopRunZ = false;
						
						int maxI = iq; // min & max row index of surrounding cells at query cell row
						int minI = iq;
						
						for (int i = 0; i < this.N; i++) // running through rows 0 to N
						{
							if (surrounding_cells.contains(i + jq * this.N + kq * zfloor)) // getting cells at query cell column (jq) and z-level (kq)
							{
								maxI = Math.max(i, maxI);	
								minI = Math.min(i, minI);
							}
						}
						
						if ((maxI - minI)*ds > 2*R) // if surrounding cells width is more than 2*R, set stop var to 'true'
							stopRunX = true;
						
						int maxJ = jq; // min & max column index of surrounding cells at celli's column
						int minJ = jq;
						
						for (int j = 0; j < this.N; j++) // running through columns 0 to N
						{
							if (surrounding_cells.contains(iq + j * this.N + kq * zfloor)) // getting cells at celli's row (i0) and z-level (k0)
							{
								maxJ = Math.max(j, maxJ);	
								minJ = Math.min(j, minJ);
							}
						}
						
						if ((maxJ - minJ)*ds > 2*R) // if surrounding cells width is more than 2*R, set stop var to 'true'
							stopRunY = true;
						
						int maxK = kq; // min & max z index of surrounding cells at celli's row (iq) and column (jq)
						int minK = kq;
						
						for (int k = 0; k < this.N; k++) // running through z-axis cells 0 to N
						{
							if (surrounding_cells.contains(iq + jq * this.N + k * zfloor)) // getting cells at celli's row (iq) and column (jq)
							{
								maxK = Math.max(k, maxK);	
								minK = Math.min(k, minK);
							}
						}
						
						if ((maxK - minK) * ds > 2 * R) // if surrounding cells width is more than 2*R, set stop var to 'true'
							stopRunZ = true;
						
						// if all stop vars are set to 'true', stop loop
						if (stopRunX && stopRunY && stopRunZ)
							runAgain = false;
					}
				} // end while (runagain)
				
				// checking for overlaps in surroundings
				for (int square: surrounding_cells)
				{
					// proceed only if cell contains any training points and skip query point cell
					if (square != intQCell && this.cell_tpoints.containsKey(String.valueOf(square)))
					{
						final String sq = String.valueOf(square);
						
						// 2d case
						if (!this.is3d)
						{
							// cell_id = j*n + i
							final int i = square % this.N;
							final int j = (square - i) / this.N;
							// get cell center coordinates
							final double cx = i * ds + ds / 2;
							final double cy = j * ds + ds / 2;
							// circle center to cell center distance
							final double centers_dist_x = Math.abs(xq - cx);
							final double centers_dist_y = Math.abs(yq - cy);
							
							// check circle - cell collision
							if (i > iq && j == jq) // to the east of query cell, same row
							{
								if (xq + R > i * ds) // checking collision with cell's west wall
									this.overlaps.add(sq); // there is collision, add cell to overlaps
							}
							else if (i > iq && j > jq) // to the north-east of query cell
							{
								if (centers_dist_x < R + ds / 2 && centers_dist_y < R + ds / 2) // if centers are close enough
									if (UtilityFunctions.square_distance(xq, yq, i * ds, j * ds) < R * R) // if also SW corner is inside circle
										this.overlaps.add(sq); // there is collision, add cell to overlaps
							}
							else if (i == iq && j > jq) // to the north of query cell, same column
							{
								if (yq + R > j * ds) // checking collision with cell's south wall
									this.overlaps.add(sq); // there is collision, add cell to overlaps
							}
							else if (i < iq && j > jq) // to the north-west of query cell
							{
								if (centers_dist_x < R + ds / 2 && centers_dist_y < R + ds / 2) // if centers are close enough
									if (UtilityFunctions.square_distance(xq, yq, (i + 1) * ds, j * ds) < R * R) // if also SE corner is inside circle
										this.overlaps.add(sq); // there is collision, add cell to overlaps
							}
							else if (i < iq && j == jq) // to the west of query cell, same row
							{
								if (xq - R < (i + 1) * ds) // checking collision with cell's east wall
									this.overlaps.add(sq); // there is collision, add cell to overlaps
							}
							else if (i < iq && j < jq) // to the south-west of query cell
							{
								if (centers_dist_x < R + ds / 2 && centers_dist_y < R + ds / 2) // if centers are close enough
									if (UtilityFunctions.square_distance(xq, yq, (i + 1) * ds, (j + 1) * ds) < R * R) // if also NE corner is inside circle
										this.overlaps.add(sq); // there is collision, add cell to overlaps
							}
							else if (i == iq && j < jq) // to the south of query cell, same column
							{
								if (yq - R < (j + 1) * ds) // checking collision with cell's north wall
									this.overlaps.add(sq); // there is collision, add cell to overlaps
							}
							else if (i > iq && j < jq) // to the south-east of query cell
							{
								if (centers_dist_x < R + ds / 2 && centers_dist_y < R + ds / 2) // if centers are close enough
									if (UtilityFunctions.square_distance(xq, yq, i * ds, (j + 1) * ds) < R * R) // if also NE corner is inside circle
										this.overlaps.add(sq); // there is collision, add cell to overlaps
							}
						}
						// 3d case
						else
						{
							// cell_id = i + j * N + k * N * N
							final int i = square % this.N;
							final int j = ((square - i) / this.N) % this.N;
							final int k = ((square - i) / this.N - j) / this.N; // get k
							// get cell center coordinates
							final double cx = i * ds + ds / 2;
							final double cy = j * ds + ds / 2;
							final double cz = k * ds + ds / 2;
							// circle/sphere center to cell center distance
							final double centers_dist_x = Math.abs(xq - cx);
							final double centers_dist_y = Math.abs(yq - cy);
							final double centers_dist_z = Math.abs(zq - cz);
														
							// check circle/sphere - cell collision
							// cells directly above or below
							if (i == iq && j == jq)
							{
								// lower z level
								if (k < kq)
								{
									if (zq - R < (k + 1) * ds) // checking collision with cell's ceiling (z-axis)
										this.overlaps.add(sq); // there is collision, add cell to overlaps
								}
								// upper z level
								if (k > kq)
								{
									if (zq + R > k * ds) // checking collision with cell's floor (z-axis)
										this.overlaps.add(sq); // there is collision, add cell to overlaps
								}
							}
							// cells to the east, same row
							if (i > iq && j == jq)
							{
								// same z level
								if (k == kq)
								{
									if (xq + R > i * ds) // checking collision with cell's west wall (x-y plane)
										this.overlaps.add(sq); // there is collision, add cell to overlaps
								}
								// lower z level
								else if (k < kq)
								{
									if ((centers_dist_x < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on x-z plane are close enough
									{
										if (UtilityFunctions.square_distance(xq, zq, i * ds, (k + 1) * ds) < R * R) // if also ceiling-west acme is inside circle
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
								// upper z level
								else // k > kq
								{
									if ((centers_dist_x < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on x-z plane are close enough
									{
										if (UtilityFunctions.square_distance(xq, zq, i * ds, k * ds) < R * R) // if also floor-west acme is inside circle
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
							}
							// cells to the north-east
							else if (i > iq && j > jq)
							{
								// same z level
								if (k == kq)
								{
									if ((centers_dist_x < R + ds / 2) && (centers_dist_y < R + ds / 2)) // if centers on x-y plane are close enough
									{
										if (UtilityFunctions.square_distance(xq, yq, i * ds, j * ds) < R * R) // if also SW acme is inside circle
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
								// lower z level
								else if (k < kq)
								{
									if ((centers_dist_x < R + ds / 2) && (centers_dist_y < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on x-y-z space are close enough
									{
										if (UtilityFunctions.square_distance(xq, yq, zq, i * ds, j * ds, (k + 1) * ds) < R * R) // if also ceiling-SW corner is inside sphere
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
								// upper z level
								else // k > kq
								{
									if ((centers_dist_x < R + ds / 2) && (centers_dist_y < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on x-y-z space are close enough
									{
										if (UtilityFunctions.square_distance(xq, yq, zq, i * ds, j * ds, k * ds) < R * R) // if also floor-SW corner is inside sphere
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
							}
							// cells to the north, same column
							else if (i == iq && j > jq)
							{
								// same z level
								if (k == kq)
								{
									if (yq + R > j * ds) // checking collision with cell's south wall (x-y plane)
										this.overlaps.add(sq); // there is collision, add cell to overlaps
								}
								// lower z level
								else if (k < kq)
								{
									if ((centers_dist_y < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on y-z plane are close enough
									{
										if (UtilityFunctions.square_distance(yq, zq, j * ds, (k + 1) * ds) < R * R) // if also ceiling-south acme is inside circle
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
								// upper z level
								else // k > kq
								{
									if ((centers_dist_y < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on y-z plane are close enough
									{
										if (UtilityFunctions.square_distance(yq, zq, j * ds, k * ds) < R * R) // if also floor-south acme is inside circle
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
							}
							// cells to the north-west
							else if (i < iq && j > jq)
							{
								// same z level
								if (k == kq)
								{
									if ((centers_dist_x < R + ds / 2) && centers_dist_y < R + ds / 2) // if centers on x-y plane are close enough
									{
										if (UtilityFunctions.square_distance(xq, yq, (i + 1) * ds, j * ds) < R * R) // if also SE acme is inside circle
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
								// lower z level
								else if (k < kq)
								{
									if ((centers_dist_x < R + ds / 2) && (centers_dist_y < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on x-y-z space are close enough
									{
										if (UtilityFunctions.square_distance(xq, yq, zq, (i + 1) * ds, j * ds, (k + 1) * ds) < R * R) // if also ceiling-SE corner is inside sphere
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
								// upper z level
								else // k > kq
								{
									if ((centers_dist_x < R + ds / 2) && (centers_dist_y < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on x-y-z space are close enough
									{
										if (UtilityFunctions.square_distance(xq, yq, zq, (i + 1) * ds, j * ds, k * ds) < R * R) // if also floor-SE corner is inside sphere
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
							}
							// cells to the west, same row
							else if (i < iq && j == jq)
							{
								// same z level
								if (k == kq)
								{
									if (xq - R < (i + 1) * ds) // checking collision with cell's east wall (x-y plane)
										this.overlaps.add(sq); // there is collision, add cell to overlaps
								}
								// lower z level
								else if (k < kq)
								{
									if ((centers_dist_x < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on x-z plane are close enough
									{
										if (UtilityFunctions.square_distance(xq, zq, (i + 1) * ds, (k + 1) * ds) < R * R) // if also ceiling-east acme is inside circle
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
								// upper z level
								else // k > kq
								{
									if ((centers_dist_x < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on x-z plane are close enough
									{
										if (UtilityFunctions.square_distance(xq, zq, (i + 1) * ds, k * ds) < R * R) // if also floor-east acme is inside circle
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
							}
							// cells to the south-west
							else if (i < iq && j < jq)
							{
								// same z level
								if (k == kq)
								{
									if ((centers_dist_x < R + ds / 2) && centers_dist_y < R + ds / 2) // if centers on x-y plane are close enough
									{
										if (UtilityFunctions.square_distance(xq, yq, (i + 1) * ds, (j + 1) * ds) < R * R) // if also NE acme is inside circle
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
								// lower z level
								else if (k < kq)
								{
									if ((centers_dist_x < R + ds / 2) && (centers_dist_y < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on x-y-z space are close enough
									{
										if (UtilityFunctions.square_distance(xq, yq, zq, (i + 1) * ds, (j + 1) * ds, (k + 1) * ds) < R * R) // if also ceiling-NE corner is inside sphere
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
								// upper z level
								else // k > kq
								{
									if ((centers_dist_x < R + ds / 2) && (centers_dist_y < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on x-y-z space are close enough
									{
										if (UtilityFunctions.square_distance(xq, yq, zq, (i + 1) * ds, (j + 1) * ds, k * ds) < R * R) // if also floor-NE corner is inside sphere
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
							}
							// cells to the south, same column
							else if (i == iq && j < jq)
							{
								// same z level
								if (k == kq)
								{
									if (yq - R < (j + 1) * ds) // checking collision with cell's north wall (x-y plane)
										this.overlaps.add(sq); // there is collision, add cell to overlaps
								}
								// lower z level
								else if (k < kq)
								{
									if ((centers_dist_y < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on y-z plane are close enough
									{
										if (UtilityFunctions.square_distance(yq, zq, (j + 1) * ds, (k + 1) * ds) < R * R) // if also ceiling-north acme is inside circle
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
								// upper z level
								else // k > kq
								{
									if ((centers_dist_y < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on y-z plane are close enough
									{
										if (UtilityFunctions.square_distance(yq, zq, (j + 1) * ds, k * ds) < R * R) // if also floor-north acme is inside circle
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
							}
							// cells to the south-east
							else if (i > iq && j < jq)
							{
								// same z level
								if (k == kq)
								{
									if ((centers_dist_x < R + ds / 2) && centers_dist_y < R + ds / 2) // if centers on x-y plane are close enough
									{
										if (UtilityFunctions.square_distance(xq, yq, i * ds, (j + 1) * ds) < R * R) // if also NE corner is inside circle
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
								// lower z level
								else if (k < kq)
								{
									if ((centers_dist_x < R + ds / 2) && (centers_dist_y < R + ds / 2) && (centers_dist_z < R + ds / 2)) // if centers on x-y-z space are close enough
									{
										if (UtilityFunctions.square_distance(xq, yq, zq, i * ds, (j + 1) * ds, (k + 1) * ds) < R * R) // if also ceiling-NW corner is inside sphere
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
								// upper z level
								else // k > kq
								{
									if ((centers_dist_x < R + ds/2) && (centers_dist_y < R + ds/2) && (centers_dist_z < R + ds/2)) // if centers on x-y-z space are close enough
									{
										if (UtilityFunctions.square_distance(xq, yq, zq, i * ds, (j + 1) * ds, k * ds) < R*R) // if also floor-NW corner is inside sphere
											this.overlaps.add(sq); // there is collision, add cell to overlaps
									}
								}
							}
						}
					}
				} // end for (traverse surroundings)
				
				// now find total training points from overlaps
				if (!this.overlaps.isEmpty())
					for (String s : this.overlaps)
						overlaps_points += this.cell_tpoints.get(s); // add this overlap's training points
				
				R += 0.5 * ds; // increase radius by half ds
				
				// if k neighbors found, run loop one more time and increase radius by the diagonal of a cell
				if (num + overlaps_points >= this.K)
				{
					loopvar++;
					if (!this.is3d) // 2d pythagorean
						R += Math.sqrt(2) * ds;
					else // 3d pythagorean
						R += Math.sqrt(3) * ds;
				}
			} // end while (loopvar)
		} // end case 2 else
	} // end getOverlapsGD
	
	// find quadtree query overlaps
	private void getOverlapsQT (String qcell, Point qpoint)
	{
		// read query point coordinates and neighbors list
    	final double xq = qpoint.getX();
    	final double yq = qpoint.getY();
    	// check if 3d
    	double zq = Double.NEGATIVE_INFINITY;
    	if (qpoint.getZ() != Double.NEGATIVE_INFINITY)
    	{
    		zq = qpoint.getZ();
    		this.is3d = true;
    	}
    	
		/* If
		 * root cell side length = L
		 * and for example
		 * cell id = 3012 (4 digits)
		 * then cell's length = L / (2 ^ 4)
		 */
    	
    	// total number of training points in this cell, 0 if null
    	final long num = this.cell_tpoints.getOrDefault(qcell, 0L);
    	
    	final double ds = 1.0 / Math.pow(2, qcell.length()); // ds = query cell width
		
    	// if neighbors list not empty, set circle/sphere radius R to the distance of farthest neighbor
    	// else half the cell width
		double R = !this.neighbors.isEmpty() ? this.neighbors.peek().getDist() : 0.5 * ds;
		
		// case 1: there are at least knn in this cell
		if (num >= this.K)
		{
			// draw circle/sphere and check for overlaps
			// 2d
			if (!this.is3d)
				rangeQuery(xq, yq, R, this.root, "");
			// 3d
			else
				rangeQuery(xq, yq, zq, R, this.root, "");
			
			// remove containing cell
			this.overlaps.remove(qcell);
			
			// remove overlaps not containing training points
			this.overlaps.removeIf(s -> !this.cell_tpoints.containsKey(s));
			
			// subcase 1: no overlaps containing any tpoints, point goes straight to next phase
			if (this.overlaps.isEmpty())
			{
				this.overlaps.add(qcell); // add qpoint cell
				this.listComplete = true; // set status to true
			}
			
			// subcase 2: there are overlaps, additional checks must be made in next phase
			// nothing to do here, overlaps are already updated
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
			
			int overlaps_points = 0; // total number of training points in overlaps
			
			int loopvar = 0; // loop control variable (runs until it finds >=k tpoints, then once more)
			while (loopvar < 2) // trying to find overlaps to fill k-nn
			{
				overlaps_points = 0; // reset value
				this.overlaps.clear(); // clear overlaps list
				
				// draw circle/sphere and check for overlaps
				// 2d
				if (!this.is3d)
					rangeQuery(xq, yq, r1, this.root, "");
				// 3d
				else
					rangeQuery(xq, yq, zq, r1, this.root, "");
				
				// remove containing cell
				this.overlaps.remove(qcell);
				
				for (String cell : this.overlaps)
					if (this.cell_tpoints.containsKey(cell)) // count points from non-empty cells
						overlaps_points += this.cell_tpoints.get(cell); // add this overlap's training points
				
				r1 += 0.1 * r1; // increase radius by 10%
				
				// if k neighbors found (first time only), run loop one more time and set r1 equal to the maximum distance from ipoint to all overlapped cells
				if ((num + overlaps_points >= this.K) && (loopvar == 0))
				{
					loopvar++;
					double maxSqrDist = 0; // square of maximum distance found so far (ipoint to cell)
					for (String cell : this.overlaps) // for every cell in overlaps
					{
						if (this.cell_tpoints.containsKey(cell)) // only non-empty cells
						{
							// 2d
							if (!this.is3d)
							{
								double x0 = 0; // cell's lower left corner coords initialization
								double y0 = 0;
								for (int i = 0; i < cell.length(); i++) // check cellname's digits
								{
									switch(cell.charAt(i))
									{
										case '0':
											y0 += 1.0 / Math.pow(2, i + 1); // if digit = 0 increase y0
											break;
										case '1':
											x0 += 1.0 / Math.pow(2, i + 1); // if digit = 1 increase x0
											y0 += 1.0 / Math.pow(2, i + 1); // and y0
											break;
										case '3':
											x0 += 1.0 / Math.pow(2, i + 1); // if digit = 3 increase x0
											break;
									}
								}
								double s = 1.0 / Math.pow(2, cell.length()); // cell side length
								/* cell's lower left corner: x0, y0
								 *        upper left corner: x0, y0 + s
								 *        upper right corner: x0 + s, y0 + s
								 *        lower right corner: x0 + s, y0
								 */
								 double maxX = Math.max(Math.abs(xq - x0), Math.abs(xq - x0 - s)); // maximum ipoint distance from cell in x-direction
								 double maxY = Math.max(Math.abs(yq - y0), Math.abs(yq - y0 - s)); // maximum ipoint distance from cell in y-direction
								 if (maxX * maxX + maxY * maxY > maxSqrDist)
									 maxSqrDist = maxX * maxX + maxY * maxY; // replace current maximum squared distance if new is bigger
							}
							// 3d
							else
							{
								double x0 = 0; // cell's floor south west corner coords initialization
								double y0 = 0;
								double z0 = 0;
								for (int i = 0; i < cell.length(); i++) // check cellname's digits
								{
									switch(cell.charAt(i))
									{
										case '0':
											y0 += 1.0 / Math.pow(2, i + 1); // if digit = 0 increase y0
											break;
										case '1':
											x0 += 1.0 / Math.pow(2, i + 1); // if digit = 1 increase x0
											y0 += 1.0 / Math.pow(2, i + 1); // and y0
											break;
										case '3':
											x0 += 1.0 / Math.pow(2, i + 1); // if digit = 3 increase x0
											break;
										case '4':
											y0 += 1.0 / Math.pow(2, i + 1); // if digit = 4 increase y0
											z0 += 1.0 / Math.pow(2, i + 1); // and z0
											break;
										case '5':
											x0 += 1.0 / Math.pow(2, i + 1); // if digit = 5 increase x0
											y0 += 1.0 / Math.pow(2, i + 1); // and y0
											z0 += 1.0 / Math.pow(2, i + 1); // and z0
											break;
										case '6':
											z0 += 1.0 / Math.pow(2, i + 1); // if digit = 6 increase z0
											break;
										case '7':
											x0 += 1.0 / Math.pow(2, i + 1); // if digit = 7 increase x0
											z0 += 1.0 / Math.pow(2, i + 1); // and z0
											break;
									}
								}
								double s = 1.0 / Math.pow(2, cell.length()); // cell side length
								/* cell's xmin = x0
								 *        xmax = x0 + s
								 *        ymin = y0
								 *        ymax = y0 + s
								 * 	      zmin = z0
								 *        zmax = z0 + s
								 */
								 double maxX = Math.max(Math.abs(xq - x0), Math.abs(xq - x0 - s)); // maximum ipoint distance from cell in x-direction
								 double maxY = Math.max(Math.abs(yq - y0), Math.abs(yq - y0 - s)); // maximum ipoint distance from cell in y-direction
								 double maxZ = Math.max(Math.abs(zq - z0), Math.abs(zq - z0 - s)); // maximum ipoint distance from cell in z-direction
								 if (maxX * maxX + maxY * maxY + maxZ * maxZ > maxSqrDist)
									 maxSqrDist = maxX * maxX + maxY * maxY + maxZ * maxZ; // replace current maximum squared distance if new is bigger
							}
						 } // end if
					 } // end for
					 r1 = Math.sqrt(maxSqrDist); // run overlaps check once more with this radius
				} // end if
				else if (loopvar == 1)
					loopvar++;
			} // end while
		}
	}
	
	private void rangeQuery (double x, double y, double r, Node node, String address)
	{
		if (node.getNW() == null) // leaf node
			this.overlaps.add(address);
		
		// internal node
		else
		{
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
	
	private boolean intersect (double x, double y, double r, Node node)
	{
		// if point is inside cell return true
		if (x >= node.getXmin() && x <= node.getXmax() && y >= node.getYmin() && y <= node.getYmax())
			return true;
		
		// check circle - cell collision
		final double ds = node.getXmax() - node.getXmin(); // cell's width
		
		// get cell center coordinates
		final double xc = (node.getXmin() + node.getXmax()) / 2;
		final double yc = (node.getYmin() + node.getYmax()) / 2;
		
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
		final double corner_dist_sq = Math.pow(centers_dist_x - ds / 2, 2) + Math.pow(centers_dist_y - ds / 2, 2);
		
		return corner_dist_sq <= r * r;
	}
	
	private void rangeQuery(double x, double y, double z, double r, Node node, String address)
	{
		if (node.getFNW() == null) // leaf node
			this.overlaps.add(address);
		// internal node
		else
		{
			if (intersect(node.getFNW(), x, y, z, r))
				rangeQuery(x, y, z, r, node.getFNW(), address + "0");
			
			if (intersect(node.getFNE(), x, y, z, r))
				rangeQuery(x, y, z, r, node.getFNE(), address + "1");
			
			if (intersect(node.getFSW(), x, y, z, r))
				rangeQuery(x, y, z, r, node.getFSW(), address + "2");
			
			if (intersect(node.getFSE(), x, y, z, r))
				rangeQuery(x, y, z, r, node.getFSE(), address + "3");
			
			if (intersect(node.getCNW(), x, y, z, r))
				rangeQuery(x, y, z, r, node.getCNW(), address + "4");
			
			if (intersect(node.getCNE(), x, y, z, r))
				rangeQuery(x, y, z, r, node.getCNE(), address + "5");
			
			if (intersect(node.getCSW(), x, y, z, r))
				rangeQuery(x, y, z, r, node.getCSW(), address + "6");
			
			if (intersect(node.getCSE(), x, y, z, r))
				rangeQuery(x, y, z, r, node.getCSE(), address + "7");
		}
	}
	
	private boolean intersect(Node node, double x, double y, double z, double r)
	{
		// if point is inside cell return true
		if (x >= node.getXmin() && x <= node.getXmax() && y >= node.getYmin() && y <= node.getYmax() && z >= node.getZmin() && z <= node.getZmax())
			return true;
		
		// check sphere - cell collision
		double ds = node.getXmax() - node.getXmin(); // cell's width
		// get cell center coordinates
		double xc = (node.getXmin() + node.getXmax()) / 2;
		double yc = (node.getYmin() + node.getYmax()) / 2;
		double zc = (node.getZmin() + node.getZmax()) / 2;
		// sphere center to cell center distance
		double centers_dist_x = Math.abs(x - xc);
		double centers_dist_y = Math.abs(y - yc);
		double centers_dist_z = Math.abs(z - zc);
		
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
		double corner_dist_sq = Math.pow(centers_dist_x - ds / 2, 2) + Math.pow(centers_dist_y - ds / 2, 2) + Math.pow(centers_dist_z - ds / 2, 2);
		return corner_dist_sq <= r * r;
	}
}
