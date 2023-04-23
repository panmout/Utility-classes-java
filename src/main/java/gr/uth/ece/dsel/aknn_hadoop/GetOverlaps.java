package gr.uth.ece.dsel.aknn_hadoop;

import gr.uth.ece.dsel.common_classes.*;

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
	private String qcell;
	private Point qpoint;
	private PriorityQueue<IdDist> neighbors;

	public GetOverlaps (HashMap<String, Integer> cell_tpoints, int k, String partitioning)
	{
		this.cell_tpoints = new HashMap<>(cell_tpoints);
		this.K = k;
		this.partitioning = partitioning;
	}

	public void initializeFields (Point qp, String qc, PriorityQueue<IdDist> phase2neighbors)
	{
		this.qpoint = qp;
		this.qcell = qc;
		this.neighbors = phase2neighbors;
	}

	public void setN (int n)
	{
		this.N = n;
	}

	public void setRoot (Node root)
	{
		this.root = root;
	}

	public HashSet<String> getOverlaps()
	{
		// call appropriate function according to partitioning and iterate query points of this cell
		if (this.partitioning.equals("gd"))
			return new GetOverlapsFunctions().getOverlapsGD(this.qcell, this.qpoint, this.K, this.N, this.cell_tpoints, this.neighbors); // find its overlaps
		else if (this.partitioning.equals("qt"))
			return new GetOverlapsFunctions().getOverlapsQT(this.qcell, this.qpoint, this.K, this.root, this.cell_tpoints, this.neighbors); // find its overlaps

		return new HashSet<String>();
	}
}