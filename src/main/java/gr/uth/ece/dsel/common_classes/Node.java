package gr.uth.ece.dsel.common_classes;

import java.io.Serializable;
import java.util.HashSet;

public final class Node implements Serializable
{
	// private int low, high; // lower, higher index of sample array
	private Node nw, ne, sw, se; // 2d children
	private Node cnw, cne, csw, cse, fnw, fne, fsw, fse; // 3d children (c: ceiling, f: floor)
	private double xmin, xmax, ymin, ymax, zmin, zmax; // node boundaries
	private HashSet<Integer> contPoints = new HashSet<>(); // points contained
	
	// 2d constructor
	public Node (double xmin, double ymin, double xmax, double ymax)
	{
		this.xmin = xmin;
		this.xmax = xmax;
		this.ymin = ymin;
		this.ymax = ymax;
	}
	
	// 3d constructor
	public Node (double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
	{
		// The node will be defined using the coordinates of two opposite corners
		// Floor South West (xmin, ymin, zmin) and
		// Ceiling North East (xmax, ymax, zmax)
		this.xmin = xmin;
		this.xmax = xmax;
		this.ymin = ymin;
		this.ymax = ymax;
		this.zmin = zmin;
		this.zmax = zmax;
	}
	/*
	public int getLow()
	{
		return this.low;
	}

	public void setLow(int low)
	{
		this.low = low;
	}

	public int getHigh()
	{
		return this.high;
	}

	public void setHigh(int high)
	{
		this.high = high;
	}
	*/
	// 2d set-get
	public double getXmin()
	{
		return this.xmin;
	}

	public void setXmin(double xmin)
	{
		this.xmin = xmin;
	}

	public double getXmax()
	{
		return this.xmax;
	}

	public void setXmax(double xmax)
	{
		this.xmax = xmax;
	}

	public double getYmin()
	{
		return this.ymin;
	}

	public void setYmin(double ymin)
	{
		this.ymin = ymin;
	}

	public double getYmax()
	{
		return this.ymax;
	}

	public void setYmax(double ymax)
	{
		this.ymax = ymax;
	}
	
	public Node getNW()
	{
		return this.nw;
	}

	public void setNW(Node nW)
	{
		this.nw = nW;
	}

	public Node getNE()
	{
		return this.ne;
	}

	public void setNE(Node nE)
	{
		this.ne = nE;
	}

	public Node getSW()
	{
		return this.sw;
	}

	public void setSW(Node sW)
	{
		this.sw = sW;
	}

	public Node getSE()
	{
		return this.se;
	}

	public void setSE(Node sE)
	{
		this.se = sE;
	}
	
	// 3d set-get
	public double getZmin()
	{
		return this.zmin;
	}

	public void setZmin(double zmin)
	{
		this.zmin = zmin;
	}

	public double getZmax()
	{
		return this.zmax;
	}

	public void setZmax(double zmax)
	{
		this.zmax = zmax;
	}
	
	public Node getFNW()
	{
		return this.fnw;
	}

	public void setFNW(Node fnw)
	{
		this.fnw = fnw;
	}

	public Node getFNE()
	{
		return this.fne;
	}

	public void setFNE(Node fne)
	{
		this.fne = fne;
	}

	public Node getFSW()
	{
		return this.fsw;
	}

	public void setFSW(Node fsw)
	{
		this.fsw = fsw;
	}

	public Node getFSE()
	{
		return this.fse;
	}

	public void setFSE(Node fse)
	{
		this.fse = fse;
	}
	
	public Node getCNW()
	{
		return this.cnw;
	}

	public void setCNW(Node cnw)
	{
		this.cnw = cnw;
	}

	public Node getCNE()
	{
		return this.cne;
	}

	public void setCNE(Node cne)
	{
		this.cne = cne;
	}

	public Node getCSW()
	{
		return this.csw;
	}

	public void setCSW(Node csw)
	{
		this.csw = csw;
	}

	public Node getCSE()
	{
		return this.cse;
	}

	public void setCSE(Node cse)
	{
		this.cse = cse;
	}
	
	// contained point methods
	public void addPoints(int i)
	{
		this.contPoints.add(i);
	}
	
	public void removePoints()
	{
		// 2d
		if (this.nw != null)
		{
			this.nw.removePoints();
			this.nw.contPoints.clear();
		}
		if (this.ne != null)
		{
			this.ne.removePoints();
			this.ne.contPoints.clear();
		}
		if (this.sw != null)
		{
			this.sw.removePoints();
			this.sw.contPoints.clear();
		}
		if (this.se != null)
		{
			this.se.removePoints();
			this.se.contPoints.clear();
		}
		
		// 3d
		if (this.cnw != null)
		{
			this.cnw.removePoints();
			this.cnw.contPoints.clear();
		}
		if (this.cne != null)
		{
			this.cne.removePoints();
			this.cne.contPoints.clear();
		}
		if (this.csw != null)
		{
			this.csw.removePoints();
			this.csw.contPoints.clear();
		}
		if (this.cse != null)
		{
			this.cse.removePoints();
			this.cse.contPoints.clear();
		}
		if (this.fnw != null)
		{
			this.fnw.removePoints();
			this.fnw.contPoints.clear();
		}
		if (this.fne != null)
		{
			this.fne.removePoints();
			this.fne.contPoints.clear();
		}
		if (this.fsw != null)
		{
			this.fsw.removePoints();
			this.fsw.contPoints.clear();
		}
		if (this.fse != null)
		{
			this.fse.removePoints();
			this.fse.contPoints.clear();
		}
	}

	public HashSet<Integer> getContPoints()
	{
		return this.contPoints;
	}

	public void setContPoints(HashSet<Integer> contPoints)
	{
		this.contPoints = new HashSet<>(contPoints);
	}
	
	public void removePoint(int i)
	{
		this.contPoints.remove(i);
	}
}
