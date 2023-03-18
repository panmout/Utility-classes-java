package gr.uth.ece.dsel.common_classes;

// class Point (int id, double x, double y, char type) with constructor and set-get methods

import java.io.Serializable;

public final class Point implements Serializable
{
    private final int id;
    private final double x;
    private final double y;
    private double z = Double.NEGATIVE_INFINITY;

    public Point(int id, double x, double y)
    {
        this.id = id;
        this.x = x;
        this.y = y;
    }

    public Point(int id, double x, double y, double z) // create 3d points
    {
        this.id = id;
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public int getId()
    {
        return this.id;
    }

    public double getX()
    {
        return this.x;
    }

    public double getY()
    {
        return this.y;
    }

    public double getZ()
    {
        return this.z;
    }

    @Override
    public String toString()
    {
        if (z == Double.NEGATIVE_INFINITY) // 2d
            return String.format("%d\t%9.8f\t%9.8f", this.getId(), this.getX(), this.getY());
        else // 3d
            return String.format("%d\t%9.8f\t%9.8f\t%9.8f", this.getId(), this.getX(), this.getY(), this.getZ());
    }

    public String stringCoords()
    {
        if (z == Double.NEGATIVE_INFINITY) // 2d
            return String.format("%9.8f\t%9.8f", this.getX(), this.getY());
        else // 3d
            return String.format("%9.8f\t%9.8f\t%9.8f", this.getX(), this.getY(), this.getZ());
    }
}