using System.Collections;
using System.Collections.Generic;
using UnityEngine;



public class Hashing
{
    private int p0 = 73856093;
    private int p1 = 19349663;
    private int p2 = 83492791;

    private int gridSize;
    private float invGridSize;
    private int tableSize;

    //Instantiate hash object
    public Hashing(int _gridSize, float _invGridSize, int _tableSize)
    {
        gridSize = _gridSize;
        invGridSize = _invGridSize;
        tableSize = _tableSize;
    }
    //Hash function
    public int Hash(Vector3 coordinate)
    {
        int x = Mathf.RoundToInt(coordinate.x * invGridSize);
        int y = Mathf.RoundToInt(coordinate.y * invGridSize);
        int z = Mathf.RoundToInt(coordinate.z * invGridSize);

        return (x*p0 ^ y*p1 ^ z*p2) % tableSize;
    }
    //Box for triangles
    public List<int> TriangleBoundingBoxHashes(Vector3 p0, Vector3 p1, Vector3 p2)
    {
        int minX = Mathf.RoundToInt(Mathf.Min(new float[3]{p0.x, p1.x, p2.x}));
        int minY = Mathf.RoundToInt(Mathf.Min(new float[3]{p0.y, p1.y, p2.y}));
        int minZ = Mathf.RoundToInt(Mathf.Min(new float[3]{p0.z, p1.z, p2.z}));

        int maxX = Mathf.RoundToInt(Mathf.Max(new float[3]{p0.x, p1.x, p2.x}));
        int maxY = Mathf.RoundToInt(Mathf.Max(new float[3]{p0.y, p1.y, p2.y}));
        int maxZ = Mathf.RoundToInt(Mathf.Max(new float[3]{p0.z, p1.z, p2.z}));

        List<int> hashes = new List<int>();
        for (int x=minX; x <= maxX; x+=gridSize)
        {
            for (int y=minY; y <= maxY; y+=gridSize)
            {
                for (int z=minZ; z <= maxZ; z+=gridSize)
                {
                    hashes.Add(Mathf.Abs(Hash(new Vector3(x, y, z))));
                }
            }
        }
        return hashes;
    }
    //Box for edges
    public int LineBoxHashes(Vector3 p0, Vector3 p1)
    {
        int minX = Mathf.RoundToInt(Mathf.Min(new float[2]{p0.x, p1.x}));
        int minY = Mathf.RoundToInt(Mathf.Min(new float[2]{p0.y, p1.y}));
        int minZ = Mathf.RoundToInt(Mathf.Min(new float[2]{p0.z, p1.z}));

        int maxX = Mathf.RoundToInt(Mathf.Max(new float[2]{p0.x, p1.x}));
        int maxY = Mathf.RoundToInt(Mathf.Max(new float[2]{p0.y, p1.y}));
        int maxZ = Mathf.RoundToInt(Mathf.Max(new float[2]{p0.z, p1.z}));

        int hashes = 0;
        for (int x=minX; x <= maxX; x+=gridSize)
        {
            for (int y=minY; y <= maxY; y+=gridSize)
            {
                for (int z=minZ; z <= maxZ; z+=gridSize)
                {
                    hashes = (Mathf.Abs(Hash(new Vector3(x, y, z))));
                }
            }
        }
        return hashes;
    }
    //Box for triangles
    public int TriangleBoxHashes(Vector3 p0, Vector3 p1, Vector3 p2)
    {
        int minX = Mathf.RoundToInt(Mathf.Min(new float[3]{p0.x, p1.x, p2.x}));
        int minY = Mathf.RoundToInt(Mathf.Min(new float[3]{p0.y, p1.y, p2.y}));
        int minZ = Mathf.RoundToInt(Mathf.Min(new float[3]{p0.z, p1.z, p2.z}));

        int maxX = Mathf.RoundToInt(Mathf.Max(new float[3]{p0.x, p1.x, p2.x}));
        int maxY = Mathf.RoundToInt(Mathf.Max(new float[3]{p0.y, p1.y, p2.y}));
        int maxZ = Mathf.RoundToInt(Mathf.Max(new float[3]{p0.z, p1.z, p2.z}));

        int hashes = 0;
        for (int x=minX; x <= maxX; x+=gridSize)
        {
            for (int y=minY; y <= maxY; y+=gridSize)
            {
                for (int z=minZ; z <= maxZ; z+=gridSize)
                {
                    hashes = (Mathf.Abs(Hash(new Vector3(x, y, z))));
                }
            }
        }
        return hashes;
    }
}
