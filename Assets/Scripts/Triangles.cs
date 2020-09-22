using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Triangles
{
    public int indexTriangle { get {return indextriangle;} set {indextriangle = value;}}

    public int indexTriA { get {return indextriA;} set {indextriA = value;}}
    public int indexTriB { get {return indextriB;} set {indextriA = value;}}
    public int indexTriC { get {return indextriC;} set {indextriA = value;}}

    public Vector3 posTriA { get {return postriA;} set {postriA = value;}}
    public Vector3 posTriB { get {return postriB;} set {postriB = value;}}
    public Vector3 posTriC { get {return postriC;} set {postriC = value;}}

    public Vector3 Edge1 { get {return edge1;} set {edge1 = value;}}
    public Vector3 Edge2 { get {return edge2;} set {edge2 = value;}}
    public Vector3 Edge3 { get {return edge3;} set {edge3 = value;}}    

    public Vector3 normTri { get {return normtri;} set {normtri = value;}}

    protected int indextriA;
    protected int indextriB;
    protected int indextriC;

    int indextriangle;

    protected Vector3 postriA;
    protected Vector3 postriB;
    protected Vector3 postriC;

    Vector3 edge1;
    Vector3 edge2;
    Vector3 edge3;

    protected Vector3 normtri;

    //Instantiate triangles
    public Triangles(int a, int b, int c, int index)
    {
        indextriA = a;
        indextriB = b;
        indextriC = c;
        indextriangle = index;
    }
    //Instantiate/change position of triangle
    public void PosTriangles(Vector3 posA, Vector3 posB, Vector3 posC)
    {
        postriA = posA;
        postriB = posB;
        postriC = posC;
        edge1 = postriB - postriA;
        edge2 = postriC - postriA;
        edge3 = postriC - postriB;
        normtri = Vector3.Cross(edge1, edge2).normalized;
    }
}
