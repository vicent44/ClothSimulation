using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Triangles
{
    public int indexTriA { get {return indextriA;} set {indextriA = value;}}
    public int indexTriB { get {return indextriB;} set {indextriA = value;}}
    public int indexTriC { get {return indextriC;} set {indextriA = value;}}

    public Vector3 posTriA { get {return postriA;} set {postriA = value;}}
    public Vector3 posTriB { get {return postriB;} set {postriB = value;}}
    public Vector3 posTriC { get {return postriC;} set {postriC = value;}}

    public Vector3 normTri { get {return normtri;} set {normtri = value;}}

    protected int indextriA;
    protected int indextriB;
    protected int indextriC;

    protected Vector3 postriA;
    protected Vector3 postriB;
    protected Vector3 postriC;

    protected Vector3 normtri;

    public Triangles(int a, int b, int c)
    {
        indextriA = a;
        indextriB = b;
        indextriC = c;
    }

    public void PosTriangles(Vector3 posA, Vector3 posB, Vector3 posC)
    {
        postriA = posA;
        postriB = posB;
        postriC = posC;
        Vector3 side1 = postriB - postriA;
        Vector3 side2 = postriC - postriA;
        normtri = Vector3.Cross(side1, side2).normalized;
    }
}
