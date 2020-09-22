using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Edges
{
    public int indexEdgeA { get {return indexedgeA;} set {indexedgeA = value;}}
    public int indexEdgeB { get {return indexedgeB;} set {indexedgeB = value;}}

    public int indexEdge { get {return indexedge;} set {indexedge = value;}}


    public Vector3 posEdgeA { get {return posedgeA;} set {posedgeA = value;}}
    public Vector3 posEdgeB { get {return posedgeB;} set {posedgeB = value;}}

    int indexedgeA;
    int indexedgeB;
    int indexedge;
    Vector3 posedgeA;
    Vector3 posedgeB;

    //Instantiate edges.
    public Edges(int a, int b, int indexed)
    {
        indexedgeA = a;
        indexedgeB = b;
        indexedge = indexed;
    }
    //Change/Instantiate position of edges
    public void PosEdges(Vector3 posA, Vector3 posB)
    {
        posedgeA = posA;
        posedgeB = posB;
    }
}
