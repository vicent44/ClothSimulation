using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Particles
{
    //public List<Springs> Connection { get { return connection; } }
    public Vector3 Force { get {return force;} set {force = value;}}
    public Vector3 Velocity { get {return velocity;} set {velocity = value;}}
    public Vector3 Position { get {return position;} set {position = value;}}
    public Vector3 Prev { get {return prev;} set {prev = value;}}
    public float Mass { get {return mass;} set {mass = value;}}
    public float Area { get {return area;} set {area = value;}}

    public bool isActive { get {return active;} set {active = value; if(active) velocity = Vector3.zero;}}

    public int I { get {return i;}}
    public int J { get {return j;}}

    public GameObject particleObject;

    protected Vector3 position;
    protected Vector3 prev;

    protected Vector3 velocity;

    protected float mass;
    protected float area;

    protected Vector3 force;

    protected bool active;

    protected GameObject obj;

    protected int i;
    protected int j;

    //protected Vector3 prevVel;

    //public Vector3 externalForces;
    //public Vector3 internalForces;

    //public Vector3 totalForces;

    //List<Springs> connection;

    public Particles(Vector3 p, float m, int ii, int jj)//, float a)
    {
        mass = m;
        //area = a;

        position = prev = p;
        velocity = Vector3.zero;
        active = true;
        //obj = obje;


        i = ii;
        j = jj;

        /*var obj = GameObject.CreatePrimitive(PrimitiveType.Sphere);
        obj.transform.localScale = Vector3.one * 0.1f;
        obj.transform.position = position;*/

        //externalForces = Vector3.zero;
        //internalForces = Vector3.zero;
        //totalForces = Vector3.zero;

        //connection = new List<Springs>();
        //Debug.Log(p);
    }

    public void AddForce(Vector3 forceAdded)
    {
        force += forceAdded;
        //Debug.Log(force);
    }

    public void ResetResultantForce()
    {
        force = Vector3.zero;
    }

    public void AddPosition(Vector3 positionMod)
    {
        position += positionMod;
    }

    public void UpdateParticle(float dt)
    {
        float deltaTimeMass = (dt * dt) / mass;
        var lastPosition = position; 
        position = position * 2 - prev + deltaTimeMass * force;
        prev = lastPosition;
        velocity = (position - prev) / dt;        
    }

    /*public void Step()
    {
        //Position Steep
        var x = position - prev;
        var nextx = position + x;

        prev = position;
        position = nextx;

        //Velocity Steep
        var v = velocity - prevVel;
        var nextv = velocity + v;

        prevVel = velocity;
        velocity = nextv;

        //internalForces = -1f * (velocity - prevVel) - 1f * (position - prev) + 1f * 0.5f;
    }*/

}
