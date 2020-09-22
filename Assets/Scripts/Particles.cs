using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Particles
{
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

    //Instantiate particles
    public Particles(Vector3 p, float m, int ii, int jj)
    {
        mass = m;

        position = prev = p;
        velocity = Vector3.zero;
        active = true;

        i = ii;
        j = jj;
    }
    //Add force to the summ of forces
    public void AddForce(Vector3 forceAdded)
    {
        force += forceAdded;
    }
    //Reset the force (done it every steep)
    public void ResetResultantForce()
    {
        force = Vector3.zero;
    }
    //Update particle when the particle is ancored and you
    //change it's position with arrows.
    public void AddPosition(Vector3 positionMod)
    {
        prev = position;
        position += positionMod;
        //velocity = (position - prev)/0.01f;
    }
    //Update particle using verlet
    public void UpdateParticle(float dt)
    {
        float deltaTimeMass = (dt * dt) / mass;
        var lastPosition = position; 
        position = position * 2f - prev + deltaTimeMass * force;
        prev = lastPosition;
        velocity = (position - prev) / dt;
    }
}
