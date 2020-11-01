using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using static System.Math;
using System;


public class Springs
{
    public float Elast { get { return elast;}}
    public float Dampi { get { return dampi;}}
    public float RestLength { get { return restLength;}}
    public int SpringType { get {return springType;}}

    public Particles particleA { get {return a;} set {a = value;}}
    public Particles particleB { get {return b;} set {b = value;}}

    protected Particles a;
    protected Particles b;

    float restLength;
    float elast;
    float dampi;
    int springType;

    //Posar un argument mes per entrar la constant elastica de la molla.
    public Springs(Particles a, Particles b, float elast, float dampi, int springtype)
    {
        this.a = a;
        this.b = b;
        this.elast = elast;
        this.dampi = dampi;
        this.springType = springtype;
        this.restLength = (a.Position - b.Position).magnitude;
    }

    //Appli the force of the spring to the particles that
    //take part of it
    public void ApplyForce()
    {
        //First I make a correction to avoid superelasticity.
        Vector3 direction1 = a.Position - b.Position;
        var dirnor = direction1.magnitude;
        var f = (dirnor - restLength) /restLength;
        /*if(a.isActive) a.Position -= 0.5f * direction1 * f;
        if(b.isActive) b.Position += 0.5f * direction1 * f;*/

        //Calculation of damping and elastic forces of the spring.
        Vector3 direction = a.Position - b.Position;
        float dist = direction.magnitude;
        direction = direction.normalized;

        float springForce = -elast * (dist - restLength);

        Vector3 deltaVelocity = a.Velocity - b.Velocity;
        float dampingForce = -dampi * Vector3.Dot(deltaVelocity, direction);

        Vector3 force = (springForce + dampingForce) * direction;

        //Add the force to the particle (is it's needed, if it's not ancored)
        if(float.IsNaN(dampingForce)) Debug.Log("nan-ForceDamping");
        if(float.IsNaN(a.Velocity.x)) Debug.Log("nan-ForceDam");
        if(a.isActive)
        {
            a.AddForce(force);
        }
        if(b.isActive)
        {
            b.AddForce(-force);
        }
    }

    public void SolveConstraints()
    {
        Vector3 direction1 = a.Position - b.Position;
        var dirnor = direction1.magnitude;
        var f = (dirnor - restLength) /dirnor;
        /*if(a.isActive && b.isActive)
        {
            a.Position -= 0.5f * direction1 * f;
            b.Position += 0.5f * direction1 * f;
        }*/
        if(a.isActive) a.Position -= 0.1f * direction1 * f;
        if(b.isActive) b.Position += 0.1f * direction1 * f;
    }
}
