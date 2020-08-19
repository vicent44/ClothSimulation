using System.Collections;
using System.Collections.Generic;
using UnityEngine;

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

    public void ApplyForce()
    {
            Vector3 direction1 = a.Position - b.Position;
            var dirnor = direction1.magnitude;
            Vector3 directionPrev = a.Prev - b.Prev;
            Vector3 len = direction1 - directionPrev;
            var f = (dirnor - restLength) / dirnor;
            if(a.isActive) a.Position -= 0.5f * direction1 * f;
            if(b.isActive) b.Position += 0.5f * direction1 * f;

        Vector3 direction = a.Position - b.Position;

        //Debug.Log(direction);
        float dist = direction.magnitude;
        direction = direction.normalized;
        //Debug.Log(direction);

        float springForce = -elast * (dist - restLength);

        Vector3 deltaVelocity = a.Velocity - b.Velocity;
        float dampingForce = -dampi * Vector3.Dot(deltaVelocity, direction);

        Vector3 force = (springForce + dampingForce) * direction;

        /*var v1 = Vector3.Dot(direction, a.Velocity);
        var v2 = Vector3.Dot(direction, b.Velocity);

        var springDamperForce = (-elast * (restLength - dist)) - (dampi * (v1 - v2));
        */
        //var forceInternala = springDamperForce * direction;
        //var forceInternalb = -forceInternala;

        //Debug.Log(forceInternala);

        /*Vector3 direction = a.Position - b.Position;
        dist = direct
        direction = direction.normalized;

        float springForce = - elast * ()*/
        //Debug.Log(force);

        if(a.isActive)
        {
            a.AddForce(force);
            //a.AddPosition(0.5f * (a.Position - b.Position) * ((dist - restLength)/dist));
        }
        if(b.isActive)
        {
            b.AddForce(-force);
            //b.AddPosition(-0.5f * (a.Position - b.Position) * ((dist - restLength)/dist));
        }
    }
}
