using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Simulate
{
    List<Particles> particles;
    List<Springs> springs;
    List<Triangles> triangles;
    Vector3 windForce;
    Transform plane;

    public Simulate(List<Particles> particles, List<Springs> springs, List<Triangles> triangles, Vector3 windforce, Transform plane)
    {
        this.particles = particles;
        this.springs = springs;
        this.triangles = triangles;
        this.windForce = windforce;
        this.plane = plane;
    }

    public void Update(float dt, List<Triangles> triangles)
    {
        this.triangles = triangles;
        ComputeTotalForces();
        //IntegratorEuler(dt);
        IntegratorVerlet(dt);
        CheckPlaneCollitions();
        //CheckSelfCollitions();
    }

    public void ComputeTotalForces()
    {
        foreach( var p in particles)
        {
            p.ResetResultantForce();
            if(p.isActive){
                //Afegir la gravetat i el vent
                var gravity = new Vector3(0, -9.8f, 0) * p.Mass;
                //var wind = 

                p.AddForce(gravity);
                p.AddForce(windForce);
            }
        }
        /*for(int i = 0; i < particles.Count; i++)
        {
            particles[i].ResetResultantForce();
            var gravity = new Vector3(0, -0.98f, 0) * particles[i].Mass;
            if(i > 15 and i < 31)
            {
                particles[i].AddForce(-gravity)
            }
            else
            {
                particles[i].AddForce(gravity);
            }
        }*/

        foreach (var s in springs)
        {
            s.ApplyForce();
        }
    }

    void IntegratorEuler(float dt)
    {
        //Euler
        foreach(var p in particles)
        {
            float deltaTimeMass = dt / p.Mass;
            p.Velocity += p.Force * deltaTimeMass;
            p.Position += p.Velocity * dt;
        }
    }

    void IntegratorVerlet(float dt)
    {

        foreach(var p in particles)
        {
            float deltaTimeMass = (dt * dt) / p.Mass;
            var lastPosition = p.Position; 
            p.Position = p.Position * 2 - p.Prev + deltaTimeMass * p.Force;
            p.Prev = lastPosition;
            p.Velocity = (p.Position - p.Prev) / dt;
        }
    }

    void CheckPlaneCollitions()
    {
        foreach(var p in particles)
        {
            /*float dx = p.Position.x - plane.position.x;
            float dy = p.Position.y - plane.position.y;
            float dz = p.Position.z - plane.position.z;*/
            Vector3 d = p.Position - plane.position;
            Vector3 normalPlane = new Vector3(0f,1f,0f);

            //float dot = normalPlane.x * dx + normalPlane.y * dy + normalPlane.z * dz;
            float dot = Vector3.Dot(normalPlane, d);
            if(dot <= 0)
            {
                p.Position -= dot * normalPlane;
                Vector3 normalVelocity = Vector3.Dot(normalPlane,p.Velocity) * p.Velocity;
                Vector3 tangencialVelocity = p.Velocity - normalVelocity;
                
                //Firts 0.5 is coeficient of friction
                //Second term is - kr * normalVelocity but kr = 0 
                //because is perfect inelastic.
                p.Velocity = Vector3.zero;//(1 - 0.9f) * tangencialVelocity;
                p.ResetResultantForce();
            }
        }
    }

    void CheckSelfCollitions()
    {
        /*for(int i = 0; i < ; i++)
        {
            for()
            {

            }
        }*/
    }

    /*void ComputeForces(Particles a, Particles b, float elast, float dampi, float length)
    {
        
        Vector3 direction = a.Position - b.Position;
        direction = direction.normalized;

        var v1 = Vector3.Dot(direction, a.Velocity);
        var v2 = Vector3.Dot(direction, b.Velocity);

        var springDamperForce = (-elast * (length - direction.magnitude)) - (dampi * (v1 - v2));
        
        var forceInternal = springDamperForce * direction;
        /*Vector3 direction = a.position - b.position;
        float currentLength = direction.magnitude;
        //Debug.Log(direction);
        direction = direction.normalized;
        //Spring force + falta multiplicar per la constant de la molla
        float springForce = (-1.0f)* elast * (currentLength - length); //rest length = 0.5
        
        //spring amortiguamiento
        Vector3 deltaVelocity = a.velocity - b.velocity;
        //Falta multiplicar per la constant de amortiguament
        float dampingForce = (-1.0f) * dampi * Vector3.Dot(deltaVelocity, direction);
        Vector3 forceInternal = (springForce + dampingForce) * direction;
        */
        //Add forces to vector forces in each particle
        //Debug.Log(direction);
        /*
        a.AddForce(forceInternal);
        b.AddForce(-forceInternal);
    }*/


        public void DrawGizmos()
        {
            for(int i = 0, n = springs.Count; i < n; i++)
            {
                var p = springs[i];
                //Gizmos.color = Color.yellow;
                //Gizmos.DrawSphere(p.position, 0.2f);

                Gizmos.color = Color.white;
                Gizmos.DrawLine(p.particleA.Position, p.particleB.Position);

                /*p.Connection.ForEach(e => {
                    var other = e.Other(p);
                    Gizmos.DrawLine(p.Position, other.Position);
                });*/
            }
        }
}









/*
//Aquest es per fer verlet perque la integracio no requereix temps


/*using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Simulate
{
    List<Particles> particles;

    public Simulate(List<Particles> particles)
    {
        this.particles = particles;
    }

    public void Simu(int iterations, float dt)
    {
        //Step();
        Solve(iterations, dt);
    }

    /*void Step()
    {
        particles.ForEach(p => {
            p.Step();
        });
    }

    void Solve(int iterations, float dt)
    {
        for(int i = 0; i < iterations; i++)
        {
            particles.ForEach(p => Solve(p));
        }
    }

    void Solve(Particles particle)
    {
        ComputeForcesExternal(particle);
        particle.Connection.ForEach(e =>
        {
            var other = e.Other(particle);
            var elast = e.Elast;
            var damp = e.Dampi;
            ComputeForces(particle, other, elast, damp, e.Length);

            //He tret de solve other i la logitud de la
            //molla perque ja esta tot calculat ara simplement
            //cauclo la nova posicio de cada particula.
            Solve(particle);
        });
    }

    void Solve(Particles a)
    {
        //Euler
        float

        var delta = a.position - b.position;
        var current = delta.magnitude;
        var f = (current - rest) / current;
        a.position -= f * 0.5f * delta;
        b.position += f * 0.5f * delta;

    }

    void ComputeForces(Particles a, Particles b, float elast, float dampi, float length)
    {
        Vector3 direction = a.position - b.position;
        float currentLength = direction.magnitude;
        direction = direction.normalized;
        //Spring force + falta multiplicar per la constant de la molla
        float springForce = (-1.0f)* elast * (currentLength - length); //rest length = 0.5
        
        //spring amortiguamiento
        Vector3 deltaVelocity = a.velocity - b.velocity;
        //Falta multiplicar per la constant de amortiguament
        float dampingForce = (-1.0f) * dampi * Vector3.Dot(deltaVelocity, direction);
        Vector3 forceInternal = (springForce + dampingForce) * direction;

        //Add forces to vector forces in each particle

        a.internalForces = forceInternal;
        b.internalForces = -forceInternal;
    }

    void ComputeForcesExternal(Particles a)
    {
        //Gravity


        //Wind

        a.externalForces = 
    }


        public void DrawGizmos()
        {
            for(int i = 0, n = particles.Count; i < n; i++)
            {
                var p = particles[i];
                //Gizmos.color = Color.yellow;
                //Gizmos.DrawSphere(p.position, 0.2f);

                Gizmos.color = Color.white;
                p.Connection.ForEach(e => {
                    var other = e.Other(p);
                    Gizmos.DrawLine(p.position, other.position);
                });
            }
        }

}

*/

