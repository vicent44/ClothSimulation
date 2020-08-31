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

    public struct SPHash
    {
        public List<int> indices;
    }
    TriangleIntersection intersec;
    Hashing hash;
    public int gridSize;

    public Simulate(List<Particles> particles, List<Springs> springs, List<Triangles> triangles, Vector3 windforce, Transform plane, int gridSize)
    {
        this.particles = particles;
        this.springs = springs;
        this.triangles = triangles;
        this.windForce = windforce;
        this.plane = plane;
        this.gridSize = gridSize;
        //intersec = new TriangleIntersection();
    }

    public void Update(float dt, List<Triangles> triangles)
    {
        this.triangles = triangles;
        ComputeTotalForces();
        //IntegratorEuler(dt);
        IntegratorVerlet(dt);
        CheckPlaneCollitions(dt);
        CheckSelfCollitions();
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
            if(p.isActive) p.UpdateParticle(dt);
            /*float deltaTimeMass = (dt * dt) / p.Mass;
            var lastPosition = p.Position; 
            p.Position = p.Position * 2 - p.Prev + deltaTimeMass * p.Force;
            p.Prev = lastPosition;
            p.Velocity = (p.Position - p.Prev) / dt;*/
        }
    }

    void CheckPlaneCollitions(float dt)
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
                Vector3 normalVelocity = Vector3.Dot(normalPlane,p.Velocity) * normalPlane;
                Vector3 tangencialVelocity = p.Velocity - normalVelocity;
                Vector3 normalForce = Vector3.Dot(p.Force, normalPlane) * normalPlane;
                Vector3 tangencialForce = p.Force - normalForce;
                //Firts 0.5 is coeficient of friction
                //Second term is - kr * normalVelocity but kr = 0 
                //because is perfect inelastic.
                p.Velocity = Vector3.zero;//tangencialVelocity; // -normalVelocity;//Vector3.zero;//(1 - 0.9f) * tangencialVelocity;
                p.AddForce(tangencialForce - 0.5f * normalForce.magnitude * (tangencialForce/tangencialForce.magnitude));
                //p.AddForce(- p.Mass * tangencialVelocity/dt);
                //p.ResetResultantForce();
            }
        }
    }

    void CheckSelfCollitions()
    {
        hash = new Hashing(gridSize, 1.0f/(float)gridSize, 1723);
        SPHash[] spHash = new SPHash[1723];

        for(int v = 0; v < particles.Count; v++)
        {
            int has = Mathf.Abs(hash.Hash(particles[v].Position));

            if(spHash[has].indices == null) spHash[has].indices = new List<int>();
            spHash[has].indices.Add(particles[v].I);
        }

        for(int t = 0; t < triangles.Count; t++)
        {
            var tri = triangles[t];
            var p0 = particles[tri.indexTriA].Position;
            var p1 = particles[tri.indexTriB].Position;
            var p2 = particles[tri.indexTriC].Position;

            float w0 = 1f/particles[tri.indexTriA].Mass;
            float w1 = 1f/particles[tri.indexTriB].Mass;
            float w2 = 1f/particles[tri.indexTriC].Mass;

            List<int> hashes = hash.TriangleBoundingBoxHashes(p0, p1, p2);

            for(int h = 0; h < hashes.Count; h++)
            {
                if(spHash[h].indices != null)
                {
                    for(int sph = 0; sph < spHash[h].indices.Count; sph++)
                    {
                        int idx = spHash[h].indices[sph];
                        if(idx != particles[tri.indexTriA].I && idx != particles[tri.indexTriB].I && idx != particles[tri.indexTriC].I)
                        {
                            Vector3 p = particles[idx].Position;
                            float w = 1f/particles[idx].Mass;

                            Vector3 corr, corr0, corr1, corr2, normalTri;
                            float val0, val1, val2;
                            if(TrianglePointDistance(
                                p, w,
                                p0, w0,
                                p1, w1,
                                p2, w2,
                                0.02f, 5000f, 0.0f,
                                out corr, out corr0, out corr1, out corr2, out normalTri,
                                out val0, out val1, out val2))
                            {
                                
                                //Vector3 normalVeloci = Vector3.Dot(triangles[t].normTri, particles[idx].Velocity) * particles[idx].Velocity;
                                //Vector3 tangencialVelocity = particles[idx].Velocity - normalVeloci;
                                
                                /*particles[idx].Velocity = Vector3.zero;
                                particles[idx].ResetResultantForce();

                                particles[tri.indexTriA].Velocity = Vector3.zero;
                                particles[tri.indexTriA].ResetResultantForce();

                                particles[tri.indexTriB].Velocity = Vector3.zero;
                                particles[tri.indexTriB].ResetResultantForce();

                                particles[tri.indexTriC].Velocity = Vector3.zero;
                                particles[tri.indexTriC].ResetResultantForce();*/
                                Debug.Log(particles[idx].Position);
                                /*particles[idx].Position += corr;
                                particles[tri.indexTriA].Position += corr0;
                                particles[tri.indexTriB].Position += corr1;
                                particles[tri.indexTriC].Position += corr2;*/

                                //Position
                                particles[idx].Position += corr;
                                particles[tri.indexTriA].Position += corr0;
                                particles[tri.indexTriB].Position += corr1;
                                particles[tri.indexTriC].Position += corr2;

                                //Velocity
                                Vector3 normalVelocity0 = Vector3.Dot(normalTri, particles[idx].Velocity) * normalTri;
                                Vector3 normalVelocity1 = Vector3.Dot(normalTri, particles[tri.indexTriA].Velocity) * normalTri;
                                Vector3 normalVelocity2 = Vector3.Dot(normalTri, particles[tri.indexTriB].Velocity) * normalTri;
                                Vector3 normalVelocity3 = Vector3.Dot(normalTri, particles[tri.indexTriC].Velocity) * normalTri;

                                Vector3 tangencialVelocity0 = particles[idx].Velocity - normalVelocity0;
                                Vector3 tangencialVelocity1 = particles[tri.indexTriA].Velocity - normalVelocity0;
                                Vector3 tangencialVelocity2 = particles[tri.indexTriB].Velocity - normalVelocity0;
                                Vector3 tangencialVelocity3 = particles[tri.indexTriC].Velocity - normalVelocity0;

                                //particles[idx].AddForce(((3f*w/4.0f) * (-normalVelocity0 - particles[idx].Velocity))/0.02f); //* ( (tangencialVelocity0 - 0.6f * normalVelocity0.magnitude * (tangencialVelocity0/tangencialVelocity0.magnitude) - normalVelocity0) - particles[idx].Velocity));
                                //particles[tri.indexTriA].AddForce(((3f*w/4.0f) * (normalVelocity1 - particles[tri.indexTriA].Velocity))/0.02f);
                                //particles[tri.indexTriB].AddForce(((3f*w/4.0f) * (normalVelocity2 - particles[tri.indexTriB].Velocity))/0.02f);
                                //particles[tri.indexTriC].AddForce(((3f*w/4.0f) * (normalVelocity3 - particles[tri.indexTriC].Velocity))/0.02f);
                                

                                particles[idx].Velocity = Vector3.zero;//-normalVelocity0; //tangencialVelocity0 - 0.5f * normalVelocity0.magnitude * (tangencialVelocity0/tangencialVelocity0.magnitude) - normalVelocity0;//corr;
                                particles[tri.indexTriA].Velocity = Vector3.zero;//-normalVelocity1;//1.0f * (tangencialVelocity1 - 0.5f * normalVelocity1.magnitude * (tangencialVelocity1/tangencialVelocity1.magnitude) - normalVelocity1);
                                particles[tri.indexTriB].Velocity = Vector3.zero;//-normalVelocity2;//1.0f * (tangencialVelocity2 - 0.5f * normalVelocity2.magnitude * (tangencialVelocity2/tangencialVelocity2.magnitude) - normalVelocity2);
                                particles[tri.indexTriC].Velocity = Vector3.zero;//-normalVelocity3;//1.0f * (tangencialVelocity3 - 0.5f * normalVelocity3.magnitude * (tangencialVelocity3/tangencialVelocity3.magnitude) - normalVelocity3);                                

                                //Force
                                Vector3 normalForce0 = Vector3.Dot(normalTri, particles[idx].Force) * normalTri;
                                Vector3 normalForce1 = Vector3.Dot(normalTri, particles[tri.indexTriA].Force) * normalTri;
                                Vector3 normalForce2 = Vector3.Dot(normalTri, particles[tri.indexTriB].Force) * normalTri;
                                Vector3 normalForce3 = Vector3.Dot(normalTri, particles[tri.indexTriC].Force) * normalTri;

                                Vector3 tangencialForce0 = particles[idx].Force - normalForce0;
                                Vector3 tangencialForce1 = particles[tri.indexTriA].Force - normalForce1;
                                Vector3 tangencialForce2 = particles[tri.indexTriB].Force - normalForce2;
                                Vector3 tangencialForce3 = particles[tri.indexTriC].Force - normalForce3;

                            
                                //particles[idx].AddForce(tangencialForce0 - 0.5f * normalForce0.magnitude * (tangencialForce0/tangencialForce0.magnitude));
                                //particles[tri.indexTriA].AddForce(1.0f * (tangencialForce1 - 0.5f * normalForce1.magnitude * (tangencialForce1/tangencialForce1.magnitude)));
                                //particles[tri.indexTriB].AddForce(1.0f * (tangencialForce2 - 0.5f * normalForce2.magnitude * (tangencialForce2/tangencialForce2.magnitude)));
                                //particles[tri.indexTriC].AddForce(1.0f * (tangencialForce3 - 0.5f * normalForce3.magnitude * (tangencialForce3/tangencialForce3.magnitude)));
                                
                            }   
                        }
                    }
                }
            }
        }
    }

    public bool TrianglePointDistance(
        Vector3 p, float w,
        Vector3 p0, float w0,
        Vector3 p1, float w1,
        Vector3 p2, float w2,
        float restDist,
        float compressionStiffness,
        float stretchStiffness,
        out Vector3 corr,
        out Vector3 corr0,
        out Vector3 corr1,
        out Vector3 corr2,
        out Vector3 normalTri,
        out float val0,
        out float val1,
        out float val2)
    {
        corr = corr0 = corr1 = corr2 = Vector3.zero;
        val0 = val1 = val2 = 0.0f;
        // find barycentric coordinates of closest point on triangle

        // for singular case
        float b0 = 1.0f / 3.0f;
        float b1 = b0;
        float b2 = b0;

        float a, b, c, d, e, f;
        float det;

        Vector3 d1 = p1 - p0;
        Vector3 d2 = p2 - p0;
        normalTri = Vector3.Cross(d1, d2);
        Vector3 pp0 = p - p0;
        a = Vector3.Dot(d1, d1);
        b = Vector3.Dot(d2, d1);
        c = Vector3.Dot(pp0, d1);
        d = b;
        e = Vector3.Dot(d2, d2);
        f = Vector3.Dot(pp0, d2);
        det = a*e - b*d;

        float s, t;
        Vector3 dist;
        float dist2;
        if (det != 0.0f)
        {
            s = (c*e - b*f) / det;
            t = (a*f - c*d) / det;
            // inside triangle
            b0 = 1.0f - s - t;
            b1 = s;
            b2 = t;
            if (b0 < 0.0)
            {
                // on edge 1-2
                dist = p2 - p1;
                dist2 = Vector3.Dot(dist, dist);
                t = (dist2 == 0.0f) ? 0.5f : Vector3.Dot(dist, (p - p1)) / dist2;
                if (t < 0.0) t = 0.0f;	// on point 1
                if (t > 1.0) t = 1.0f;	// on point 2
                b0 = 0.0f;
                b1 = (1.0f - t);
                b2 = t;
            }
            else if (b1 < 0.0)
            {
                // on edge 2-0
                dist = p0 - p2;
                dist2 = Vector3.Dot(dist, dist);
                t = (dist2 == 0.0f) ? 0.5f : Vector3.Dot(dist, (p - p2)) / dist2;
                if (t < 0.0) t = 0.0f;	// on point 2
                if (t > 1.0) t = 1.0f; // on point 0
                b1 = 0.0f;
                b2 = (1.0f - t);
                b0 = t;
            }
            else if (b2 < 0.0)
            {
                // on edge 0-1
                dist = p1 - p0;
                dist2 = Vector3.Dot(dist, dist);
                t = (dist2 == 0.0f) ? 0.5f : Vector3.Dot(dist, (p - p0)) / dist2;
                if (t < 0.0) t = 0.0f;	// on point 0
                if (t > 1.0) t = 1.0f;	// on point 1
                b2 = 0.0f;
                b0 = (1.0f - t);
                b1 = t;
            }
        }
        Vector3 q = p0 * b0 + p1 * b1 + p2 * b2;
        //normalTri = q;
        Vector3 n = p - q;
        float l = n.magnitude;
        Vector3.Normalize(n);
        float C = l - restDist;
        Vector3 grad = n;
        Vector3 grad0 = -n * b0;
        Vector3 grad1 = -n * b1;
        Vector3 grad2 = -n * b2;

        s = w + w0 * b0*b0 + w1 * b1*b1 + w2 * b2*b2;
        if (s == 0.0f)
        return false;

        s = C / s;
        if (C < 0.0f) s *= compressionStiffness;
        else s *= stretchStiffness;

        if (s == 0.0f) return false;

        corr = -s * w * grad;
        corr0 = -s * w0 * grad0;
        corr1 = -s * w1 * grad1;
        corr2 = -s * w2 * grad2;
        /*corr = -n;
        corr0 = n * b0;
        corr1 = n * b1;
        corr2 = n * b2;*/

        val0 = b0;
        val1 = b1;
        val2 = b2;

        return true;
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

