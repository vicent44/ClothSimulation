﻿using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;

public class Simulate
{
    List<Particles> particles;
    List<Springs> springs;
    List<Triangles> triangles;
    Vector3 winddirectiondensity;
    Transform plane;
    Transform secondPlane;
    Vector3 directionNormalSecondPlane;

    public struct SPHash
    {
        public List<int> indices;
    }
    Hashing hash;
    public int gridSize;
    public float frictionConstPlane;
    public float dissipationConstPlane;
    public float frictionConstCloth;
    public float dissipationConstCloth;
    public bool drawSprings;

    public Simulate(List<Particles> particles, List<Springs> springs, List<Triangles> triangles, Vector3 winddirectiondensity, Transform plane, Transform secondPlane, Vector3 directionNormalSecondPlane, int gridSize, float frictionConstPlane, float dissipationConstPlane, float frictionConstCloth, float dissipationConstCloth, bool drawSprings)
    {
        this.particles = particles;
        this.springs = springs;
        this.triangles = triangles;
        this.winddirectiondensity = winddirectiondensity;
        this.plane = plane;
        this.secondPlane = secondPlane;
        this.directionNormalSecondPlane = directionNormalSecondPlane;
        this.gridSize = gridSize;
        this.frictionConstPlane = frictionConstPlane;
        this.dissipationConstPlane = dissipationConstPlane;
        this.frictionConstCloth = frictionConstCloth;
        this.dissipationConstCloth = dissipationConstCloth;
        this.drawSprings = drawSprings;
    }

    public void Update(float dt, List<Triangles> triangles)
    {
        this.triangles = triangles; //Update new triangles for the windForce calculation
        ComputeTotalForces();
        WindForce();
        IntegratorVerlet(dt);
        CheckPlaneCollitions(dt, frictionConstPlane, dissipationConstPlane);
        CheckSelfCollitions(dt, frictionConstCloth, dissipationConstCloth);
    }

    public void ComputeTotalForces()
    {
        //First add the gravity force
        foreach( var p in particles)
        {
            p.ResetResultantForce();
            if(p.isActive){
                var gravity = new Vector3(0, -9.8f, 0) * p.Mass;
                p.AddForce(gravity);
            }
        }
        //Second add the hook force and damping force
        foreach (var s in springs)
        {
            s.ApplyForce();
        }
    }
    //Euler integrator
    void IntegratorEuler(float dt)
    {
        foreach(var p in particles)
        {
            float deltaTimeMass = dt / p.Mass;
            p.Velocity += p.Force * deltaTimeMass;
            p.Position += p.Velocity * dt;
        }
    }
    //Verlet integrator called for each particle the verlet integrator
    void IntegratorVerlet(float dt)
    {
        foreach(var p in particles)
        {
            if(p.isActive) p.UpdateParticle(dt);
        }
    }

    void CheckPlaneCollitions(float dt, float frictionConstPlane, float dissipationConstPlane)
    {
        foreach(var p in particles)
        {
            Vector3 d = p.Position - (plane.position + new Vector3(0f,0.02f,0f));
            Vector3 normalPlane = new Vector3(0f,1f,0f);

            float dot = Vector3.Dot(normalPlane, d);
            if(dot <= 0)
            {
                p.Position -= dot * normalPlane;
                Vector3 normalVelocity = Vector3.Dot(normalPlane,p.Velocity) * normalPlane;
                Vector3 tangencialVelocity = p.Velocity - normalVelocity;
                Vector3 normalForce = Vector3.Dot(p.Force, normalPlane) * normalPlane;
                Vector3 tangencialForce = p.Force - normalForce;

                p.Position = p.Position - dt * (tangencialVelocity-frictionConstPlane*normalVelocity.magnitude*(tangencialVelocity/tangencialVelocity.magnitude)-dissipationConstPlane*normalVelocity);
            }
        }
    }

    void CheckSecondPlaneCollitions(float dt, float frictionConstPlane, float dissipationConstPlane)
    {
        foreach(var p in particles)
        {
            Vector3 normalPlane = directionNormalSecondPlane;
            Vector3 d = p.Position - (secondPlane.position + 0.02f * normalPlane);

            float dot = Vector3.Dot(normalPlane, d);
            if(dot <= 0)
            {
                p.Position -= dot * normalPlane;
                Vector3 normalVelocity = Vector3.Dot(normalPlane,p.Velocity) * normalPlane;
                Vector3 tangencialVelocity = p.Velocity - normalVelocity;
                Vector3 normalForce = Vector3.Dot(p.Force, normalPlane) * normalPlane;
                Vector3 tangencialForce = p.Force - normalForce;

                p.Position = p.Position - dt * (tangencialVelocity-frictionConstPlane*normalVelocity.magnitude*(tangencialVelocity/tangencialVelocity.magnitude)-dissipationConstPlane*normalVelocity);
            }
        }
    }

    void CheckSelfCollitions(float dt, float frictionConstCloth, float dissipationConstCloth)
    {
        hash = new Hashing(gridSize, 1.0f/(float)gridSize, 1723);
        SPHash[] spHash = new SPHash[1723];

        for(int v = 0; v < particles.Count; v++)
        {
            int has = Mathf.Abs(hash.Hash(particles[v].Prev));

            if(spHash[has].indices == null) spHash[has].indices = new List<int>();
            spHash[has].indices.Add(particles[v].I);
        }

        for(int t = 0; t < triangles.Count; t++)
        {
            var tri = triangles[t];
            var p0 = particles[tri.indexTriA].Prev;
            var p1 = particles[tri.indexTriB].Prev;
            var p2 = particles[tri.indexTriC].Prev;

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
                            Vector3 p = particles[idx].Prev;
                            float w = 1f/particles[idx].Mass;

                            Vector3 corr, corr0, corr1, corr2, normalTri;
                            float val0, val1, val2;
                            if(TrianglePointDistance(
                                p, w,
                                p0, w0,
                                p1, w1,
                                p2, w2,
                                0.05f, 500f, 0.0f,
                                out corr, out corr0, out corr1, out corr2, out normalTri,
                                out val0, out val1, out val2))
                            {
                                //Debug.Log(particles[idx].Position);
                                //particles[idx].Position += corr;
                                //particles[tri.indexTriA].Position += corr0;
                                //particles[tri.indexTriB].Position += corr1;
                                //particles[tri.indexTriC].Position += corr2;

                                /*Vector3 vel = particles[idx].Velocity + corr;
                                Vector3 vel0 = particles[tri.indexTriA].Velocity + corr0;
                                Vector3 vel1 = particles[tri.indexTriB].Velocity + corr1;
                                Vector3 vel2 = particles[tri.indexTriC].Velocity + corr2;
                                
                                //Position
                                particles[idx].Position = particles[idx].Prev + 0.02f * vel;
                                particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev + 0.02f * vel0; 
                                particles[tri.indexTriB].Position = particles[tri.indexTriB].Prev + 0.02f * vel1;
                                particles[tri.indexTriC].Position = particles[tri.indexTriC].Prev + 0.02f * vel2;*/

                                Vector3 d0 = particles[idx].Position - corr;
                                float dot0 = Vector3.Dot(normalTri, d0);
                                //particles[idx].Position -= dot0 * normalTri;
                                //if(dot0<0)particles[idx].AddPosition(-dot0 * normalTri); 
                                //if(dot0>0)particles[idx].AddPosition(dot0 * normalTri);

                                //Velocity
                                Vector3 normalVelocity0 = Vector3.Dot(normalTri, particles[idx].Velocity) * normalTri;
                                Vector3 normalVelocity1 = Vector3.Dot(normalTri, particles[tri.indexTriA].Velocity) * normalTri;
                                Vector3 normalVelocity2 = Vector3.Dot(normalTri, particles[tri.indexTriB].Velocity) * normalTri;
                                Vector3 normalVelocity3 = Vector3.Dot(normalTri, particles[tri.indexTriC].Velocity) * normalTri;
                                
                                /*Vector3 normalVelocity0 = Vector3.Dot(corr0, particles[idx].Velocity) * corr0;
                                Vector3 normalVelocity1 = Vector3.Dot(corr0, particles[tri.indexTriA].Velocity) * corr0;
                                Vector3 normalVelocity2 = Vector3.Dot(corr0, particles[tri.indexTriB].Velocity) * corr0;
                                Vector3 normalVelocity3 = Vector3.Dot(corr0, particles[tri.indexTriC].Velocity) * corr0;*/
                                
                                Vector3 tangencialVelocity0 = particles[idx].Velocity - normalVelocity0;
                                Vector3 tangencialVelocity1 = particles[tri.indexTriA].Velocity - normalVelocity0;
                                Vector3 tangencialVelocity2 = particles[tri.indexTriB].Velocity - normalVelocity0;
                                Vector3 tangencialVelocity3 = particles[tri.indexTriC].Velocity - normalVelocity0;

                                Vector3 velocity0 = - normalVelocity0;
                                Vector3 velocity1 = - normalVelocity1;
                                Vector3 velocity2 = - normalVelocity2;
                                Vector3 velocity3 = - normalVelocity3;

                                Vector3 velocityBest0 = tangencialVelocity0 - frictionConstCloth*normalVelocity0.magnitude*(tangencialVelocity0/tangencialVelocity0.magnitude) - dissipationConstCloth*normalVelocity0;
                                Vector3 velocityBest1 = tangencialVelocity1 - frictionConstCloth*normalVelocity1.magnitude*(tangencialVelocity1/tangencialVelocity1.magnitude) - dissipationConstCloth*normalVelocity1;
                                Vector3 velocityBest2 = tangencialVelocity2 - frictionConstCloth*normalVelocity2.magnitude*(tangencialVelocity2/tangencialVelocity2.magnitude) - dissipationConstCloth*normalVelocity2;
                                Vector3 velocityBest3 = tangencialVelocity3 - frictionConstCloth*normalVelocity3.magnitude*(tangencialVelocity3/tangencialVelocity3.magnitude) - dissipationConstCloth*normalVelocity3;

                                //Position
                                //particles[idx].AddPosition(-dt * velocityBest0);
                                particles[idx].Position = particles[idx].Prev + dt * velocityBest0;
                                //particles[tri.indexTriA].AddPosition(- dt * velocityBest1 * val0);
                                particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev - dt * velocityBest1 * val0;
                                //particles[tri.indexTriB].AddPosition(- dt * velocityBest2 * val1);
                                particles[tri.indexTriB].Position = particles[tri.indexTriB].Prev - dt * velocityBest2 * val1;
                                //particles[tri.indexTriC].AddPosition(- dt * velocityBest3 * val2);
                                particles[tri.indexTriC].Position = particles[tri.indexTriC].Prev - dt * velocityBest3 * val2;


                                /*particles[idx].Position = particles[idx].Prev - normalVelocity0 * 0.02f;
                                particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev + normalVelocity1 * 0.02f; 
                                particles[tri.indexTriB].Position = particles[tri.indexTriA].Prev + normalVelocity2 * 0.02f;
                                particles[tri.indexTriC].Position = particles[tri.indexTriA].Prev + normalVelocity3 * 0.02f;*/

                                //particles[idx].Velocity = -normalVelocity0; //tangencialVelocity0 - 0.5f * normalVelocity0.magnitude * (tangencialVelocity0/tangencialVelocity0.magnitude) - normalVelocity0;//corr;
                                //particles[tri.indexTriA].Velocity = -normalVelocity1;//1.0f * (tangencialVelocity1 - 0.5f * normalVelocity1.magnitude * (tangencialVelocity1/tangencialVelocity1.magnitude) - normalVelocity1);
                                //particles[tri.indexTriB].Velocity = -normalVelocity2;//1.0f * (tangencialVelocity2 - 0.5f * normalVelocity2.magnitude * (tangencialVelocity2/tangencialVelocity2.magnitude) - normalVelocity2);
                                //particles[tri.indexTriC].Velocity = -normalVelocity3;//1.0f * (tangencialVelocity3 - 0.5f * normalVelocity3.magnitude * (tangencialVelocity3/tangencialVelocity3.magnitude) - normalVelocity3);                                

                                //Force
                                Vector3 normalForce0 = Vector3.Dot(normalTri, particles[idx].Force) * normalTri;
                                Vector3 normalForce1 = Vector3.Dot(normalTri, particles[tri.indexTriA].Force) * normalTri;
                                Vector3 normalForce2 = Vector3.Dot(normalTri, particles[tri.indexTriB].Force) * normalTri;
                                Vector3 normalForce3 = Vector3.Dot(normalTri, particles[tri.indexTriC].Force) * normalTri;

                                Vector3 tangencialForce0 = particles[idx].Force - normalForce0;
                                Vector3 tangencialForce1 = particles[tri.indexTriA].Force - normalForce1;
                                Vector3 tangencialForce2 = particles[tri.indexTriB].Force - normalForce2;
                                Vector3 tangencialForce3 = particles[tri.indexTriC].Force - normalForce3;
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
        Vector3 q = p0 * b0 + p1 * b1 + p2 * b2; //Point of baricentric

        Vector3 n = p - q; //Distance between point baricentric and node
        float l = n.magnitude; //Distance in magnitude
        //Vector3.Normalize(n); //Direction where the point is from the surface of the triangle
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

        corr = q;//-s * w * grad;
        corr0 = n;//-s * w0 * grad0;
        corr1 = -s * w1 * grad1;
        corr2 = -s * w2 * grad2;

        val0 = b0;
        val1 = b1;
        val2 = b2;

        return true;
    }

    //Adding the wind force taking in account every triangle area
    void WindForce()
    {
        foreach(var t in triangles)
        {
            Vector3 inter = Vector3.Cross(t.posTriA - t.posTriC, t.posTriB - t.posTriC); //Calculate the area of square
            float areatriangle = 0.5f * (float)Math.Sqrt(Vector3.Dot(inter, inter)); //Area triangles is half than square
            Vector3 forceTriangle = areatriangle * winddirectiondensity; //Finally compute the total forçe once we know the area
            //Add one third of the force to each node
            particles[t.indexTriA].AddForce(forceTriangle/3.0f);
            particles[t.indexTriB].AddForce(forceTriangle/3.0f);
            particles[t.indexTriC].AddForce(forceTriangle/3.0f);
        }
    }
    //Draw the lines between particles
    public void DrawGizmos()
    {
        if(drawSprings)
        {
            for(int i = 0, n = springs.Count; i < n; i++)
            {
                var p = springs[i];
                Gizmos.color = Color.white;
                Gizmos.DrawLine(p.particleA.Position, p.particleB.Position);
            }
        }
    }
}