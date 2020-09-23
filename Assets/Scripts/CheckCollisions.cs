using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;

public class CheckCollisions
{

    List<Particles> particles;
    List<Springs> springs;
    List<Triangles> triangles;
    List<Edges> edges;
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

    public CheckCollisions(float dt, float frictionConstCloth, float dissipationConstCloth)
    {

    }

    /*public Update(float dt, List<Triangles> triangles_p, List<Edges> edges_p, List<Particles> particles_p, List<Triangles> triangles_c, List<Edges> edges_c, List<Particles> particles_c)
    {

    }*/

    //Check self collisions edge-triangle
    void CheckSelfCollitionsEdge(float dt, float frictionConstCloth, float dissipationConstCloth)
    {
        //First I create the hesh table.
        hash = new Hashing(gridSize, 1.0f/(float)gridSize, 1723);
        SPHash[] spHash = new SPHash[1723];

        //Second I put the edges on a minimun boxes and
        //I put it into the hash table
        for(int e = 0; e < edges.Count; e++)
        {
            var ed = edges[e];
            var p0 = particles[ed.indexEdgeA].Position;
            var p1 = particles[ed.indexEdgeB].Position;

            int has = hash.LineBoxHashes(p0, p1);

            if(spHash[has].indices == null) spHash[has].indices = new List<int>();
            spHash[has].indices.Add(ed.indexEdge);
        }

        //Now I put every triangle in a hash table (well, I cehck
        //the hash number of it) and if in this place on the edge
        //table there are any edge I check if there is a collition
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
                        var a0 = particles[edges[idx].indexEdgeA].I;
                        var a1 = particles[edges[idx].indexEdgeB].I;

                        var b0 = particles[tri.indexTriA].I;
                        var b1 = particles[tri.indexTriB].I;
                        var b2 = particles[tri.indexTriC].I;
                        //Here to check that is not the same edge
                        if(a0 != b0 && a0 != b1 && a0 != b2 && a1 != b0 && a1 != b1 && a1 != b2)
                        {
                            var ed = edges[idx];
                            var u0 = particles[ed.indexEdgeA].Position;
                            var u1 = particles[ed.indexEdgeB].Position;

                            Vector3 out0, out1, out2, out3, normalTri;
                            Vector3 corrP1, corrP2, corrT1, corrT2, corrT3;
                            int test;
                            int situation;
                            float valP, val1, val2, val3;
                            //Finally check if there is a collision
                            if(EdgeTriangleIntersect(u0, u1, p0, p1, p2, 0.2f, 50f, 0f,out corrP1, out corrP2, out corrT1, out corrT2, out corrT3, out normalTri, out situation,
                                out valP, out val1, out val2, out val3))
                            {
                                //Now I know that there is a collision so I have to
                                //make a correction to avoid that collision.
                                //The technique is the same than in the plane case
                                //but here the plane is the triangle and its normal.
                                //However here is posible for a particle to stay on the
                                //negative part of the trignel normal. So in edge case,
                                //if there is a collision I check where the collision
                                //have just happened if is on the 2 extrems (the points)
                                //just no correction however if the edge is cut I check
                                //where have happened and I where the the short distance
                                //penetrated is I correct that point because it have to
                                //be the last that have just entered.
                                //Then I do a projection depending on the side of the plane
                                //the dot product will be positive or negative.
                                //After that I modify the velocity to be able to have
                                //inelastic collisions.

                                //Position
                                Vector3 dA = particles[ed.indexEdgeA].Position - (corrT1 + 0.02f*normalTri);
                                Vector3 dB = particles[ed.indexEdgeB].Position - (corrT1 + 0.02f*normalTri);

                                float dotA = Vector3.Dot(normalTri, dA);
                                float dotB = Vector3.Dot(normalTri, dB);
                                
                                if(situation == 3)
                                {
                                    if(dotA < 0f)
                                    {
                                        particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Prev - normalTri*dotA;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev + normalTri*dotA * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Prev + normalTri*dotA * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Prev + normalTri*dotA * val3;                                    
                                    }
                                    if(dotA > 0f)
                                    {
                                        particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Prev + normalTri*dotA;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev - normalTri*dotA * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Prev - normalTri*dotA * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Prev - normalTri*dotA * val3;                                    
                                    }
                                }
                                if(situation == 4)
                                {
                                    if(dotB < 0f)
                                    {
                                        particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Prev - normalTri*dotB;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev + normalTri*dotB * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Prev + normalTri*dotB * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Prev + normalTri*dotB * val3;
                                    }
                                    if(dotB > 0f)
                                    {
                                        particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Prev + normalTri*dotB;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev - normalTri*dotB * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Prev - normalTri*dotB * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Prev - normalTri*dotB * val3;
                                    }
                                }

                                //Velocity
                                Vector3 normalVelocityA = Vector3.Dot(normalTri, particles[ed.indexEdgeA].Velocity) * normalTri;
                                Vector3 normalVelocityB = Vector3.Dot(normalTri, particles[ed.indexEdgeB].Velocity) * normalTri;
                                Vector3 normalVelocity1 = Vector3.Dot(normalTri, particles[tri.indexTriA].Velocity) * normalTri;
                                Vector3 normalVelocity2 = Vector3.Dot(normalTri, particles[tri.indexTriB].Velocity) * normalTri;
                                Vector3 normalVelocity3 = Vector3.Dot(normalTri, particles[tri.indexTriC].Velocity) * normalTri;
                                
                                Vector3 tangencialVelocityA = particles[ed.indexEdgeA].Velocity - normalVelocityA;
                                Vector3 tangencialVelocityB = particles[ed.indexEdgeB].Velocity - normalVelocityB;
                                Vector3 tangencialVelocity1 = particles[tri.indexTriA].Velocity - normalVelocity1;
                                Vector3 tangencialVelocity2 = particles[tri.indexTriB].Velocity - normalVelocity2;
                                Vector3 tangencialVelocity3 = particles[tri.indexTriC].Velocity - normalVelocity3;

                                Vector3 velocityBestA = tangencialVelocityA - frictionConstCloth*normalVelocityA.magnitude*(tangencialVelocityA/tangencialVelocityA.magnitude) - dissipationConstCloth*normalVelocityA;
                                Vector3 velocityBestB = tangencialVelocityB - frictionConstCloth*normalVelocityB.magnitude*(tangencialVelocityB/tangencialVelocityB.magnitude) - dissipationConstCloth*normalVelocityB;
                                Vector3 velocityBest1 = tangencialVelocity1 - frictionConstCloth*normalVelocity1.magnitude*(tangencialVelocity1/tangencialVelocity1.magnitude) - dissipationConstCloth*normalVelocity1;
                                Vector3 velocityBest2 = tangencialVelocity2 - frictionConstCloth*normalVelocity2.magnitude*(tangencialVelocity2/tangencialVelocity2.magnitude) - dissipationConstCloth*normalVelocity2;
                                Vector3 velocityBest3 = tangencialVelocity3 - frictionConstCloth*normalVelocity3.magnitude*(tangencialVelocity3/tangencialVelocity3.magnitude) - dissipationConstCloth*normalVelocity3;
                                
                                /*if(situation == 3)
                                {
                                    if(dotA < 0f)
                                    {
                                        particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Position - dt * velocityBestA;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Position + dt * velocityBest1 * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Position + dt * velocityBest2 * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Position + dt * velocityBest3 * val3;                                    
                                    }
                                    if(dotA > 0f)
                                    {
                                        particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Position + dt * velocityBestA;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Position - dt * velocityBest1 * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Position - dt * velocityBest2 * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Position - dt * velocityBest3 * val3;                                    
                                    }
                                }
                                if(situation == 4)
                                {
                                    if(dotB < 0f)
                                    {
                                        particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Position - dt * velocityBestB;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Position + dt * velocityBest1 * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Position + dt * velocityBest2 * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Position + dt * velocityBest3 * val3;
                                    }
                                    if(dotB > 0f)
                                    {
                                        particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Position + dt * velocityBestB;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Position - dt * velocityBest1 * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Position - dt * velocityBest2 * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Position - dt * velocityBest3 * val3;
                                    }
                                }*/                               
                            }   
                        }
                    }
                }
            }
        }
    }

    //Check self collisions point-triangle
    void CheckSelfCollitionsPoint(float dt, float frictionConstCloth, float dissipationConstCloth)
    {
        //Here I use the same metodology than before (edge-triangle collisions)
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

            List<int> hashes = hash.TriangleBoundingBoxHashes(p0, p1, p2);
    
            for(int h = 0; h < hashes.Count; h++)
            {
                if(spHash[h].indices != null)
                {
                    for(int sph = 0; sph < spHash[h].indices.Count; sph++)
                    {
                        int idx = spHash[h].indices[sph];
                        //Here check that the point don't takes part of the triangle
                        if(idx != particles[tri.indexTriA].I && idx != particles[tri.indexTriB].I && idx != particles[tri.indexTriC].I)
                        {
                            Vector3 p = particles[idx].Position;

                            Vector3 corr, corr0, corr1, corr2, normalTri;
                            float val0, val1, val2;
                            if(TrianglePointDistance(
                                p,
                                p0,
                                p1,
                                p2,
                                0.02f, 500f, 0.0f,
                                out corr, out corr0, out corr1, out corr2, out normalTri,
                                out val0, out val1, out val2))
                            {
                                //Here same metodology than plane check collisions
                                //but here is posible to have the particle on the 
                                //negative triangle normal plane. If the particle
                                //has to be on the negative side correct the position
                                //and the other way rround. Also correction velocity
                                //to simulate inelastic collitions.
                                Vector3 d0 = particles[idx].Position - corr;
                                float dot0 = Vector3.Dot(normalTri, d0);
                                if(dot0 < 0f)
                                {
                                    particles[idx].Position = particles[idx].Prev - dot0 * normalTri;
                                    particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev + dot0 * normalTri * val0; 
                                    particles[tri.indexTriB].Position = particles[tri.indexTriB].Prev + dot0 * normalTri * val1;
                                    particles[tri.indexTriC].Position = particles[tri.indexTriC].Prev + dot0 * normalTri * val2;
                                }

                                if(dot0 > 0f)
                                {
                                    particles[idx].Position = particles[idx].Prev + dot0 * normalTri;
                                    particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev - dot0 * normalTri * val0; 
                                    particles[tri.indexTriB].Position = particles[tri.indexTriB].Prev - dot0 * normalTri * val1;
                                    particles[tri.indexTriC].Position = particles[tri.indexTriC].Prev - dot0 * normalTri * val2;
                                }

                                //Velocity
                                Vector3 normalVelocity0 = Vector3.Dot(normalTri, particles[idx].Velocity) * normalTri;
                                Vector3 normalVelocity1 = Vector3.Dot(normalTri, particles[tri.indexTriA].Velocity) * normalTri;
                                Vector3 normalVelocity2 = Vector3.Dot(normalTri, particles[tri.indexTriB].Velocity) * normalTri;
                                Vector3 normalVelocity3 = Vector3.Dot(normalTri, particles[tri.indexTriC].Velocity) * normalTri;
                                
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

                                //Velocity

                                /*if(dot0 < 0f)
                                {
                                    particles[idx].Position = particles[idx].Position - dt * velocityBest0;
                                    particles[tri.indexTriA].Position = particles[tri.indexTriA].Position + dt * velocityBest1 * val0; 
                                    particles[tri.indexTriB].Position = particles[tri.indexTriB].Position + dt * velocityBest2 * val1;
                                    particles[tri.indexTriC].Position = particles[tri.indexTriC].Position + dt * velocityBest3 * val2;
                                }

                                if(dot0 > 0f)
                                {
                                    particles[idx].Position = particles[idx].Position + dt * velocityBest0;
                                    particles[tri.indexTriA].Position = particles[tri.indexTriA].Position - dt * velocityBest1 * val0; 
                                    particles[tri.indexTriB].Position = particles[tri.indexTriB].Position - dt * velocityBest2 * val1;
                                    particles[tri.indexTriC].Position = particles[tri.indexTriC].Position - dt * velocityBest3 * val2;
                                }*/
                            }   
                        }
                    }
                }
            }
        }
    }

    public bool TrianglePointDistance(
        Vector3 p,
        Vector3 p0,
        Vector3 p1,
        Vector3 p2, float restDist, float compressionStiffness, float stretchStiffness,
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
        Vector3.Normalize(n); //Direction where the point is from the surface of the triangle
        float C = l - restDist;
        Vector3 grad = n;
        Vector3 grad0 = -n * b0;
        Vector3 grad1 = -n * b1;
        Vector3 grad2 = -n * b2;

        s = 1f + b0*b0 + b1*b1 + b2*b2;
        if (s == 0.0f)
        return false;

        s = C / s;
        if (C < 0.0f) s *= compressionStiffness;
        else s *= stretchStiffness;

        if (s == 0.0f) return false;

        corr = q;//-s * w * grad;
        corr0 = n;//-s * w0 * grad0;
        corr1 = -s * grad1;
        corr2 = -s * grad2;

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

    public bool EdgeTriangleIntersect( 
        Vector3 orig, Vector3 fin, 
        Vector3 v0, Vector3 v1, Vector3 v2, float restDist, float compressionStiffness, float stretchStiffness,
        out Vector3 corrP1, out Vector3 corrP2, 
        out Vector3 corrT1, out Vector3 corrT2, out Vector3 corrT3,
        out Vector3 normalTri, out int situation,
        out float valP, out float val1, out float val2, out float val3) 
    { 
        corrP1 = corrP2 = corrT1 = corrT2 = corrT3 = normalTri = Vector3.zero;
        //float t, u, v, uv;
        valP = val1 = val2 = val3 = 0f;
        situation = 0;
        float eps = 1e-6f;

        //Calculte the baricentric point of the triangle and
        //calculate the point of collision.
        Vector3 edge1 = v1 - v0;
        Vector3 edge2 = v2 - v0;
        normalTri = Vector3.Cross(edge1, edge2);
        //Vector3 n = Vector3.Cross(edge1, edge2);
        Vector3 dir = fin - orig;
        Vector3 h = Vector3.Cross(dir, edge2);
        float a = Vector3.Dot(edge1, h);

        if(a > -eps && a < eps) return false;
        float f = 1f/a;
        Vector3 s = orig - v0;
        float u = f * Vector3.Dot(s, h);

        if(u < 0f || u > 1f) return false;
        Vector3 q = Vector3.Cross(s, edge1);
        float v = f * Vector3.Dot(dir, q);

        if(v < 0f || u + v > 1f) return false;
        float t = f * Vector3.Dot(edge2, q);

        if(t < 0f || t > 1f) return false;

        float uv = 1f - u - v;
        Vector3 p = uv * v0 + u * v1 + v * v2;
        Vector3 e = orig + t * (fin - orig);
        Vector3 nn = e - p;
        float dist = nn.magnitude;
        Vector3.Normalize(nn);
        float C = dist - restDist;
        Vector3 gradA = nn * u;
        Vector3 gradB = nn * v;
        Vector3 grad1 = -nn * u;
        Vector3 grad2 = -nn * v;
        Vector3 grad3 = -nn * uv;

        float se = u*u + v*v + uv*uv;

        //Here I check where the collision has done,
        //if the collisions is between any of the
        //vertices of the edge just nothing, there
        //isn't penetration but if the collision
        //is detected on the edge then check
        //the short side of the edge and push the
        //edge to that side (the limit of the sides
        //is the center of the edge).
        if(se == 0f) return false;

            Vector3 vvv = e - p;
            float disti = Vector3.Dot(vvv, normalTri);
            Vector3 projectedpoint = e - disti*normalTri;

        //Point triangle collision (nothing done)
        if(orig == e)
        {
            situation = 1;
            Vector3 normOrig = fin - orig;
            float saberRespecteNormalTri = Vector3.Dot(normOrig, normalTri);
            if(saberRespecteNormalTri > 0f) //el edge ve desde la normal
            {
                //Sha daplicar una forca a la particula
                // orig en la direccio normal del triangle
                //per tirarla fora
                corrP1 = Vector3.Normalize(normalTri) * normOrig.magnitude;
            }
            else
            {
                //Sha d'aplicar una froca e al particula
                // orig en la direccio contraria de la normal
                //del triangle per tirarla fora
                corrP1 = -Vector3.Normalize(normalTri) * normOrig.magnitude;
            }
        }
        if(fin == e)
        {
            situation = 2;
            Vector3 normFin = orig - fin;
            float saberRespecteNormalTri = Vector3.Dot(normFin, normalTri);
            if(saberRespecteNormalTri > 0f) //el edge ve desde la normal
            {
                corrP2 = Vector3.Normalize(normalTri) * normFin.magnitude;
            }
            else
            {
                corrP2 = -Vector3.Normalize(normalTri) * normFin.magnitude;
            }
        }
        //Busquem el centre del segement
        Vector3 centerEdge = (fin + orig)/2f;
        Vector3 tri1 = centerEdge - orig;
        Vector3 tri2 = centerEdge - fin;
        if(tri1.magnitude < tri2.magnitude)
        {
            situation = 3;
        }
        if(tri2.magnitude < tri1.magnitude)
        {
            situation = 4;
        }

        //corrP1 = -se * gradA;
        //corrP2 = -se * gradB;
        corrT1 = p;//-se * grad1;
        corrT2 = projectedpoint;//-se * grad2;
        corrT3 = -se * grad3;
        valP = t;
        val1 = uv;
        val2 = u;
        val3 = v;
    
        return true; // this ray hits the triangle 
    }

}
