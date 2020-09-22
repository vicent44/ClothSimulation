using System.Collections;
using System.Collections.Generic;
using UnityEngine;

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


    void CheckSelfCollitions(float dt, float frictionConstCloth, float dissipationConstCloth)
    {
        hash = new Hashing(gridSize, 1.0f/(float)gridSize, 1723);
        SPHash[] spHash = new SPHash[1723];

        /*for(int v = 0; v < particles.Count; v++)
        {
            int has = Mathf.Abs(hash.Hash(particles[v].Position));

            if(spHash[has].indices == null) spHash[has].indices = new List<int>();
            spHash[has].indices.Add(particles[v].I);
        }*/

        /*for(int t = 0; t < triangles.Count; t++)
        {
            var tri = triangles[t];
            var p0 = particles[tri.indexTriA].Position;
            var p1 = particles[tri.indexTriB].Position;
            var p2 = particles[tri.indexTriC].Position;

            int has = hash.TriangleBoxHashes(p0, p1, p2);

            if(spHash[has].indices == null) spHash[has].indices = new List<int>();
            spHash[has].indices.Add(tri.indexTriangle);
        }*/

        for(int e = 0; e < edges.Count; e++)
        {
            var ed = edges[e];
            var p0 = particles[ed.indexEdgeA].Position;
            var p1 = particles[ed.indexEdgeB].Position;

            int has = hash.LineBoxHashes(p0, p1);

            if(spHash[has].indices == null) spHash[has].indices = new List<int>();
            spHash[has].indices.Add(ed.indexEdge);
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
            //Debug.Log(hashes.Count);
            for(int h = 0; h < hashes.Count; h++)
            {
                if(spHash[h].indices != null)
                {
                    for(int sph = 0; sph < spHash[h].indices.Count; sph++)
                    {
                        int idx = spHash[h].indices[sph];
                        var a0 = particles[edges[idx].indexEdgeA].I;
                        var a1 = particles[edges[idx].indexEdgeB].I;

                        //var a2 = particles[triangles[idx].indexTriC].I;
                        var b0 = particles[tri.indexTriA].I;
                        var b1 = particles[tri.indexTriB].I;
                        var b2 = particles[tri.indexTriC].I;
                        if(a0 != b0 && a0 != b1 && a0 != b2 && a1 != b0 && a1 != b1 && a1 != b2)//idx != tri.indexTriangle)//particles[tri.indexTriA].I && idx != particles[tri.indexTriB].I && idx != particles[tri.indexTriC].I)
                        {
                            var ed = edges[idx];
                            var u0 = particles[ed.indexEdgeA].Position;
                            var u1 = particles[ed.indexEdgeB].Position;
                            //Vector3 p = particles[idx].Position;
                            //float w = 1f/particles[idx].Mass;
                            //var tr = triangles[idx];
                            //var u0 = particles[tr.indexTriA].Position;
                            //var u1 = particles[tr.indexTriB].Position;
                            //var u2 = particles[tr.indexTriC].Position; 
                            //Debug.Log(idx);
                            //Debug.Log(hashes[0]);
                            Vector3 out0, out1, out2, out3, normalTri;
                            Vector3 corrP1, corrP2, corrT1, corrT2, corrT3;
                            int test;
                            int situation;
                            float valP, val1, val2, val3;
                            //if(PointTriangleIntersect2(particles[ed.indexEdgeA], particles[ed.indexEdgeB], particles[tri.indexTriA], particles[tri.indexTriB], particles[tri.indexTriC], 0.2f, 50f, 0f,out corrP1, out corrP2, out corrT1, out corrT2, out corrT3, out normalTri, out situation,
                            //    out valP, out val1, out val2, out val3))
                            if(PointTriangleIntersect(u0, u1, p0, p1, p2, 0.2f, 50f, 0f,out corrP1, out corrP2, out corrT1, out corrT2, out corrT3, out normalTri, out situation,
                                out valP, out val1, out val2, out val3))
                            //if(EdgeforTriangle(u0, u1, p0, p1, p2, 0.05f, 50f, 0f, out out0, out out1, out out2, out out3, out test))
                            //if(TriTriIntersect(u0, u1, u2, p0, p1, p2))
                            //Vector3 corr, corr0, corr1, corr2, normalTri;
                            //float val0, val1, val2;
                            /*if(TrianglePointDistance(
                                p, w,
                                p0, w0,
                                p1, w1,
                                p2, w2,
                                0.05f, 500f, 0.0f,
                                out corr, out corr0, out corr1, out corr2, out normalTri,
                                out val0, out val1, out val2))*/
                            {
                                /*Vector3 corr, corr0, corr1, corr2, normalTri;
                                float val0, val1, val2;
                                TrianglePointDistance(p0, u0, u1, u2, 0f, 5f, 0f, out corr, out corr0, out corr1, out corr2, out normalTri,
                                out val0, out val1, out val2);
                                Vector3 n1 = corr0;
                                TrianglePointDistance(p1, u0, u1, u2, 0f, 5f, 0f, out corr, out corr0, out corr1, out corr2, out normalTri,
                                out val0, out val1, out val2);
                                Vector3 n2 = corr0;
                                TrianglePointDistance(p2, u0, u1, u2, 0f, 5f, 0f, out corr, out corr0, out corr1, out corr2, out normalTri,
                                out val0, out val1, out val2);
                                Vector3 n3 = corr0;*/



                                Vector3 dA = particles[ed.indexEdgeA].Position - (corrT1 + 0.02f*normalTri);
                                Vector3 dB = particles[ed.indexEdgeB].Position - (corrT1 + 0.02f*normalTri);

                                float dotA = Vector3.Dot(normalTri, dA);
                                float dotB = Vector3.Dot(normalTri, dB);
                                //particles[ed.indexEdgeA].Position -= dotA * normalTri;
                                //particles[ed.indexEdgeB].Position -= dotB * normalTri;
                                if(situation == 3)
                                {
                                    //if(dotA < 0f) particles[ed.indexEdgeA].Position -= dotA * normalTri;
                                    if(dotA < 0f)
                                    {
                                        particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Prev - normalTri*dotA;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev + normalTri*dotA * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Prev + normalTri*dotA * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Prev + normalTri*dotA * val3;                                    
                                    }
                                    //if(dotA > 0f) particles[ed.indexEdgeA].Position += dotA * normalTri;
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
                                    //if(dotB < 0f) particles[ed.indexEdgeB].Position -= dotB * normalTri;
                                    if(dotB < 0f)
                                    {
                                        particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Prev - normalTri*dotB;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev + normalTri*dotB * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Prev + normalTri*dotB * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Prev + normalTri*dotB * val3;
                                    }
                                    //if(dotB > 0f) particles[ed.indexEdgeB].Position += dotB * normalTri;
                                    if(dotB > 0f)
                                    {
                                        particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Prev + normalTri*dotB;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev - normalTri*dotB * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Prev - normalTri*dotB * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Prev - normalTri*dotB * val3;
                                    }
                                }

                                //particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Prev;
                                //particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Prev;
                                //particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev;
                                //particles[tri.indexTriB].Position = particles[tri.indexTriB].Prev;
                                //particles[tri.indexTriC].Position = particles[tri.indexTriC].Prev;                               


                                //Debug.Log("Collision");
                                //particles[ed.indexEdgeA].Position += corrP1;
                                //particles[ed.indexEdgeB].Position += corrP2;
                                //particles[tri.indexTriA].Position += corrT1;
                                //particles[tri.indexTriB].Position += corrT2;
                                //particles[tri.indexTriC].Position += corrT3;

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

                                //Vector3 d0 = particles[idx].Position - corr;
                                //float dot0 = Vector3.Dot(normalTri, d0);
                                //particles[idx].Position -= dot0 * normalTri;
                                //if(dot0<0)particles[idx].AddPosition(-dot0 * normalTri); 
                                //if(dot0>0)particles[idx].AddPosition(dot0 * normalTri);

                                //Velocity
                                
                                Vector3 normalVelocityA = Vector3.Dot(normalTri, particles[ed.indexEdgeA].Velocity) * normalTri;
                                Vector3 normalVelocityB = Vector3.Dot(normalTri, particles[ed.indexEdgeB].Velocity) * normalTri;
                                Vector3 normalVelocity1 = Vector3.Dot(normalTri, particles[tri.indexTriA].Velocity) * normalTri;
                                Vector3 normalVelocity2 = Vector3.Dot(normalTri, particles[tri.indexTriB].Velocity) * normalTri;
                                Vector3 normalVelocity3 = Vector3.Dot(normalTri, particles[tri.indexTriC].Velocity) * normalTri;
                                
                                /*Vector3 normalVelocity0 = Vector3.Dot(corr0, particles[idx].Velocity) * corr0;
                                Vector3 normalVelocity1 = Vector3.Dot(corr0, particles[tri.indexTriA].Velocity) * corr0;
                                Vector3 normalVelocity2 = Vector3.Dot(corr0, particles[tri.indexTriB].Velocity) * corr0;
                                Vector3 normalVelocity3 = Vector3.Dot(corr0, particles[tri.indexTriC].Velocity) * corr0;*/
                                
                                Vector3 tangencialVelocityA = particles[ed.indexEdgeA].Velocity - normalVelocityA;
                                Vector3 tangencialVelocityB = particles[ed.indexEdgeB].Velocity - normalVelocityB;
                                Vector3 tangencialVelocity1 = particles[tri.indexTriA].Velocity - normalVelocity1;
                                Vector3 tangencialVelocity2 = particles[tri.indexTriB].Velocity - normalVelocity2;
                                Vector3 tangencialVelocity3 = particles[tri.indexTriC].Velocity - normalVelocity3;

                                Vector3 velocityA = - normalVelocityA;
                                Vector3 velocityB = - normalVelocityB;
                                Vector3 velocity1 = - normalVelocity1;
                                Vector3 velocity2 = - normalVelocity2;
                                Vector3 velocity3 = - normalVelocity3;

                                Vector3 velocityBestA = tangencialVelocityA - frictionConstCloth*normalVelocityA.magnitude*(tangencialVelocityA/tangencialVelocityA.magnitude) - dissipationConstCloth*normalVelocityA;
                                Vector3 velocityBestB = tangencialVelocityB - frictionConstCloth*normalVelocityB.magnitude*(tangencialVelocityB/tangencialVelocityB.magnitude) - dissipationConstCloth*normalVelocityB;
                                Vector3 velocityBest1 = tangencialVelocity1 - frictionConstCloth*normalVelocity1.magnitude*(tangencialVelocity1/tangencialVelocity1.magnitude) - dissipationConstCloth*normalVelocity1;
                                Vector3 velocityBest2 = tangencialVelocity2 - frictionConstCloth*normalVelocity2.magnitude*(tangencialVelocity2/tangencialVelocity2.magnitude) - dissipationConstCloth*normalVelocity2;
                                Vector3 velocityBest3 = tangencialVelocity3 - frictionConstCloth*normalVelocity3.magnitude*(tangencialVelocity3/tangencialVelocity3.magnitude) - dissipationConstCloth*normalVelocity3;
                                
                                if(situation == 3)
                                {
                                    //if(dotA < 0f) particles[ed.indexEdgeA].Position -= dotA * normalTri;
                                    if(dotA < 0f)
                                    {
                                        particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Position - dt * velocityBestA;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Position + dt * velocityBest1 * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Position + dt * velocityBest2 * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Position + dt * velocityBest3 * val3;                                    
                                    }
                                    //if(dotA > 0f) particles[ed.indexEdgeA].Position += dotA * normalTri;
                                    if(dotA > 0f)
                                    {
                                        particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Position + dt * velocityBestA;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Position - dt * velocityBest1 * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Position - dt * velocityBest2 * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Position - dt * velocityBest3 * val3;                                    
                                    }
                                    //particles[ed.indexEdgeA].Position = corrT2;
                                    //particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Prev - normalTri*dotA;
                                }
                                if(situation == 4)
                                {
                                    //if(dotB < 0f) particles[ed.indexEdgeB].Position -= dotB * normalTri;
                                    if(dotB < 0f)
                                    {
                                        particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Position - dt * velocityBestB;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Position + dt * velocityBest1 * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Position + dt * velocityBest2 * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Position + dt * velocityBest3 * val3;
                                    }
                                    //if(dotB > 0f) particles[ed.indexEdgeB].Position += dotB * normalTri;
                                    if(dotB > 0f)
                                    {
                                        particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Position + dt * velocityBestB;
                                        particles[tri.indexTriA].Position = particles[tri.indexTriA].Position - dt * velocityBest1 * val1;
                                        particles[tri.indexTriB].Position = particles[tri.indexTriB].Position - dt * velocityBest2 * val2;
                                        particles[tri.indexTriC].Position = particles[tri.indexTriC].Position - dt * velocityBest3 * val3;
                                    }
                                    //particles[ed.indexEdgeB].Position = corrT2;
                                    //particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Prev - normalTri*dotB;
                                }                               
                                /*if(situation == 3)
                                {
                                    //if(dotA < 0f) particles[ed.indexEdgeA].Position -= dotA * normalTri;
                                    particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Position - dt * velocityBestA;
                                    //if(dotA > 0f) particles[ed.indexEdgeA].Position += dotA * normalTri;
                                    //if(dotB > 0f) particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Position + dt * velocityBestA;
                                    //particles[ed.indexEdgeA].Position = corrT2;
                                    //particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Prev - normalTri*dotA;
                                }
                                if(situation == 4)
                                {
                                    //if(dotB < 0f) particles[ed.indexEdgeB].Position -= dotB * normalTri;
                                    particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Position - dt * velocityBestB;
                                    //if(dotB > 0f) particles[ed.indexEdgeB].Position += dotB * normalTri;
                                    //if(dotB > 0f) particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Position + dt * velocityBestB;
                                    //particles[ed.indexEdgeB].Position = corrT2;
                                    //particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Prev - normalTri*dotB;
                                }*/
                                /*if(situation == 3)
                                {
                                    particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Position - dt * velocityBestA;
                                }
                                if(situation == 4)
                                {
                                    particles[ed.indexEdgeB].Position = particles[ed.indexEdgeB].Position - dt * velocityBestB;
                                }*/

                                //Position
                                //particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Prev + dt * velocityBestA;
                                //particles[ed.indexEdgeA].Position = particles[ed.indexEdgeA].Prev + dt * velocityBestB;
                                //particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev - dt * velocityBest1 * val1;
                                //particles[tri.indexTriB].Position = particles[tri.indexTriB].Prev - dt * velocityBest2 * val2;
                                //particles[tri.indexTriC].Position = particles[tri.indexTriC].Prev - dt * velocityBest3 * val3;


                                /*particles[idx].Position = particles[idx].Prev - normalVelocity0 * 0.02f;
                                particles[tri.indexTriA].Position = particles[tri.indexTriA].Prev + normalVelocity1 * 0.02f; 
                                particles[tri.indexTriB].Position = particles[tri.indexTriA].Prev + normalVelocity2 * 0.02f;
                                particles[tri.indexTriC].Position = particles[tri.indexTriA].Prev + normalVelocity3 * 0.02f;*/

                                //particles[idx].Velocity = -normalVelocity0; //tangencialVelocity0 - 0.5f * normalVelocity0.magnitude * (tangencialVelocity0/tangencialVelocity0.magnitude) - normalVelocity0;//corr;
                                //particles[tri.indexTriA].Velocity = -normalVelocity1;//1.0f * (tangencialVelocity1 - 0.5f * normalVelocity1.magnitude * (tangencialVelocity1/tangencialVelocity1.magnitude) - normalVelocity1);
                                //particles[tri.indexTriB].Velocity = -normalVelocity2;//1.0f * (tangencialVelocity2 - 0.5f * normalVelocity2.magnitude * (tangencialVelocity2/tangencialVelocity2.magnitude) - normalVelocity2);
                                //particles[tri.indexTriC].Velocity = -normalVelocity3;//1.0f * (tangencialVelocity3 - 0.5f * normalVelocity3.magnitude * (tangencialVelocity3/tangencialVelocity3.magnitude) - normalVelocity3);                                

                                //Force
                                /*
                                Vector3 normalForce0 = Vector3.Dot(normalTri, particles[idx].Force) * normalTri;
                                Vector3 normalForce1 = Vector3.Dot(normalTri, particles[tri.indexTriA].Force) * normalTri;
                                Vector3 normalForce2 = Vector3.Dot(normalTri, particles[tri.indexTriB].Force) * normalTri;
                                Vector3 normalForce3 = Vector3.Dot(normalTri, particles[tri.indexTriC].Force) * normalTri;

                                Vector3 tangencialForce0 = particles[idx].Force - normalForce0;
                                Vector3 tangencialForce1 = particles[tri.indexTriA].Force - normalForce1;
                                Vector3 tangencialForce2 = particles[tri.indexTriB].Force - normalForce2;
                                Vector3 tangencialForce3 = particles[tri.indexTriC].Force - normalForce3;
                                */
                            }   
                        }
                    }
                }
            }
        }
    }

    public bool PointTriangleIntersect( 
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
        //Debug.Log(e);
        //if(C < 0) Debug.Log(C);
        //Debug.Log("First: " + p + "   Dist: " + nn+ "   Second: " + e + "  Orig:" + orig + "  Fin:" + fin);
        if(se == 0f) return false;

            Vector3 vvv = e - p;
            float disti = Vector3.Dot(vvv, normalTri);
            Vector3 projectedpoint = e - disti*normalTri;

        //Point triangle collision
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

        

        /*se = C/se;
        if(C < 0f) se *= compressionStiffness;
        else se *=stretchStiffness;*/

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
