using System.Collections;
using System.Collections.Generic;
using UnityEngine;





namespace PositionBasedDynamics
{
  public class PBD
  {
    // epsilon (a value that determine if the changes is too small to change the position)
    public static float eps = 1e-6f;

    public static void UpdatePosition(MeshData md)
    {
      for (int i=0; i < md.particles.Length; i++)
      {
        md.particles[i].velocity = md.particles[i].predictedPos - md.particles[i].pos;
      }
    }

    public static bool ExternalForce(
      float dt,
      Vector3 p, Vector3 v, float w,
      Vector3 force,
      float damping,
      out Vector3 corr)
    {
      corr = Vector3.zero;
      if (w == 0.0f) return false;

      v += dt * w * force * damping;
      corr = dt * v;
      return true;
    }

    public static bool DistanceConstraint(
      Vector3 p0, float w0,
      Vector3 p1, float w1,
      float restLength,
      float stretchStiffness,
      float compressionStiffness,
      out Vector3 corr0,
      out Vector3 corr1)
    {
      corr0 = corr1 = Vector3.zero;
      float wSum = w0 + w1;
      if (wSum == 0.0f) return false;

      Vector3 n = p0 - p1;
      float d = n.magnitude;

      Vector3 corr;
      if (d < restLength)
      {
        corr = compressionStiffness * n * (d - restLength) / wSum;
      } else
      {
        corr = stretchStiffness * n * (d - restLength) / wSum;
      }

      corr0 = -w0 * corr;
      corr1 = w1 * corr;
      return true;
    }

    public static bool DihedralConstraint(
      Vector3 p0, float w0,
      Vector3 p1, float w1,
      Vector3 p2, float w2,
      Vector3 p3, float w3,
      float restAngle,
      float stiffness,
      out Vector3 corr0,
      out Vector3 corr1,
      out Vector3 corr2,
      out Vector3 corr3)
    {
      corr0 = corr1 = corr2 = corr3 = Vector3.zero;
      if (w0 == 0.0f && w1 == 0.0f) return false;

      Vector3 e = p3 - p2;
      float elen = e.magnitude;
      if (elen < eps) return false;

      float invElen = 1 / elen;

      Vector3 n1 = Vector3.Cross((p2 - p0), (p3 - p0)); n1 /= n1.sqrMagnitude;
      Vector3 n2 = Vector3.Cross((p3 - p1), (p2 - p1)); n2 /= n2.sqrMagnitude;

      Vector3 d0 = elen * n1;
      Vector3 d1 = elen * n2;
      Vector3 d2 = Vector3.Dot((p0 - p3), e) * invElen * n1 + Vector3.Dot((p1-p3), e) * invElen * n2;
      Vector3 d3 = Vector3.Dot((p2 - p0), e) * invElen * n1 + Vector3.Dot((p2-p1), e) * invElen * n2;

      Vector3.Normalize(n1);
      Vector3.Normalize(n2);
      float dot = Vector3.Dot(n1, n2);

      if (dot < -1.0f) dot = -1.0f;
      if (dot >  1.0f) dot =  1.0f;
      float phi = Mathf.Acos(dot);	

      // float phi = (-0.6981317 * dot * dot - 0.8726646) * dot + 1.570796;	// fast approximation

      float lambda = 
        w0 * d0.sqrMagnitude +
        w1 * d1.sqrMagnitude +
        w2 * d2.sqrMagnitude +
        w3 * d3.sqrMagnitude;

      if (lambda == 0.0f) return false;	

      // stability
      // 1.5 is the largest magic number I found to be stable in all cases :-)
      //if (stiffness > 0.5 && fabs(phi - b.restAngle) > 1.5)		
      //	stiffness = 0.5;

      lambda = (phi - restAngle) / lambda * stiffness;

      if (Vector3.Dot(Vector3.Cross(n1, n2), e) > 0.0f) lambda = -lambda;	

      corr0 = - w0 * lambda * d0;
      corr1 = - w1 * lambda * d1;
      corr2 = - w2 * lambda * d2;
      corr3 = - w3 * lambda * d3;
      return true;
    }

    public static bool VolumeConstraint(
      Vector3 p0, float w0,
      Vector3 p1, float w1,
      Vector3 p2, float w2,
      Vector3 p3, float w3,
      float restVolume,
      float negVolumeStiffness,
      float posVolumeStiffness,
      out Vector3 corr0,
      out Vector3 corr1,
      out Vector3 corr2,
      out Vector3 corr3)
    {
      corr0 = corr1 = corr2 = corr3 = Vector3.zero;
      float volume = (1 / 6) * Vector3.Dot(Vector3.Cross((p1 - p0), (p2 - p0)), (p3 - p0));

      if (posVolumeStiffness == 0.0f && volume > 0.0f) return false;
      if (negVolumeStiffness == 0.0f && volume < 0.0f) return false;

      Vector3 grad0 = Vector3.Cross((p1 - p2), (p3 - p2));
      Vector3 grad1 = Vector3.Cross((p2 - p0), (p3 - p0));
      Vector3 grad2 = Vector3.Cross((p0 - p1), (p3 - p1));
      Vector3 grad3 = Vector3.Cross((p1 - p0), (p2 - p0));

      float lambda = 
        w0 * grad0.sqrMagnitude +
        w1 * grad1.sqrMagnitude +
        w2 * grad2.sqrMagnitude +
        w3 * grad3.sqrMagnitude;

      if (Mathf.Abs(lambda) < eps) return false;

      if (volume < 0.0f) lambda = negVolumeStiffness * (volume - restVolume) / lambda;
      else lambda = posVolumeStiffness * (volume - restVolume) / lambda;

      corr0 = -lambda * w0 * grad0;
      corr1 = -lambda * w1 * grad1;
      corr2 = -lambda * w2 * grad2;
      corr3 = -lambda * w3 * grad3;

      return true;
    }

    public static bool EdgePointDistanceConstraint(
      Vector3 p, float w,
      Vector3 p0, float w0,
      Vector3 p1, float w1,
      float restDist,
      float compressionStiffness,
      float stretchStiffness,
      out Vector3 corr,
      out Vector3 corr0,
      out Vector3 corr1)
    {
      corr = corr0 = corr1 = Vector3.zero;

      Vector3 d = p1 - p0;
      float t;
      if ((p0 - p1).sqrMagnitude < eps * eps) t = 0.5f;
      else
      {
        float d2 = Vector3.Dot(d, d);
        t = Vector3.Dot(d, (p - p1)) / d2;
        if (t < 0.0f) t = 0.0f;
        else if (t > 1.0f) t = 1.0f;
      }

      // closest point on edge
      Vector3 q = p0 + d*t;
      Vector3 n = p - q;
      float dist = n.magnitude;
      Vector3.Normalize(n);
      float C = dist - restDist;
      float b0 = 1.0f - t;
      float b1 = t;

      Vector3 grad = n;
      Vector3 grad0 = -n * b0;
      Vector3 grad1 = -n * b1;

      float s = w + w0 * b0 * b0 + w1 * b1 * b1;
      if (s == 0.0f) return false;

      s = C / s;
      if (C < 0.0f) s *= compressionStiffness;
      else s *= stretchStiffness;

      if (s == 0.0f) return false;

      corr = -s * w * grad;
      corr0 = -s * w0 * grad0;
      corr1 = -s * w1 * grad1;

      return true;
    }

    public static bool TrianglePointDistanceConstraint(
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
      out Vector3 corr2)
    {
      corr = corr0 = corr1 = corr2 = Vector3.zero;
      // find barycentric coordinates of closest point on triangle

      // for singular case
      float b0 = 1.0f / 3.0f;
      float b1 = b0;
      float b2 = b0;

      float a, b, c, d, e, f;
      float det;

      Vector3 d1 = p1 - p0;
      Vector3 d2 = p2 - p0;
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

      return true;
    }

    public static bool EdgeEdgeDistanceConstraint(
      Vector3 p0, float w0,
      Vector3 p1, float w1,
      Vector3 p2, float w2,
      Vector3 p3, float w3,
      float restDist,
      float compressionStiffness,
      float stretchStiffness,
      out Vector3 corr0,
      out Vector3 corr1,
      out Vector3 corr2,
      out Vector3 corr3)
    {
      corr0 = corr1 = corr2 = corr3 = Vector3.zero;
      Vector3 d0 = p1 - p0;
      Vector3 d1 = p3 - p2;

      float a, b, c, d, e, f;

      a = d0.sqrMagnitude;
      b = -Vector3.Dot(d0, d1);
      c = Vector3.Dot(d0, d1);
      d = -d1.sqrMagnitude;
      e = Vector3.Dot((p2 - p0), d0);
      f = Vector3.Dot((p2 - p0), d1);
      float det = a*d - b*c;
      float s, t;
      if (det != 0.0f)
      {
        det = 1.0f / det;
        s = (e*d - b*f) * det;
        t = (a*f - e*c) * det;
      }
      else
      {
        // d0 and d1 parallel
        float s0 = Vector3.Dot(p0, d0);
        float s1 = Vector3.Dot(p1, d0);
        float t0 = Vector3.Dot(p2, d0);
        float t1 = Vector3.Dot(p3, d0);
        bool flip0 = false;
        bool flip1 = false;

        if (s0 > s1) {f = s0; s0 = s1; s1 = f; flip0 = true;}
        if (t0 > t1) {f = t0; t0 = t1; t1 = f; flip1 = true;}

        if (s0 >= t1)
        {
          s = !flip0 ? 0.0f : 1.0f;
          t = !flip1 ? 1.0f : 0.0f;
        } else if (t0 >= s1)
        {
          s = !flip0 ? 1.0f : 0.0f;
          t = !flip1 ? 0.0f : 1.0f;
        } else
        {
          // overlap
          float mid = (s0 > t0) ? (s0 + t1) * 0.5f : (t0 + s1) * 0.5f;
          s = (s0 == s1) ? 0.5f : (mid - s0) / (s1 - s0);
          t = (t0 == t1) ? 0.5f : (mid - t0) / (t1 - t0);
        }
      }
      if (s < 0.0) s = 0.0f;
      if (s > 1.0) s = 1.0f;
      if (t < 0.0) t = 0.0f;
      if (t > 1.0) t = 1.0f;

      float b0 = 1.0f - s;
      float b1 = s;
      float b2 = 1.0f - t;
      float b3 = t;

      Vector3 q0 = p0 * b0 + p1 * b1;
      Vector3 q1 = p2 * b2 + p3 * b3;
      Vector3 n = q0 - q1;
      float dist = n.magnitude;
      Vector3.Normalize(n);
      float C = dist - restDist;
      Vector3 grad0 = n * b0;
      Vector3 grad1 = n * b1;
      Vector3 grad2 = -n * b2;
      Vector3 grad3 = -n * b3;

      s = w0 * b0*b0 + w1 * b1*b1 + w2 * b2*b2 + w3 * b3*b3;
      if (s == 0.0) return false;

      s = C / s;
      if (C < 0.0) s *= compressionStiffness;
      else s *= stretchStiffness;

      if (s == 0.0) return false;

      corr0 = -s * w0 * grad0;
      corr1 = -s * w1 * grad1;
      corr2 = -s * w2 * grad2;
      corr3 = -s * w3 * grad3;

      return true;
    }

  }
}


















 /*public class TriangleIntersection
 {
     private static void Sort(Vector2 v)
     {
         if (v.x > v.y)
         {
             float c;
             c = v.x;
             v.x = v.y;
             v.y = c;
         }
     }
     
     /// <summary>
     /// This edge to edge test is based on Franlin Antonio's gem: "Faster Line Segment Intersection", in Graphics Gems III, pp. 199-202 
     /// </summary>
     private static bool EdgeEdgeTest(Vector3 v0, Vector3 v1, Vector3 u0, Vector3 u1, int i0, int i1)
     {
         float Ax, Ay, Bx, By, Cx, Cy, e, d, f;
         Ax = v1[i0] - v0[i0];
         Ay = v1[i1] - v0[i1];
 
         Bx = u0[i0] - u1[i0];
         By = u0[i1] - u1[i1];
         Cx = v0[i0] - u0[i0];
         Cy = v0[i1] - u0[i1];
         f = Ay * Bx - Ax * By;
         d = By * Cx - Bx * Cy;
         if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f))
         {
             e = Ax * Cy - Ay * Cx;
             if (f > 0)
             {
                 if (e >= 0 && e <= f) { return true; }
             }
             else
             {
                 if (e <= 0 && e >= f) { return true; }
             }
         }
 
         return false;
     }
 
     private static bool EdgeAgainstTriEdges(Vector3 v0, Vector3 v1, Vector3 u0, Vector3 u1, Vector3 u2, short i0, short i1)
     {
         // test edge u0,u1 against v0,v1
         if (EdgeEdgeTest(v0, v1, u0, u1, i0, i1)) { return true; }
 
         // test edge u1,u2 against v0,v1 
         if (EdgeEdgeTest(v0, v1, u1, u2, i0, i1)) { return true; }
 
         // test edge u2,u1 against v0,v1 
         if (EdgeEdgeTest(v0, v1, u2, u0, i0, i1)) { return true; }
 
         return false;
     }
 
     private static bool PointInTri(Vector3 v0, Vector3 u0, Vector3 u1, Vector3 u2, short i0, short i1)
     {
         float a, b, c, d0, d1, d2;
 
         // is T1 completly inside T2?
         // check if v0 is inside tri(u0,u1,u2)
         a = u1[i1] - u0[i1];
         b = -(u1[i0] - u0[i0]);
         c = -a * u0[i0] - b * u0[i1];
         d0 = a * v0[i0] + b * v0[i1] + c;
 
         a = u2[i1] - u1[i1];
         b = -(u2[i0] - u1[i0]);
         c = -a * u1[i0] - b * u1[i1];
         d1 = a * v0[i0] + b * v0[i1] + c;
 
         a = u0[i1] - u2[i1];
         b = -(u0[i0] - u2[i0]);
         c = -a * u2[i0] - b * u2[i1];
         d2 = a * v0[i0] + b * v0[i1] + c;
 
         if (d0 * d1 > 0.0f)
         {
             if (d0 * d2 > 0.0f) { return true; }
         }
 
         return false;
     }
 
     private static bool TriTriCoplanar(Vector3 N, Vector3 v0, Vector3 v1, Vector3 v2, Vector3 u0, Vector3 u1, Vector3 u2)
     {
         float[] A = new float[3];
         short i0, i1;
 
         // first project onto an axis-aligned plane, that maximizes the area
         // of the triangles, compute indices: i0,i1. 
         A[0] = Mathf.Abs(N[0]);
         A[1] = Mathf.Abs(N[1]);
         A[2] = Mathf.Abs(N[2]);
         if (A[0] > A[1])
         {
             if (A[0] > A[2])
             {
                 // A[0] is greatest
                 i0 = 1;      
                 i1 = 2;
             }
             else
             {
                 // A[2] is greatest
                 i0 = 0;      
                 i1 = 1;
             }
         }
         else  
         {
             if (A[2] > A[1])
             {
                 // A[2] is greatest 
                 i0 = 0;      
                 i1 = 1;
             }
             else
             {
                 // A[1] is greatest 
                 i0 = 0;      
                 i1 = 2;
             }
         }
 
         // test all edges of triangle 1 against the edges of triangle 2 
         if (EdgeAgainstTriEdges(v0, v1, u0, u1, u2, i0, i1)) { return true; }
         if (EdgeAgainstTriEdges(v1, v2, u0, u1, u2, i0, i1)) { return true; }
         if (EdgeAgainstTriEdges(v2, v0, u0, u1, u2, i0, i1)) { return true; }
 
         // finally, test if tri1 is totally contained in tri2 or vice versa 
         if (PointInTri(v0, u0, u1, u2, i0, i1)) { return true; }
         if (PointInTri(u0, v0, v1, v2, i0, i1)) { return true; }
 
         return false;
     }
 
 
 
     private static bool ComputeIntervals(float VV0, float VV1, float VV2,
                               float D0, float D1, float D2, float D0D1, float D0D2,
                               ref float A, ref float B, ref float C, ref float X0, ref float X1)
     {
         if (D0D1 > 0.0f)
         {
             // here we know that D0D2<=0.0 
             // that is D0, D1 are on the same side, D2 on the other or on the plane 
             A = VV2; B = (VV0 - VV2) * D2; C = (VV1 - VV2) * D2; X0 = D2 - D0; X1 = D2 - D1;
         }
         else if (D0D2 > 0.0f)
         {
             // here we know that d0d1<=0.0 
             A = VV1; B = (VV0 - VV1) * D1; C = (VV2 - VV1) * D1; X0 = D1 - D0; X1 = D1 - D2;
         }
         else if (D1 * D2 > 0.0f || D0 != 0.0f)
         {
             // here we know that d0d1<=0.0 or that D0!=0.0 
             A = VV0; B = (VV1 - VV0) * D0; C = (VV2 - VV0) * D0; X0 = D0 - D1; X1 = D0 - D2;
         }
         else if (D1 != 0.0f)
         {
             A = VV1; B = (VV0 - VV1) * D1; C = (VV2 - VV1) * D1; X0 = D1 - D0; X1 = D1 - D2;
         }
         else if (D2 != 0.0f)
         {
             A = VV2; B = (VV0 - VV2) * D2; C = (VV1 - VV2) * D2; X0 = D2 - D0; X1 = D2 - D1;
         }
         else
         {
             return true;
         }
 
         return false;
     }
 
     /// <summary>
     /// Checks if the triangle V(v0, v1, v2) intersects the triangle U(u0, u1, u2).
     /// </summary>
     /// <param name="v0">Vertex 0 of V</param>
     /// <param name="v1">Vertex 1 of V</param>
     /// <param name="v2">Vertex 2 of V</param>
     /// <param name="u0">Vertex 0 of U</param>
     /// <param name="u1">Vertex 1 of U</param>
     /// <param name="u2">Vertex 2 of U</param>
     /// <returns>Returns <c>true</c> if V intersects U, otherwise <c>false</c></returns>
     public bool TriTriIntersect(Vector3 v0, Vector3 v1, Vector3 v2, Vector3 u0, Vector3 u1, Vector3 u2)
     {
         Vector3 e1, e2;
         Vector3 n1, n2;    
         Vector3 dd;
         Vector2 isect1 = Vector2.zero, isect2 = Vector2.zero;
 
         float du0, du1, du2, dv0, dv1, dv2, d1, d2;
         float du0du1, du0du2, dv0dv1, dv0dv2;
         float vp0, vp1, vp2;
         float up0, up1, up2;
         float bb, cc, max;
 
         short index;
 
         // compute plane equation of triangle(v0,v1,v2) 
         e1 = v1 - v0;
         e2 = v2 - v0;
         n1 = Vector3.Cross(e1, e2);
         d1 = -Vector3.Dot(n1, v0);
         // plane equation 1: N1.X+d1=0 */
 /*
         // put u0,u1,u2 into plane equation 1 to compute signed distances to the plane
         du0 = Vector3.Dot(n1, u0) + d1;
         du1 = Vector3.Dot(n1, u1) + d1;
         du2 = Vector3.Dot(n1, u2) + d1;
 
         // coplanarity robustness check 
         if (Mathf.Abs(du0) < Mathf.Epsilon) { du0 = 0.0f; }
         if (Mathf.Abs(du1) < Mathf.Epsilon) { du1 = 0.0f; }
         if (Mathf.Abs(du2) < Mathf.Epsilon) { du2 = 0.0f; }
 
         du0du1 = du0 * du1;
         du0du2 = du0 * du2;
 
         // same sign on all of them + not equal 0 ? 
         if (du0du1 > 0.0f && du0du2 > 0.0f)
         {
             // no intersection occurs
             return false;                    
         }
 
         // compute plane of triangle (u0,u1,u2)
         e1 = u1 - u0;
         e2 = u2 - u0;
         n2 = Vector3.Cross(e1, e2);
         d2 = -Vector3.Dot(n2, u0);
 
         // plane equation 2: N2.X+d2=0 
         // put v0,v1,v2 into plane equation 2
         dv0 = Vector3.Dot(n2, v0) + d2;
         dv1 = Vector3.Dot(n2, v1) + d2;
         dv2 = Vector3.Dot(n2, v2) + d2;
 
         if (Mathf.Abs(dv0) < Mathf.Epsilon) { dv0 = 0.0f; }
         if (Mathf.Abs(dv1) < Mathf.Epsilon) { dv1 = 0.0f; }
         if (Mathf.Abs(dv2) < Mathf.Epsilon) { dv2 = 0.0f; }
 
 
         dv0dv1 = dv0 * dv1;
         dv0dv2 = dv0 * dv2;
 
         // same sign on all of them + not equal 0 ? 
         if (dv0dv1 > 0.0f && dv0dv2 > 0.0f)
         {
             // no intersection occurs
             return false;                    
         }
 
         // compute direction of intersection line 
         dd = Vector3.Cross(n1, n2);
 
         // compute and index to the largest component of D 
         max = (float)Mathf.Abs(dd[0]);
         index = 0;
         bb = (float)Mathf.Abs(dd[1]);
         cc = (float)Mathf.Abs(dd[2]);
         if (bb > max) { max = bb; index = 1; }
         if (cc > max) { max = cc; index = 2; }
 
         // this is the simplified projection onto L
         vp0 = v0[index];
         vp1 = v1[index];
         vp2 = v2[index];
 
         up0 = u0[index];
         up1 = u1[index];
         up2 = u2[index];
 
         // compute interval for triangle 1 
         float a = 0, b = 0, c = 0, x0 = 0, x1 = 0;
         if (ComputeIntervals(vp0, vp1, vp2, dv0, dv1, dv2, dv0dv1, dv0dv2, ref a, ref b, ref c, ref x0, ref x1)) 
         { 
             return TriTriCoplanar(n1, v0, v1, v2, u0, u1, u2); 
         }
 
         // compute interval for triangle 2 
         float d = 0, e = 0, f = 0, y0 = 0, y1 = 0;
         if (ComputeIntervals(up0, up1, up2, du0, du1, du2, du0du1, du0du2, ref d, ref e, ref f, ref y0, ref y1)) 
         { 
             return TriTriCoplanar(n1, v0, v1, v2, u0, u1, u2); 
         }
 
         float xx, yy, xxyy, tmp;
         xx = x0 * x1;
         yy = y0 * y1;
         xxyy = xx * yy;
 
         tmp = a * xxyy;
         isect1[0] = tmp + b * x1 * yy;
         isect1[1] = tmp + c * x0 * yy;
 
         tmp = d * xxyy;
         isect2[0] = tmp + e * xx * y1;
         isect2[1] = tmp + f * xx * y0;
 
         Sort(isect1);
         Sort(isect2);
 
         return !(isect1[1] < isect2[0] || isect2[1] < isect1[0]);
     }
 }*/
