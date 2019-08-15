using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using Unity.Mathematics;
using static Unity.Mathematics.math;
using Nebukam.Geom;


namespace Nebukam.Splines
{
    
    /// <summary>
    /// Managed path allow tight point management as well as custom Point types, and as such may be used for tool authoring
    /// while still retaining path-related functionalities.
    /// </summary>
    public abstract class VertexPath<V> : VertexGroup<V>
        where V : Vertex, IVertex, new()
    {
        
        public bool loop { get; set; } = true;

        #region Catmull-Rom Spline

        #region Interpolation

        /// <summary>
        ///Get a point on a Catmull-Rom spline.
        ///The percentage is in range 0 to 1, which starts at the second control point and ends at the second last control point. 
        ///The array cPoints should contain all control points. The minimum amount of control points should be 4. 
        ///Based on : https://forum.unity.com/threads/waypoints-and-constant-variable-speed-problems.32954/#post-213942
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        public float3 CRInterp(float t)
        {
            int numSections = m_vertices.Count - 3;
            int currPt = Mathf.Min(Mathf.FloorToInt(t * (float)numSections), numSections - 1);
            float u = t * (float)numSections - (float)currPt, uu = u*u, uuu = uu * u;

            float3 a = m_vertices[currPt].pos, 
                b = m_vertices[currPt + 1].pos,
                c = m_vertices[currPt + 2].pos,
                d = m_vertices[currPt + 3].pos;

            float ax = a.x, ay = a.y, az = a.z, amx = -ax, amy = -ay, amz = -az, bx = b.x, by = b.y, bz = b.z, cx = c.x, cy = c.y, cz = c.z, dx = d.x, dy = d.y, dz = d.z;

            return float3(
                .5f * ((amx + 3f * bx - 3f * cx + dx) * uuu + (2f * ax - 5f * bx + 4f * cx - dx) * uu + (amx + cx) * u + 2f * bx),
                .5f * ((amy + 3f * by - 3f * cy + dy) * uuu + (2f * ay - 5f * by + 4f * cy - dy) * uu + (amy + cy) * u + 2f * by),
                .5f * ((amz + 3f * bz - 3f * cz + dz) * uuu + (2f * az - 5f * bz + 4f * cz - dz) * uu + (amz + cz) * u + 2f * bz)
                );
        }

        /// <summary>
        ///Get a point on a Catmull-Rom spline, clamped between two nodes ( index range 1:n-1 ).
        ///The percentage is in range 0 to 1, which starts at the second control point and ends at the second last control point.
        /// </summary>
        /// <param name="from">index range 1:n-1</param>
        /// <param name="t"></param>
        /// <returns></returns>
        public float3 CRInterpClamped(int from, float t)
        {

            float u = t, uu = u * u, uuu = uu * u;

            float3 a = m_vertices[from-1].pos,
                b = m_vertices[from].pos,
                c = m_vertices[from + 1].pos,
                d = m_vertices[from + 2].pos;

            float ax = a.x, ay = a.y, az = a.z, amx = -ax, amy = -ay, amz = -az, bx = b.x, by = b.y, bz = b.z, cx = c.x, cy = c.y, cz = c.z, dx = d.x, dy = d.y, dz = d.z;

            return float3(
                .5f * ((amx + 3f * bx - 3f * cx + dx) * uuu + (2f * ax - 5f * bx + 4f * cx - dx) * uu + (amx + cx) * u + 2f * bx),
                .5f * ((amy + 3f * by - 3f * cy + dy) * uuu + (2f * ay - 5f * by + 4f * cy - dy) * uu + (amy + cy) * u + 2f * by),
                .5f * ((amz + 3f * bz - 3f * cz + dz) * uuu + (2f * az - 5f * bz + 4f * cz - dz) * uu + (amz + cz) * u + 2f * bz)
                );
        }

        #endregion

        #region Velocity

        /// <summary>
        /// Catmul-Rom Spline velocity solver
        /// Require at least 4 points
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        public float3 CRVelocity(float t)
        {
            int numSections = m_vertices.Count - 3;
            int currPt = Mathf.Min(Mathf.FloorToInt(t * (float)numSections), numSections - 1);
            float u = t * (float)numSections - (float)currPt, uu = u * u;

            float3 a = m_vertices[currPt].pos,
                b = m_vertices[currPt + 1].pos,
                c = m_vertices[currPt + 2].pos,
                d = m_vertices[currPt + 3].pos;

            float ax = a.x, ay = a.y, az = a.z, bx = b.x, by = b.y, bz = b.z, cx = c.x, cy = c.y, cz = c.z, dx = d.x, dy = d.y, dz = d.z;

            return float3(
                1.5f * (-ax + 3f * bx - 3f * cx + dx) * uu + (2f * ax - 5f * bx + 4f * cx - d.x) * u + .5f * cx - .5f * ax,
                1.5f * (-ay + 3f * by - 3f * cy + dy) * uu + (2f * ay - 5f * by + 4f * cy - d.y) * u + .5f * cy - .5f * ay,
                1.5f * (-az + 3f * bz - 3f * cz + dz) * uu + (2f * az - 5f * bz + 4f * cz - d.z) * u + .5f * cz - .5f * az
                );
            
        }

        /// <summary>
        /// Catmul-Rom Spline velocity solver, clamped between two nodes ( index range 1:n-1 ).
        /// Require at least 4 points
        /// </summary>
        /// <param name="from">index range 1:n-1</param>
        /// <param name="t"></param>
        /// <returns></returns>
        public float3 CRVelocityClamped(int from, float t)
        {

            float u = t, uu = u * u;

            float3 a = m_vertices[from - 1].pos,
                b = m_vertices[from].pos,
                c = m_vertices[from + 1].pos,
                d = m_vertices[from + 2].pos;

            float ax = a.x, ay = a.y, az = a.z, bx = b.x, by = b.y, bz = b.z, cx = c.x, cy = c.y, cz = c.z, dx = d.x, dy = d.y, dz = d.z;

            return float3(
                1.5f * (-ax + 3f * bx - 3f * cx + dx) * uu + (2f * ax - 5f * bx + 4f * cx - d.x) * u + .5f * cx - .5f * ax,
                1.5f * (-ay + 3f * by - 3f * cy + dy) * uu + (2f * ay - 5f * by + 4f * cy - d.y) * u + .5f * cy - .5f * ay,
                1.5f * (-az + 3f * bz - 3f * cz + dz) * uu + (2f * az - 5f * bz + 4f * cz - d.z) * u + .5f * cz - .5f * az
                );

        }

        #endregion

        #region Closest point on Catmull-Rom Spline

        #endregion

        #endregion

    }

}
