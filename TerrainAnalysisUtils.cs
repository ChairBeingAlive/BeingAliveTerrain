using System;
using System.Collections.Generic;
using Rhino.Geometry;

namespace BeingAliveTerrain {
  /// <summary>
  /// Utility class for terrain analysis operations shared across components
  /// </summary>
  public static class TerrainAnalysisUtils {
    /// <summary>
    /// Get the height of the mesh at a specific XY location using ray casting
    /// </summary>
    /// <param name="mesh">Terrain mesh</param>
    /// <param name="point">Point to sample (XY coordinates used, Z
    /// ignored)</param> <returns>Height (Z coordinate) at the point, or NaN if
    /// outside mesh</returns>
    public static double GetMeshHeightAtPoint(Mesh mesh, Point3d point) {
      // Create a ray pointing downward from high above the point
      Point3d rayStart = new Point3d(point.X, point.Y, 10000);
      Vector3d rayDirection = new Vector3d(0, 0, -1);
      Ray3d ray = new Ray3d(rayStart, rayDirection);

      // Find intersection with mesh
      double rayParameter = Rhino.Geometry.Intersect.Intersection.MeshRay(mesh, ray);

      if (rayParameter >= 0) {
        Point3d intersectionPoint = ray.PointAt(rayParameter);
        return intersectionPoint.Z;
      }

      return double.NaN;
    }

    /// <summary>
    /// Project a point onto the terrain mesh using ray casting
    /// </summary>
    /// <param name="mesh">Terrain mesh</param>
    /// <param name="point">Point to project</param>
    /// <returns>Projected point on mesh surface, or point with NaN Z if
    /// projection fails</returns>
    public static Point3d ProjectPointOntoTerrain(Mesh mesh, Point3d point) {
      // Try from above first
      Point3d rayStart = new Point3d(point.X, point.Y, 10000);
      Vector3d rayDirection = new Vector3d(0, 0, -1);
      Ray3d ray = new Ray3d(rayStart, rayDirection);

      double rayParameter = Rhino.Geometry.Intersect.Intersection.MeshRay(mesh, ray);

      if (rayParameter >= 0) {
        return ray.PointAt(rayParameter);
      }

      // If no intersection found, try from below
      rayStart = new Point3d(point.X, point.Y, -10000);
      rayDirection = new Vector3d(0, 0, 1);
      ray = new Ray3d(rayStart, rayDirection);
      rayParameter = Rhino.Geometry.Intersect.Intersection.MeshRay(mesh, ray);

      if (rayParameter >= 0) {
        return ray.PointAt(rayParameter);
      }

      // Return invalid point if projection failed
      return new Point3d(point.X, point.Y, double.NaN);
    }

    /// <summary>
    /// Calculate the terrain gradient at a point using central difference method
    /// Falls back to one-sided differences at mesh edges where samples are unavailable
    /// </summary>
    /// <param name="mesh">Terrain mesh</param>
    /// <param name="point">Point to calculate gradient at</param>
    /// <param name="sampleDistance">Distance to sample in each direction for gradient
    /// calculation</param>
    /// <returns>Gradient vector (dz/dx, dz/dy, 0), or zero vector if
    /// calculation fails</returns>
    public static Vector3d CalculateGradient(Mesh mesh, Point3d point, double sampleDistance) {
      return CalculateGradientWithStatus(mesh, point, sampleDistance, out _);
    }

    /// <summary>
    /// Calculate the terrain gradient at a point with adaptive method selection
    /// Uses mesh normals near edges and finite differences in the interior
    /// </summary>
    /// <param name="mesh">Terrain mesh</param>
    /// <param name="point">Point to calculate gradient at</param>
    /// <param name="sampleDistance">Distance to sample in each direction</param>
    /// <param name="calculationFailed">Output: true if gradient calculation completely failed,
    /// false if successful or genuinely flat</param> <returns>Gradient vector (dz/dx, dz/dy,
    /// 0)</returns>
    private static Vector3d CalculateGradientWithStatus(Mesh mesh, Point3d point,
                                                        double sampleDistance,
                                                        out bool calculationFailed) {
      calculationFailed = false;

      // Sample heights in 4 cardinal directions
      double heightCenter = GetMeshHeightAtPoint(mesh, point);
      double heightEast = GetMeshHeightAtPoint(mesh, point + new Vector3d(sampleDistance, 0, 0));
      double heightWest = GetMeshHeightAtPoint(mesh, point + new Vector3d(-sampleDistance, 0, 0));
      double heightNorth = GetMeshHeightAtPoint(mesh, point + new Vector3d(0, sampleDistance, 0));
      double heightSouth = GetMeshHeightAtPoint(mesh, point + new Vector3d(0, -sampleDistance, 0));

      // Count how many samples are valid
      int validSamples = 0;
      if (!double.IsNaN(heightCenter))
        validSamples++;
      if (!double.IsNaN(heightEast))
        validSamples++;
      if (!double.IsNaN(heightWest))
        validSamples++;
      if (!double.IsNaN(heightNorth))
        validSamples++;
      if (!double.IsNaN(heightSouth))
        validSamples++;

      // If we have very few valid samples, this point is close to mesh edge
      // Use mesh normal method as primary approach
      if (validSamples <= 2) {
        Vector3d normalBasedGradient = GetGradientFromMeshNormal(mesh, point);
        if (normalBasedGradient != Vector3d.Zero) {
          return normalBasedGradient;
        }
        // If mesh normal method also fails, mark as failed
        calculationFailed = true;
        return Vector3d.Zero;
      }

      // Calculate gradient components using adaptive difference method
      double dzdx = 0;
      double dzdy = 0;
      bool xGradientValid = false;
      bool yGradientValid = false;

      // X-direction gradient (East-West)
      if (!double.IsNaN(heightEast) && !double.IsNaN(heightWest)) {
        // Central difference (best accuracy)
        dzdx = (heightEast - heightWest) / (2.0 * sampleDistance);
        xGradientValid = true;
      } else if (!double.IsNaN(heightEast) && !double.IsNaN(heightCenter)) {
        // Forward difference (at west edge)
        dzdx = (heightEast - heightCenter) / sampleDistance;
        xGradientValid = true;
      } else if (!double.IsNaN(heightWest) && !double.IsNaN(heightCenter)) {
        // Backward difference (at east edge)
        dzdx = (heightCenter - heightWest) / sampleDistance;
        xGradientValid = true;
      }

      // Y-direction gradient (North-South)
      if (!double.IsNaN(heightNorth) && !double.IsNaN(heightSouth)) {
        // Central difference (best accuracy)
        dzdy = (heightNorth - heightSouth) / (2.0 * sampleDistance);
        yGradientValid = true;
      } else if (!double.IsNaN(heightNorth) && !double.IsNaN(heightCenter)) {
        // Forward difference (at south edge)
        dzdy = (heightNorth - heightCenter) / sampleDistance;
        yGradientValid = true;
      } else if (!double.IsNaN(heightSouth) && !double.IsNaN(heightCenter)) {
        // Backward difference (at north edge)
        dzdy = (heightCenter - heightSouth) / sampleDistance;
        yGradientValid = true;
      }

      // If we have partial gradient info, combine it with mesh normal
      if (!xGradientValid || !yGradientValid) {
        Vector3d normalBasedGradient = GetGradientFromMeshNormal(mesh, point);

        if (normalBasedGradient != Vector3d.Zero) {
          // Use finite difference for the valid direction, mesh normal for the invalid direction
          if (!xGradientValid && yGradientValid) {
            dzdx = normalBasedGradient.X;
            xGradientValid = true;
          }
          if (!yGradientValid && xGradientValid) {
            dzdy = normalBasedGradient.Y;
            yGradientValid = true;
          }
          // If both are invalid, use the mesh normal entirely
          if (!xGradientValid && !yGradientValid) {
            return normalBasedGradient;
          }
        }
      }

      // If both gradients failed to calculate, mark as failed
      if (!xGradientValid && !yGradientValid) {
        calculationFailed = true;
      }

      return new Vector3d(dzdx, dzdy, 0);
    }

    /// <summary>
    /// Calculate the slope angle at a point in degrees
    /// </summary>
    /// <param name="mesh">Terrain mesh</param>
    /// <param name="point">Point to calculate slope at</param>
    /// <param name="sampleDistance">Distance to sample for gradient calculation</param>
    /// <returns>Slope angle in degrees (0 = flat, 90 = vertical), or NaN if calculation
    /// fails</returns>
    public static double CalculateSlope(Mesh mesh, Point3d point, double sampleDistance) {
      Vector3d gradient = CalculateGradient(mesh, point, sampleDistance);
      return CalculateSlopeFromGradient(gradient);
    }

    /// <summary>
    /// Calculate the aspect (direction of slope) at a point
    /// </summary>
    /// <param name="mesh">Terrain mesh</param>
    /// <param name="point">Point to calculate aspect at</param>
    /// <param name="sampleDistance">Distance to sample for gradient
    /// calculation</param> <returns>Aspect angle in degrees (0 = North, 90 =
    /// East, 180 = South, 270 = West), or NaN if flat or calculation
    /// fails</returns>
    public static double CalculateAspect(Mesh mesh, Point3d point, double sampleDistance) {
      Vector3d gradient = CalculateGradient(mesh, point, sampleDistance);

      if (gradient == Vector3d.Zero) {
        return double.NaN;
      }

      // Aspect = atan2(-dy/dz, dx/dz) converted to compass bearing
      // Note: negative dy because we want North = 0°
      double aspectRadians = Math.Atan2(-gradient.Y, gradient.X);
      double aspectDegrees = aspectRadians * (180.0 / Math.PI);

      // Convert to compass bearing (0 = North, clockwise)
      aspectDegrees = 90.0 - aspectDegrees;
      if (aspectDegrees < 0) {
        aspectDegrees += 360.0;
      }

      return aspectDegrees;
    }

    /// <summary>
    /// Get the steepest descent direction at a point
    /// </summary>
    /// <param name="mesh">Terrain mesh</param>
    /// <param name="point">Point to calculate direction at</param>
    /// <param name="sampleDistance">Distance to sample for gradient
    /// calculation</param> <returns>Normalized direction vector pointing downhill
    /// (XY plane), or zero if flat/invalid</returns>
    public static Vector3d GetSteepestDescentDirection(Mesh mesh, Point3d point,
                                                       double sampleDistance) {
      Vector3d gradient = CalculateGradient(mesh, point, sampleDistance);

      if (gradient == Vector3d.Zero) {
        return Vector3d.Zero;
      }

      // Steepest descent is in the negative gradient direction (downhill)
      Vector3d descentDirection = new Vector3d(-gradient.X, -gradient.Y, 0);

      // Normalize to unit vector
      descentDirection.Unitize();

      return descentDirection;
    }

    /// <summary>
    /// Pre-compute a gradient map for the entire terrain on a regular grid
    /// Useful for components that need to query gradients frequently
    /// </summary>
    /// <param name="mesh">Terrain mesh</param>
    /// <param name="bbox">Bounding box of area to compute</param>
    /// <param name="gridSize">Grid cell size</param>
    /// <param name="sampleDistance">Distance to sample for gradient calculation</param>
    /// <returns>Gradient map with grid information</returns>
    public static GradientMap PrecomputeGradientMap(Mesh mesh, BoundingBox bbox, double gridSize,
                                                    double sampleDistance) {
      int xCount = (int)Math.Ceiling((bbox.Max.X - bbox.Min.X) / gridSize);
      int yCount = (int)Math.Ceiling((bbox.Max.Y - bbox.Min.Y) / gridSize);

      Vector3d[,] gradients = new Vector3d[xCount, yCount];
      double[,] slopes = new double[xCount, yCount];
      double[,] heights = new double[xCount, yCount];

      int meshNormalUsedCount = 0;
      int meshNormalFallbackCount = 0;
      int flatTerrainCount = 0;

      for (int i = 0; i < xCount; i++) {
        for (int j = 0; j < yCount; j++) {
          double x = bbox.Min.X + (i + 0.5) * gridSize;
          double y = bbox.Min.Y + (j + 0.5) * gridSize;
          Point3d point = new Point3d(x, y, 0);

          heights[i, j] = GetMeshHeightAtPoint(mesh, point);

          // Skip if point is outside mesh
          if (double.IsNaN(heights[i, j])) {
            gradients[i, j] = Vector3d.Zero;
            slopes[i, j] = double.NaN;
            continue;
          }

          // Update point with actual height
          point = new Point3d(x, y, heights[i, j]);

          // Calculate gradient with status information
          bool calculationFailed;
          Vector3d gradient =
              CalculateGradientWithStatus(mesh, point, sampleDistance, out calculationFailed);

          // Check if gradient is from mesh normal by sampling again with small distance
          if (gradient != Vector3d.Zero && !calculationFailed) {
            // Count valid samples to determine if mesh normal was likely used
            double testDistance = sampleDistance;
            int validCount = 0;
            if (!double.IsNaN(GetMeshHeightAtPoint(mesh, point)))
              validCount++;
            if (!double.IsNaN(GetMeshHeightAtPoint(mesh, point + new Vector3d(testDistance, 0, 0))))
              validCount++;
            if (!double.IsNaN(
                    GetMeshHeightAtPoint(mesh, point + new Vector3d(-testDistance, 0, 0))))
              validCount++;
            if (!double.IsNaN(GetMeshHeightAtPoint(mesh, point + new Vector3d(0, testDistance, 0))))
              validCount++;
            if (!double.IsNaN(
                    GetMeshHeightAtPoint(mesh, point + new Vector3d(0, -testDistance, 0))))
              validCount++;

            if (validCount <= 2) {
              meshNormalUsedCount++;
            }
          }

          gradients[i, j] = gradient;

          // If gradient calculation completely failed, it's a fallback case
          if (calculationFailed) {
            meshNormalFallbackCount++;
          } else if (gradient == Vector3d.Zero) {
            // Gradient is zero but calculation succeeded - genuinely flat terrain
            flatTerrainCount++;
          }

          slopes[i, j] = CalculateSlopeFromGradient(gradients[i, j]);
        }
      }

      var gradientMap = new GradientMap { BoundingBox = bbox,
                                          GridSize = gridSize,
                                          XCount = xCount,
                                          YCount = yCount,
                                          Gradients = gradients,
                                          Slopes = slopes,
                                          Heights = heights,
                                          MeshNormalUsedCount = meshNormalUsedCount,
                                          MeshNormalFallbackCount = meshNormalFallbackCount,
                                          FlatTerrainCount = flatTerrainCount };

      return gradientMap;
    }

    /// <summary>
    /// Get gradient estimate from mesh normal at a point (primary method for edge cases)
    /// </summary>
    /// <param name="mesh">Terrain mesh</param>
    /// <param name="point">Point on mesh surface</param>
    /// <returns>Gradient vector estimated from mesh normal, or zero if failed</returns>
    private static Vector3d GetGradientFromMeshNormal(Mesh mesh, Point3d point) {
      // Find closest point on mesh
      MeshPoint meshPoint = mesh.ClosestMeshPoint(point, 0.0);
      if (meshPoint == null) {
        return Vector3d.Zero;
      }

      // Get the face index from the mesh point
      int faceIndex = meshPoint.FaceIndex;
      if (faceIndex < 0 || faceIndex >= mesh.Faces.Count) {
        return Vector3d.Zero;
      }

      // Get face normal
      Vector3d normal = mesh.FaceNormals[faceIndex];

      // If normal is vertical (flat terrain), return zero gradient
      if (Math.Abs(normal.Z) < 0.001) {
        return Vector3d.Zero;
      }

      // Convert normal to gradient
      // Normal = (-dz/dx, -dz/dy, 1) normalized
      // So gradient = (dz/dx, dz/dy) = (-normal.X/normal.Z, -normal.Y/normal.Z)
      double dzdx = -normal.X / normal.Z;
      double dzdy = -normal.Y / normal.Z;

      return new Vector3d(dzdx, dzdy, 0);
    }

    /// <summary>
    /// Calculate slope angle from a gradient vector
    /// </summary>
    /// <param name="gradient">Gradient vector (dz/dx, dz/dy, 0)</param>
    /// <returns>Slope angle in degrees</returns>
    private static double CalculateSlopeFromGradient(Vector3d gradient) {
      if (gradient == Vector3d.Zero) {
        return 0.0;  // Flat terrain has 0° slope
      }

      // Slope = arctan(sqrt(dz/dx² + dz/dy²))
      double gradientMagnitude = Math.Sqrt(gradient.X * gradient.X + gradient.Y * gradient.Y);
      double slopeRadians = Math.Atan(gradientMagnitude);
      double slopeDegrees = slopeRadians * (180.0 / Math.PI);

      return slopeDegrees;
    }
  }

  /// <summary>
  /// Pre-computed gradient map for efficient terrain analysis
  /// </summary>
  public class GradientMap {
    public BoundingBox BoundingBox { get; set; }
    public double GridSize { get; set; }
    public int XCount { get; set; }
    public int YCount { get; set; }
    public Vector3d[,] Gradients { get; set; }
    public double[,] Slopes { get; set; }
    public double[,] Heights { get; set; }

    /// <summary>
    /// Number of points where mesh normal was used as primary method (near edges with few valid
    /// samples)
    /// </summary>
    public int MeshNormalUsedCount { get; set; }

    /// <summary>
    /// Number of points where mesh normal was used as fallback after gradient calculation failed
    /// </summary>
    public int MeshNormalFallbackCount { get; set; }

    /// <summary>
    /// Number of points with genuinely flat terrain (zero gradient)
    /// </summary>
    public int FlatTerrainCount { get; set; }

    /// <summary>
    /// Get gradient at a specific point using bilinear interpolation
    /// </summary>
    public Vector3d GetGradientAt(Point3d point) {
      // Convert world coordinates to grid indices
      double x = (point.X - BoundingBox.Min.X) / GridSize - 0.5;
      double y = (point.Y - BoundingBox.Min.Y) / GridSize - 0.5;

      int i0 = (int)Math.Floor(x);
      int j0 = (int)Math.Floor(y);
      int i1 = i0 + 1;
      int j1 = j0 + 1;

      // Check bounds
      if (i0 < 0 || i1 >= XCount || j0 < 0 || j1 >= YCount) {
        return Vector3d.Zero;
      }

      // Bilinear interpolation weights
      double wx = x - i0;
      double wy = y - j0;

      Vector3d g00 = Gradients[i0, j0];
      Vector3d g10 = Gradients[i1, j0];
      Vector3d g01 = Gradients[i0, j1];
      Vector3d g11 = Gradients[i1, j1];

      Vector3d gx0 = g00 * (1 - wx) + g10 * wx;
      Vector3d gx1 = g01 * (1 - wx) + g11 * wx;
      Vector3d result = gx0 * (1 - wy) + gx1 * wy;

      return result;
    }

    /// <summary>
    /// Get slope at a specific point using bilinear interpolation
    /// </summary>
    public double GetSlopeAt(Point3d point) {
      double x = (point.X - BoundingBox.Min.X) / GridSize - 0.5;
      double y = (point.Y - BoundingBox.Min.Y) / GridSize - 0.5;

      int i0 = (int)Math.Floor(x);
      int j0 = (int)Math.Floor(y);
      int i1 = i0 + 1;
      int j1 = j0 + 1;

      if (i0 < 0 || i1 >= XCount || j0 < 0 || j1 >= YCount) {
        return double.NaN;
      }

      double wx = x - i0;
      double wy = y - j0;

      double s00 = Slopes[i0, j0];
      double s10 = Slopes[i1, j0];
      double s01 = Slopes[i0, j1];
      double s11 = Slopes[i1, j1];

      double sx0 = s00 * (1 - wx) + s10 * wx;
      double sx1 = s01 * (1 - wx) + s11 * wx;
      double result = sx0 * (1 - wy) + sx1 * wy;

      return result;
    }
  }
}  // namespace BeingAliveTerrain
