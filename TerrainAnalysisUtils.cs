using System;
using System.Collections.Generic;
using Rhino.Geometry;

namespace BeingAliveTerrain {
/// <summary>
/// Utility class for terrain analysis operations shared across components
/// </summary>
public
static class TerrainAnalysisUtils {
  /// <summary>
  /// Get the height of the mesh at a specific XY location using ray casting
  /// </summary>
  /// <param name="mesh">Terrain mesh</param>
  /// <param name="point">Point to sample (XY coordinates used, Z
  /// ignored)</param> <returns>Height (Z coordinate) at the point, or NaN if
  /// outside mesh</returns>
 public
  static double GetMeshHeightAtPoint(Mesh mesh, Point3d point) {
    // Create a ray pointing downward from high above the point
    Point3d rayStart = new Point3d(point.X, point.Y, 10000);
    Vector3d rayDirection = new Vector3d(0, 0, -1);
    Ray3d ray = new Ray3d(rayStart, rayDirection);

    // Find intersection with mesh
    double rayParameter =
        Rhino.Geometry.Intersect.Intersection.MeshRay(mesh, ray);

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
 public
  static Point3d ProjectPointOntoTerrain(Mesh mesh, Point3d point) {
    // Try from above first
    Point3d rayStart = new Point3d(point.X, point.Y, 10000);
    Vector3d rayDirection = new Vector3d(0, 0, -1);
    Ray3d ray = new Ray3d(rayStart, rayDirection);

    double rayParameter =
        Rhino.Geometry.Intersect.Intersection.MeshRay(mesh, ray);

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
  /// </summary>
  /// <param name="mesh">Terrain mesh</param>
  /// <param name="point">Point to calculate gradient at</param>
  /// <param name="sampleDistance">Distance to sample in each direction for
  /// gradient calculation</param> <returns>Gradient vector (dz/dx, dz/dy, 0),
  /// or zero vector if calculation fails</returns>
 public
  static Vector3d CalculateGradient(Mesh mesh, Point3d point,
                                    double sampleDistance) {
    // Sample heights in 4 cardinal directions using central difference
    double heightEast =
        GetMeshHeightAtPoint(mesh, point + new Vector3d(sampleDistance, 0, 0));
    double heightWest =
        GetMeshHeightAtPoint(mesh, point + new Vector3d(-sampleDistance, 0, 0));
    double heightNorth =
        GetMeshHeightAtPoint(mesh, point + new Vector3d(0, sampleDistance, 0));
    double heightSouth =
        GetMeshHeightAtPoint(mesh, point + new Vector3d(0, -sampleDistance, 0));

    // Check if any samples are invalid
    if (double.IsNaN(heightEast) || double.IsNaN(heightWest) ||
        double.IsNaN(heightNorth) || double.IsNaN(heightSouth)) {
      return Vector3d.Zero;
    }

    // Calculate gradient using central difference: (f(x+h) - f(x-h)) / (2h)
    double dzdx = (heightEast - heightWest) / (2.0 * sampleDistance);
    double dzdy = (heightNorth - heightSouth) / (2.0 * sampleDistance);

    return new Vector3d(dzdx, dzdy, 0);
  }

  /// <summary>
  /// Calculate the slope angle at a point in degrees
  /// </summary>
  /// <param name="mesh">Terrain mesh</param>
  /// <param name="point">Point to calculate slope at</param>
  /// <param name="sampleDistance">Distance to sample for gradient
  /// calculation</param> <returns>Slope angle in degrees (0 = flat, 90 =
  /// vertical), or NaN if calculation fails</returns>
 public
  static double CalculateSlope(Mesh mesh, Point3d point,
                               double sampleDistance) {
    Vector3d gradient = CalculateGradient(mesh, point, sampleDistance);

    if (gradient == Vector3d.Zero) {
      return double.NaN;
    }

    // Slope = arctan(sqrt(dz/dx² + dz/dy²))
    double gradientMagnitude =
        Math.Sqrt(gradient.X * gradient.X + gradient.Y * gradient.Y);
    double slopeRadians = Math.Atan(gradientMagnitude);
    double slopeDegrees = slopeRadians * (180.0 / Math.PI);

    return slopeDegrees;
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
 public
  static double CalculateAspect(Mesh mesh, Point3d point,
                                double sampleDistance) {
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
 public
  static Vector3d GetSteepestDescentDirection(Mesh mesh, Point3d point,
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
  /// <param name="sampleDistance">Distance to sample for gradient
  /// calculation</param> <returns>Gradient map with grid information</returns>
 public
  static GradientMap PrecomputeGradientMap(Mesh mesh, BoundingBox bbox,
                                           double gridSize,
                                           double sampleDistance) {
    int xCount = (int)Math.Ceiling((bbox.Max.X - bbox.Min.X) / gridSize);
    int yCount = (int)Math.Ceiling((bbox.Max.Y - bbox.Min.Y) / gridSize);

    Vector3d[, ] gradients = new Vector3d[xCount, yCount];
    double[, ] slopes = new double[xCount, yCount];
    double[, ] heights = new double[xCount, yCount];

    for (int i = 0; i < xCount; i++) {
      for (int j = 0; j < yCount; j++) {
        double x = bbox.Min.X + (i + 0.5) * gridSize;
        double y = bbox.Min.Y + (j + 0.5) * gridSize;
        Point3d point = new Point3d(x, y, 0);

        heights[i, j] = GetMeshHeightAtPoint(mesh, point);
        gradients[i, j] = CalculateGradient(mesh, point, sampleDistance);
        slopes[i, j] = CalculateSlope(mesh, point, sampleDistance);
      }
    }

    return new GradientMap{BoundingBox = bbox,    GridSize = gridSize,
                           XCount = xCount,       YCount = yCount,
                           Gradients = gradients, Slopes = slopes,
                           Heights = heights};
  }
}

/// <summary>
/// Pre-computed gradient map for efficient terrain analysis
/// </summary>
public class GradientMap {
 public
  BoundingBox BoundingBox {
    get;
    set;
  }
 public
  double GridSize {
    get;
    set;
  }
 public
  int XCount {
    get;
    set;
  }
 public
  int YCount {
    get;
    set;
  }
 public
  Vector3d[, ] Gradients {
    get;
    set;
  }
 public
  double[, ] Slopes {
    get;
    set;
  }
 public
  double[, ] Heights {
    get;
    set;
  }

  /// <summary>
  /// Get gradient at a specific point using bilinear interpolation
  /// </summary>
 public
  Vector3d GetGradientAt(Point3d point) {
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
 public
  double GetSlopeAt(Point3d point) {
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
