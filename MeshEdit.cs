using System;
using System.Collections.Generic;
using System.Linq;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace BeingAliveTerrain {
  // Helper class for smooth interpolation functions
  public static class SmoothInterpolation {
    /// <summary>
    /// Smooth step function (3t² - 2t³) - provides smooth transition with zero derivatives at edges
    /// </summary>
    /// <param name="t">Input value between 0 and 1</param>
    /// <returns>Smoothly interpolated value between 0 and 1</returns>
    public static double SmoothStep(double t) {
      t = Math.Max(0, Math.Min(1, t));  // Clamp to [0,1]
      return t * t * (3.0 - 2.0 * t);
    }

    /// <summary>
    /// Smoother step function (6t⁵ - 15t⁴ + 10t³) - even smoother with zero first and second
    /// derivatives at edges
    /// </summary>
    /// <param name="t">Input value between 0 and 1</param>
    /// <returns>Very smoothly interpolated value between 0 and 1</returns>
    public static double SmootherStep(double t) {
      t = Math.Max(0, Math.Min(1, t));  // Clamp to [0,1]
      return t * t * t * (t * (t * 6.0 - 15.0) + 10.0);
    }

    /// <summary>
    /// Cosine interpolation - uses cosine function for smooth transition
    /// /// </summary>
    /// <param name="t">Input value between 0 and 1</param>
    /// <returns>Cosine-interpolated value between 0 and 1</returns>
    public static double CosineInterpolation(double t) {
      t = Math.Max(0, Math.Min(1, t));  // Clamp to [0,1]
      return (1.0 - Math.Cos(t * Math.PI)) * 0.5;
    }
  }

  public class PointEdit : GH_Component {
    public PointEdit()
        : base("TerrainPointEdit", "batPtTerrain", "Using a point to pull/push the mesh terrain.",
               "BAT", "Edit") {}

    public override GH_Exposure Exposure => GH_Exposure.primary;
    public override Guid ComponentGuid => new Guid("6e3e71c3-92f1-4158-8b45-0cdb4efb52da");
    protected override System.Drawing.Bitmap Icon => Properties.Resources.pointEdit;

    protected override void RegisterInputParams(GH_InputParamManager pManager) {
      pManager.AddMeshParameter("Mesh", "M", "Mesh to edit", GH_ParamAccess.item);
      pManager.AddPointParameter("Point", "P", "Point(s) to use for mesh editing",
                                 GH_ParamAccess.list);  // Changed to list access
      pManager.AddNumberParameter(
          "Strength", "S",
          "Editing strength (-100 to 100). Positive values push up, negative values pull down.",
          GH_ParamAccess.item, 0);
      pManager.AddNumberParameter(
          "BlurDistance", "B",
          "Distance for gradual transition between edited and original mesh (world units)",
          GH_ParamAccess.item, 10);

      pManager[2].Optional = true;
      pManager[3].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager) {
      pManager.AddMeshParameter("EditedMesh", "EM", "The edited mesh", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA) {
      Mesh inputMesh = null;
      List<Point3d> editPoints = new List<Point3d>();
      double strength = 0;
      double blurDistance = 10;

      // Get input data
      if (!DA.GetData("Mesh", ref inputMesh))
        return;
      if (!DA.GetDataList("Point", editPoints))  // Changed to GetDataList
        return;
      if (!DA.GetData("Strength", ref strength))
        return;
      if (!DA.GetData("BlurDistance", ref blurDistance))
        return;

      // Validate inputs
      if (inputMesh == null) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input mesh cannot be null");
        return;
      }

      if (editPoints == null || editPoints.Count == 0) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "At least one edit point is required");
        return;
      }

      if (strength < -100 || strength > 100) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                          "Strength should be between -100 and 100");
        strength = Math.Max(-100, Math.Min(100, strength));
      }

      if (blurDistance < 0) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Blur distance should be positive");
        blurDistance = Math.Max(0, blurDistance);
      }

      // Create a copy of the input mesh
      Mesh editedMesh = inputMesh.DuplicateMesh();

      // Apply point-based mesh editing with multiple points
      EditMeshWithPoints(editedMesh, editPoints, strength, blurDistance);

      // Output the edited mesh
      DA.SetData("EditedMesh", editedMesh);
    }

    private void EditMeshWithPoints(Mesh mesh, List<Point3d> editPoints, double strength,
                                    double blurDistance) {
      if (strength == 0 || editPoints.Count == 0)
        return;

      // Work directly with the mesh vertices to maintain order and structure
      for (int i = 0; i < mesh.Vertices.Count; i++) {
        Point3d vertex = mesh.Vertices[i];
        double maxAbsoluteDisplacement = 0;
        double finalDisplacement = 0;

        // Check influence from each edit point
        foreach (Point3d editPoint in editPoints) {
          // Calculate 2D distance (XY plane) from vertex to edit point
          double distance =
              Math.Sqrt(Math.Pow(vertex.X - editPoint.X, 2) + Math.Pow(vertex.Y - editPoint.Y, 2));

          // Check if vertex is within the blur radius
          if (distance <= blurDistance) {
            // Calculate normalized distance (0 = at center, 1 = at edge)
            double normalizedDistance = distance / blurDistance;

            // Apply smooth interpolation for influence (1 = full influence, 0 = no influence)
            double influence = 1.0 - SmoothInterpolation.SmootherStep(normalizedDistance);

            // Calculate target displacement based on strength
            double displacement = 0;
            if (strength > 0) {
              // Positive strength: move toward editPoint.Z
              displacement = (editPoint.Z - vertex.Z) * (strength / 100.0) * influence;
            } else {
              // Negative strength: move away from editPoint.Z
              displacement = -(editPoint.Z - vertex.Z) * (Math.Abs(strength) / 100.0) * influence;
            }

            // Keep the displacement with the largest absolute value
            if (Math.Abs(displacement) > Math.Abs(maxAbsoluteDisplacement)) {
              maxAbsoluteDisplacement = displacement;
              finalDisplacement = displacement;
            }
          }
        }

        // Apply the final displacement if any influence was found
        if (finalDisplacement != 0) {
          double targetZ = vertex.Z + finalDisplacement;
          // Update the vertex in-place, keeping X and Y the same, only changing Z
          mesh.Vertices[i] = new Point3f((float)vertex.X, (float)vertex.Y, (float)targetZ);
        }
      }

      // Recompute normals after vertex modification
      mesh.Normals.ComputeNormals();
    }
  }

  public class CurveEdit : GH_Component {
    public CurveEdit()
        : base("TerrainCurveEdit", "batCrvTerrain", "Using a curve to pull/push the mesh terrain.",
               "BAT", "Edit") {}

    public override GH_Exposure Exposure => GH_Exposure.primary;
    public override Guid ComponentGuid => new Guid("7e4e81d3-93f2-4259-9b46-1cdb5efb62db");
    protected override System.Drawing.Bitmap Icon => Properties.Resources.curveEdit;

    protected override void RegisterInputParams(GH_InputParamManager pManager) {
      pManager.AddMeshParameter("Mesh", "M", "Mesh to edit", GH_ParamAccess.item);
      pManager.AddCurveParameter("Curve", "C", "Curve(s) to use for mesh editing",
                                 GH_ParamAccess.list);  // Changed to list access
      pManager.AddNumberParameter(
          "Strength", "S",
          "Editing strength (-100 to 100). Positive values push up, negative values pull down.",
          GH_ParamAccess.item, 0);
      pManager.AddNumberParameter(
          "BlurDistance", "B",
          "Distance for gradual transition between edited and original mesh (world units)",
          GH_ParamAccess.item, 10);

      pManager[2].Optional = true;
      pManager[3].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager) {
      pManager.AddMeshParameter("EditedMesh", "EM", "The edited mesh", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA) {
      Mesh inputMesh = null;
      List<Curve> editCurves = new List<Curve>();
      double strength = 0;
      double blurDistance = 10;

      // Get input data
      if (!DA.GetData("Mesh", ref inputMesh))
        return;
      if (!DA.GetDataList("Curve", editCurves))  // Changed to GetDataList
        return;
      if (!DA.GetData("Strength", ref strength))
        return;
      if (!DA.GetData("BlurDistance", ref blurDistance))
        return;

      // Validate inputs
      if (inputMesh == null) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input mesh cannot be null");
        return;
      }

      if (editCurves == null || editCurves.Count == 0) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "At least one edit curve is required");
        return;
      }

      // Check for null curves
      for (int i = 0; i < editCurves.Count; i++) {
        if (editCurves[i] == null) {
          AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Input curve {i} cannot be null");
          return;
        }
      }

      if (strength < -100 || strength > 100) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                          "Strength should be between -100 and 100");
        strength = Math.Max(-100, Math.Min(100, strength));
      }

      if (blurDistance < 0) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Blur distance should be positive");
        blurDistance = Math.Max(0, blurDistance);
      }

      // Create a copy of the input mesh
      Mesh editedMesh = inputMesh.DuplicateMesh();

      // Apply curve-based mesh editing with multiple curves
      EditMeshWithCurves(editedMesh, editCurves, strength, blurDistance);

      // Output the edited mesh
      DA.SetData("EditedMesh", editedMesh);
    }

    private void EditMeshWithCurves(Mesh mesh, List<Curve> editCurves, double strength,
                                    double blurDistance) {
      if (strength == 0 || editCurves.Count == 0)
        return;

      // Work directly with the mesh vertices to maintain order and structure
      for (int i = 0; i < mesh.Vertices.Count; i++) {
        Point3d vertex = mesh.Vertices[i];
        Point3d vertexXY = new Point3d(vertex.X, vertex.Y, 0);
        double maxAbsoluteDisplacement = 0;
        double finalDisplacement = 0;

        // Check influence from each edit curve
        foreach (Curve editCurve in editCurves) {
          // Find closest point on curve to vertex (projected to XY plane)
          double t;
          editCurve.ClosestPoint(vertexXY, out t);
          Point3d closestPoint = editCurve.PointAt(t);

          // Calculate 2D distance (XY plane) from vertex to closest point on curve
          double distance = Math.Sqrt(Math.Pow(vertex.X - closestPoint.X, 2) +
                                      Math.Pow(vertex.Y - closestPoint.Y, 2));

          // Check if vertex is within the blur distance
          if (distance <= blurDistance) {
            // Calculate normalized distance (0 = at curve, 1 = at edge of influence)
            double normalizedDistance = distance / blurDistance;

            // Apply smooth interpolation for influence (1 = full influence, 0 = no influence)
            double influence = 1.0 - SmoothInterpolation.SmootherStep(normalizedDistance);

            // Calculate target displacement based on strength (similar to point editing)
            // Use the Z coordinate of the closest point on the curve as target
            double displacement = 0;
            if (strength > 0) {
              // Positive strength: move toward closestPoint.Z
              displacement = (closestPoint.Z - vertex.Z) * (strength / 100.0) * influence;
            } else {
              // Negative strength: move away from closestPoint.Z
              displacement =
                  -(closestPoint.Z - vertex.Z) * (Math.Abs(strength) / 100.0) * influence;
            }

            // Keep the displacement with the largest absolute value
            if (Math.Abs(displacement) > Math.Abs(maxAbsoluteDisplacement)) {
              maxAbsoluteDisplacement = displacement;
              finalDisplacement = displacement;
            }
          }
        }

        // Apply the final displacement if any influence was found
        if (finalDisplacement != 0) {
          double targetZ = vertex.Z + finalDisplacement;
          // Update the vertex in-place, keeping X and Y the same, only changing Z
          mesh.Vertices[i] = new Point3f((float)vertex.X, (float)vertex.Y, (float)targetZ);
        }
      }

      // Recompute normals after vertex modification
      mesh.Normals.ComputeNormals();
    }
  }

  public class AreaEdit : GH_Component {
    public AreaEdit()
        : base("TerrainAreaEdit", "batAreaTerrain",
               "Using a closed curve to pull/push the mesh terrain within the area.", "BAT",
               "Edit") {}

    public override GH_Exposure Exposure => GH_Exposure.primary;
    public override Guid ComponentGuid => new Guid("8e5e91e3-94f3-4360-9c47-2cdb6efb73dc");
    protected override System.Drawing.Bitmap Icon => Properties.Resources.areaEdit;

    protected override void RegisterInputParams(GH_InputParamManager pManager) {
      pManager.AddMeshParameter("Mesh", "M", "Mesh to edit", GH_ParamAccess.item);
      pManager.AddCurveParameter("Area", "A", "Closed curve(s) defining the area(s) to edit",
                                 GH_ParamAccess.list);
      pManager.AddNumberParameter(
          "Strength", "S",
          "Editing strength (-100 to 100). Positive values push up, negative values pull down.",
          GH_ParamAccess.item, 0);
      pManager.AddNumberParameter(
          "BlurDistance", "B",
          "Distance for gradual transition between edited and original mesh (world units)",
          GH_ParamAccess.item, 10);
      pManager.AddBooleanParameter(
          "Flat", "F",
          "If true, averages all vertices inside the area to create a flat surface. If false, moves vertices toward border height.",
          GH_ParamAccess.item, false);

      pManager[2].Optional = true;
      pManager[3].Optional = true;
      pManager[4].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager) {
      pManager.AddMeshParameter("EditedMesh", "EM", "The edited mesh", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA) {
      Mesh inputMesh = null;
      List<Curve> areaCurves = new List<Curve>();
      double strength = 0;
      double blurDistance = 10;
      bool flat = false;

      // Get input data
      if (!DA.GetData("Mesh", ref inputMesh))
        return;
      if (!DA.GetDataList("Area", areaCurves))
        return;
      if (!DA.GetData("Strength", ref strength))
        return;
      if (!DA.GetData("BlurDistance", ref blurDistance))
        return;
      if (!DA.GetData("Flat", ref flat))
        return;

      // Validate inputs
      if (inputMesh == null) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input mesh cannot be null");
        return;
      }

      if (areaCurves == null || areaCurves.Count == 0) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "At least one area curve is required");
        return;
      }

      // Check for null curves and closed status
      for (int i = 0; i < areaCurves.Count; i++) {
        if (areaCurves[i] == null) {
          AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Input area curve {i} cannot be null");
          return;
        }
        if (!areaCurves[i].IsClosed) {
          AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Area curve {i} must be closed");
          return;
        }
        if (!areaCurves[i].IsPlanar()) {
          AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                            $"Area curve {i} should be planar for best results");
        }
      }

      if (strength < -100 || strength > 100) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                          "Strength should be between -100 and 100");
        strength = Math.Max(-100, Math.Min(100, strength));
      }

      if (blurDistance < 0) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Blur distance should be positive");
        blurDistance = Math.Max(0, blurDistance);
      }

      // Create a copy of the input mesh
      Mesh editedMesh = inputMesh.DuplicateMesh();

      // Apply area-based mesh editing with flattening behavior
      EditMeshWithAreas(editedMesh, areaCurves, strength, blurDistance, flat);

      // Output the edited mesh
      DA.SetData("EditedMesh", editedMesh);
    }

    private void EditMeshWithAreas(Mesh mesh, List<Curve> areaCurves, double strength,
                                   double blurDistance, bool flat) {
      if (strength == 0 || areaCurves.Count == 0)
        return;

      // If flat mode is enabled, we need to calculate average heights first
      Dictionary<int, double> averageHeights = new Dictionary<int, double>();

      if (flat) {
        // Calculate average height for each area curve
        for (int curveIndex = 0; curveIndex < areaCurves.Count; curveIndex++) {
          Curve areaCurve = areaCurves[curveIndex];
          List<double> insideHeights = new List<double>();

          // Collect all heights of vertices inside or on the curve
          for (int i = 0; i < mesh.Vertices.Count; i++) {
            Point3d vertex = mesh.Vertices[i];
            Point3d vertexXY = new Point3d(vertex.X, vertex.Y, 0);
            PointContainment containment = areaCurve.Contains(vertexXY, Plane.WorldXY, 0.01);

            if (containment == PointContainment.Inside ||
                containment == PointContainment.Coincident) {
              insideHeights.Add(vertex.Z);
            }
          }

          // Calculate average height
          if (insideHeights.Count > 0) {
            averageHeights[curveIndex] = insideHeights.Average();
          }
        }
      }

      // Work directly with the mesh vertices to maintain order and structure
      for (int i = 0; i < mesh.Vertices.Count; i++) {
        Point3d vertex = mesh.Vertices[i];
        Point3d vertexXY = new Point3d(vertex.X, vertex.Y, 0);
        double maxAbsoluteDisplacement = 0;
        double finalDisplacement = 0;

        // Check influence from each area curve
        for (int curveIndex = 0; curveIndex < areaCurves.Count; curveIndex++) {
          Curve areaCurve = areaCurves[curveIndex];

          // Check if point is inside the closed curve
          PointContainment containment = areaCurve.Contains(vertexXY, Plane.WorldXY, 0.01);

          double influence = 0;
          double targetHeight = 0;

          if (containment == PointContainment.Inside) {
            // Point is inside - move directly to target height (full influence, no blur)
            influence = 1.0;

            // Use border height
            double t;
            areaCurve.ClosestPoint(vertexXY, out t);
            Point3d closestBoundaryPoint = areaCurve.PointAt(t);
            targetHeight = closestBoundaryPoint.Z;
          } else if (containment == PointContainment.Coincident) {
            // Point is on the boundary - also full influence
            influence = 1.0;

            // Use border height
            double t;
            areaCurve.ClosestPoint(vertexXY, out t);
            Point3d closestBoundaryPoint = areaCurve.PointAt(t);
            targetHeight = closestBoundaryPoint.Z;
          } else {
            // Point is outside - apply blur effect only here
            double t;
            areaCurve.ClosestPoint(vertexXY, out t);
            Point3d closestBoundaryPoint = areaCurve.PointAt(t);
            double distanceToBoundary = Math.Sqrt(Math.Pow(vertex.X - closestBoundaryPoint.X, 2) +
                                                  Math.Pow(vertex.Y - closestBoundaryPoint.Y, 2));

            if (distanceToBoundary <= blurDistance) {
              // Calculate normalized distance (0 = at boundary, 1 = at edge of influence)
              double normalizedDistance = distanceToBoundary / blurDistance;
              // Apply smooth interpolation for fade outside
              influence = 1.0 - SmoothInterpolation.SmootherStep(normalizedDistance);

              // Use border height
              targetHeight = closestBoundaryPoint.Z;
            }
          }

          // Apply Z displacement based on strength and influence
          if (influence > 0) {
            double displacement = 0;

            if (flat && averageHeights.ContainsKey(curveIndex) &&
                (containment == PointContainment.Inside ||
                 containment == PointContainment.Coincident)) {
              // Flat mode for inside/coincident vertices:
              // Calculate adjusted average height based on strength moving toward border height
              double averageHeight = averageHeights[curveIndex];
              double borderHeight = targetHeight;
              double adjustedAverageHeight;

              if (strength > 0) {
                // Positive strength: move average height toward border height
                adjustedAverageHeight =
                    averageHeight + (borderHeight - averageHeight) * (strength / 100.0);
              } else {
                // Negative strength: move average height away from border height
                adjustedAverageHeight =
                    averageHeight - (borderHeight - averageHeight) * (Math.Abs(strength) / 100.0);
              }

              // All inside vertices move to this adjusted average height
              displacement = adjustedAverageHeight - vertex.Z;
            } else if (flat && averageHeights.ContainsKey(curveIndex)) {
              // Flat mode for outside vertices: transition toward adjusted average height
              double averageHeight = averageHeights[curveIndex];
              double borderHeight = targetHeight;
              double adjustedAverageHeight;

              if (strength > 0) {
                adjustedAverageHeight =
                    averageHeight + (borderHeight - averageHeight) * (strength / 100.0);
              } else {
                adjustedAverageHeight =
                    averageHeight - (borderHeight - averageHeight) * (Math.Abs(strength) / 100.0);
              }

              // Outside vertices transition toward adjusted average height with influence falloff
              displacement = (adjustedAverageHeight - vertex.Z) * influence;
            } else {
              // Normal mode: move toward target height with strength
              if (strength > 0) {
                displacement = (targetHeight - vertex.Z) * (strength / 100.0) * influence;
              } else {
                displacement =
                    -(targetHeight - vertex.Z) * (Math.Abs(strength) / 100.0) * influence;
              }
            }

            // Keep the displacement with the largest absolute value
            if (Math.Abs(displacement) > Math.Abs(maxAbsoluteDisplacement)) {
              maxAbsoluteDisplacement = displacement;
              finalDisplacement = displacement;
            }
          }
        }

        // Apply the final displacement if any influence was found
        if (finalDisplacement != 0) {
          double targetZ = vertex.Z + finalDisplacement;
          // Update the vertex in-place, keeping X and Y the same, only changing Z
          mesh.Vertices[i] = new Point3f((float)vertex.X, (float)vertex.Y, (float)targetZ);
        }
      }

      // Recompute normals after vertex modification
      mesh.Normals.ComputeNormals();
    }
  }
}
