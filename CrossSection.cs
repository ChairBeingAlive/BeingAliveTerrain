using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace BeingAliveTerrain {
public class CrossSection : GH_Component {
  public CrossSection()
      : base("CrossSection",
             "batXSec",
             "Compute cross-section profiles of a mesh terrain along a given curve.",
             "BAT",
             "Analysis") {}

  public override GH_Exposure Exposure => GH_Exposure.primary;
  public override Guid ComponentGuid => new Guid("6c512b4b-2503-4ba0-b49e-896614bbacf4");

  protected override System.Drawing.Bitmap Icon => Properties.Resources.crossSection;

  protected override void RegisterInputParams(GH_InputParamManager pManager) {
    pManager.AddMeshParameter("Mesh", "M", "Existing terrain mesh to analyze.", GH_ParamAccess.item);
    pManager.AddCurveParameter("Curve", "C", "Curve along which to compute cross-sections.", GH_ParamAccess.item);
    pManager.AddIntegerParameter("Count",
                                 "N",
                                 "Number of cross-section sampling points along the curve (evenly distributed).",
                                 GH_ParamAccess.item,
                                 10);
    pManager.AddNumberParameter(
        "Width", "W", "Width of each cross-section profile line (perpendicular to curve).", GH_ParamAccess.item, 10.0);

    pManager[2].Optional = true;
    pManager[3].Optional = true;
  }

  protected override void RegisterOutputParams(GH_OutputParamManager pManager) {
    pManager.AddCurveParameter("Sections", "S", "Cross-section curves from the terrain cuts.", GH_ParamAccess.list);
    pManager.AddPlaneParameter("Planes", "P", "Corresponding planes of each cross-section cut.", GH_ParamAccess.list);
  }

  protected override void SolveInstance(IGH_DataAccess DA) {
    Mesh terrainMesh = null;
    Curve guideCurve = null;
    int sampleCount = 10;
    double sectionWidth = 10.0;

#region Data Input& Validation
    if (!DA.GetData("Mesh", ref terrainMesh))
      return;
    if (!DA.GetData("Curve", ref guideCurve))
      return;
    DA.GetData("Count", ref sampleCount);
    DA.GetData("Width", ref sectionWidth);

    // Validate inputs
    if (terrainMesh == null) {
      AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input mesh cannot be null.");
      return;
    }

    if (guideCurve == null) {
      AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input curve cannot be null.");
      return;
    }

    if (sampleCount < 2) {
      AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Sample count should be at least 2. Using 2.");
      sampleCount = 2;
    }

    if (sectionWidth <= 0) {
      AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Section width must be positive. Using 10.0.");
      sectionWidth = 10.0;
    }
#endregion

    List<Curve> sectionCurves = new List<Curve>();
    List<Plane> sectionPlanes = new List<Plane>();

    // Get the bounding box of the mesh to determine ray casting height
    BoundingBox meshBbox = terrainMesh.GetBoundingBox(false);
    double rayStartHeight = meshBbox.Max.Z + 1000;

    // Sample points along the curve at evenly distributed parameters
    double curveLength = guideCurve.GetLength();
    double halfWidth = sectionWidth / 2.0;

    for (int i = 0; i < sampleCount; i++) {
      // Calculate parameter along curve (normalized from 0 to 1)
      double t = (double)i / (sampleCount - 1);

      // Get the point on the curve at this parameter
      double curveParam;
      guideCurve.LengthParameter(t * curveLength, out curveParam);
      Point3d pointOnCurve = guideCurve.PointAt(curveParam);

      // Get the tangent vector at this point
      Vector3d tangent = guideCurve.TangentAt(curveParam);
      tangent.Unitize();

      // Calculate the perpendicular direction (cross section direction)
      // Cross product of tangent with Z-axis gives perpendicular in XY plane
      Vector3d perpendicular = Vector3d.CrossProduct(tangent, Vector3d.ZAxis);
      perpendicular.Unitize();

      // Sample points along the cross-section width to build the profile
      List<Point3d> profilePoints = new List<Point3d>();
      int widthSamples = Math.Max(20, (int)(sectionWidth / 0.5));  // Sample every 0.5 units minimum

      for (int j = 0; j <= widthSamples; j++) {
        // Calculate position along the section width (-halfWidth to +halfWidth)
        double widthParam = (double)j / widthSamples;
        double offset = -halfWidth + widthParam * sectionWidth;

        // Calculate the 3D position for ray casting
        Point3d samplePoint = pointOnCurve + perpendicular * offset;

        // Cast ray downward to find mesh intersection
        Point3d rayStart = new Point3d(samplePoint.X, samplePoint.Y, rayStartHeight);
        Ray3d ray = new Ray3d(rayStart, -Vector3d.ZAxis);

        double rayParam = Rhino.Geometry.Intersect.Intersection.MeshRay(terrainMesh, ray);

        if (rayParam >= 0) {
          // Found intersection - add point at terrain height
          Point3d intersectionPoint = ray.PointAt(rayParam);
          profilePoints.Add(intersectionPoint);
        }
      }

      // Create a polyline curve from the profile points if we have enough points
      if (profilePoints.Count >= 2) {
        Polyline profilePolyline = new Polyline(profilePoints);
        Curve profileCurve = profilePolyline.ToNurbsCurve();
        sectionCurves.Add(profileCurve);

        // Create the section plane with origin at the midpoint of the profile curve (on terrain)
        // This makes it easier to reorient sections elsewhere using the plane as reference
        Point3d profileMidpoint = profileCurve.PointAtNormalizedLength(0.5);
        Plane sectionPlane = new Plane(profileMidpoint, perpendicular, Vector3d.ZAxis);
        sectionPlanes.Add(sectionPlane);
      } else {
        // Add a null or warning if no valid profile could be created
        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                          $"Could not create profile at sample {i + 1} - insufficient terrain intersection points.");
      }
    }

#region Data Output
    DA.SetDataList("Sections", sectionCurves);
    DA.SetDataList("Planes", sectionPlanes);
#endregion
  }
}
}
