using System;
using System.Collections.Generic;
using System.Linq;

using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace BeingAliveTerrain {
  public class CutAndFill : GH_Component {
    public CutAndFill()
        : base("CutFill", "batCut",
               "Compute the 'Cut and Fill' difference between two given terrains.", "BAT",
               "MeshAnalysis") {}

    public override GH_Exposure Exposure => GH_Exposure.primary;
    public override Guid ComponentGuid => new Guid("db9f6fd9-572d-4304-9fd9-230ee00c29fd");

    protected override void RegisterInputParams(GH_InputParamManager pManager) {
      pManager.AddCurveParameter(
          "Boundary", "B", "One or a set of boundaries of the analyzed area.", GH_ParamAccess.list);
      pManager.AddMeshParameter("Existing", "E",
                                "Existing terrain mesh covering the analyzed area.",
                                GH_ParamAccess.item);
      pManager.AddMeshParameter("Proposed", "P",
                                "Proposed terrain mesh covering the analyzed area.",
                                GH_ParamAccess.item);
      pManager.AddIntegerParameter("GridSize", "G",
                                   "The size of a single grid cell to sample the terrain meshes. " +
                                       "Smaller grid sizes lead to more accurate results, " +
                                       "but also increase computation time.",
                                   GH_ParamAccess.item, 10);
      pManager[3].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager) {
      pManager.AddMeshParameter("AnalysisMesh", "M", "Mesh showing the cut and fill analysis.",
                                GH_ParamAccess.item);
      pManager.AddNumberParameter(
          "CutVolume", "C", "Volume of material that needs to be excavated.", GH_ParamAccess.list);
      pManager.AddNumberParameter("FillVolume", "F", "Volume of material that needs to be added.",
                                  GH_ParamAccess.list);
      pManager.AddNumberParameter(
          "NetVolume", "N",
          "Net volume of material that needs to be added or excavated. " +
              "Positive values indicate a net fill, negative values a net cut.",
          GH_ParamAccess.list);
    }

    protected override void SolveInstance(IGH_DataAccess DA) {
      List<Curve> boundaries = new List<Curve>();
      Mesh existingMesh = null;
      Mesh proposedMesh = null;
      int gridSize = 10;

#region data input& validation
      if (!DA.GetDataList("Boundary", boundaries))
        return;
      if (!DA.GetData("Existing", ref existingMesh))
        return;
      if (!DA.GetData("Proposed", ref proposedMesh))
        return;
      if (!DA.GetData("GridSize", ref gridSize))
        return;

      // Validate inputs
      if (gridSize <= 0) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Grid size must be positive");
        return;
      }

      if (existingMesh == null || proposedMesh == null) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input meshes cannot be null");
        return;
      }
#endregion

      // Create combined bounding box from both meshes
      BoundingBox existingBbox = existingMesh.GetBoundingBox(false);
      BoundingBox proposedBbox = proposedMesh.GetBoundingBox(false);
      BoundingBox combinedBbox = BoundingBox.Union(existingBbox, proposedBbox);

      // Create global grid covering the entire area
      List<Point3d> gridCenters = CreateGridCenters(combinedBbox, gridSize);

      // Sample heights at all grid points
      List<double> existingHeights = new List<double>();
      List<double> proposedHeights = new List<double>();
      List<double> heightDifferences = new List<double>();

      foreach (Point3d center in gridCenters) {
        double existingHeight = GetMeshHeightAtPoint(existingMesh, center);
        double proposedHeight = GetMeshHeightAtPoint(proposedMesh, center);

        existingHeights.Add(existingHeight);
        proposedHeights.Add(proposedHeight);

        // Calculate height difference (0 if either height is invalid)
        double heightDiff = 0;
        if (!double.IsNaN(existingHeight) && !double.IsNaN(proposedHeight)) {
          heightDiff = proposedHeight - existingHeight;
        }
        heightDifferences.Add(heightDiff);
      }

      // Calculate volumes for each boundary
      List<double> cutVolumes = new List<double>();
      List<double> fillVolumes = new List<double>();
      List<double> netVolumes = new List<double>();

      foreach (Curve boundary in boundaries) {
        if (boundary == null || !boundary.IsPlanar()) {
          AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Boundary curve must be planar.");
          cutVolumes.Add(0);
          fillVolumes.Add(0);
          netVolumes.Add(0);
          continue;
        }

        double boundaryCutVolume = 0;
        double boundaryFillVolume = 0;
        double cellArea = gridSize * gridSize;

        for (int i = 0; i < gridCenters.Count; i++) {
          Point3d center = gridCenters[i];
          PointContainment containment = boundary.Contains(center, Plane.WorldXY, 0.01);

          if (containment == PointContainment.Inside ||
              containment == PointContainment.Coincident) {
            double heightDiff = heightDifferences[i];
            double cellVolume = heightDiff * cellArea;

            if (cellVolume > 0) {
              boundaryFillVolume += cellVolume;
            } else {
              boundaryCutVolume += Math.Abs(cellVolume);
            }
          }
        }

        cutVolumes.Add(boundaryCutVolume);
        fillVolumes.Add(boundaryFillVolume);
        netVolumes.Add(boundaryFillVolume - boundaryCutVolume);
      }

      // Create analysis mesh with RedWhiteGreen colormap (white center for zero cut/fill)
      Mesh analysisMesh = CreateAnalysisMesh(gridCenters, existingHeights, heightDifferences,
                                             boundaries, gridSize, combinedBbox, "RedWhiteGreen");

#region data output
      DA.SetData("AnalysisMesh", analysisMesh);
      DA.SetDataList("CutVolume", cutVolumes);
      DA.SetDataList("FillVolume", fillVolumes);
      DA.SetDataList("NetVolume", netVolumes);
#endregion
    }

    private List<Point3d> CreateGridCenters(BoundingBox bbox, double gridSize) {
      List<Point3d> centers = new List<Point3d>();

      int xCount = (int)Math.Ceiling((bbox.Max.X - bbox.Min.X) / gridSize);
      int yCount = (int)Math.Ceiling((bbox.Max.Y - bbox.Min.Y) / gridSize);

      for (int i = 0; i < xCount; i++) {
        for (int j = 0; j < yCount; j++) {
          double x = bbox.Min.X + (i + 0.5) * gridSize;
          double y = bbox.Min.Y + (j + 0.5) * gridSize;
          centers.Add(new Point3d(x, y, 0));
        }
      }

      return centers;
    }

    private double GetMeshHeightAtPoint(Mesh mesh, Point3d point) {
      // Create a ray pointing downward from high above the point
      Point3d rayStart = new Point3d(point.X, point.Y, 10000);
      Vector3d rayDirection = new Vector3d(0, 0, -1);
      Ray3d ray = new Ray3d(rayStart, rayDirection);

      // Find intersection with mesh
      double rayParameter = Rhino.Geometry.Intersect.Intersection.MeshRay(mesh, ray);

      if (rayParameter >= 0) {
        // Return the Z coordinate of the intersection point
        Point3d intersectionPoint = ray.PointAt(rayParameter);
        return intersectionPoint.Z;
      }

      return double.NaN;
    }

    private Mesh CreateAnalysisMesh(List<Point3d> gridCenters, List<double> existingHeights,
                                    List<double> heightDifferences, List<Curve> boundaries,
                                    double gridSize, BoundingBox bbox, string colorMapName) {
      Mesh analysisMesh = new Mesh();

      // Get the specified color map
      IColorMap colorMap = ColorMapHelper.GetColorMap(colorMapName);

      // Calculate grid dimensions
      int xCount = (int)Math.Ceiling((bbox.Max.X - bbox.Min.X) / gridSize);
      int yCount = (int)Math.Ceiling((bbox.Max.Y - bbox.Min.Y) / gridSize);

      // Find min and max height differences for color mapping (excluding zeros)
      List<double> validDifferences = heightDifferences.Where(d => d != 0).ToList();
      double minDiff = validDifferences.Count > 0 ? validDifferences.Min() : -1;
      double maxDiff = validDifferences.Count > 0 ? validDifferences.Max() : 1;

      // Ensure symmetric range around zero for better visualization
      double range = Math.Max(Math.Abs(minDiff), Math.Abs(maxDiff));
      minDiff = -range;
      maxDiff = range;

      // Calculate base Z level (average of bounding box Z values or 0)
      double baseZ = (bbox.Min.Z + bbox.Max.Z) / 2.0;
      if (double.IsNaN(baseZ))
        baseZ = 0;

      // Add vertices and colors
      for (int i = 0; i < gridCenters.Count; i++) {
        Point3d center = gridCenters[i];
        double existingHeight = existingHeights[i];
        double heightDiff = heightDifferences[i];

        // Calculate vertex Z based on height difference
        // Use height difference to create 3D visualization
        double vertexZ = baseZ + heightDiff;

        analysisMesh.Vertices.Add(center.X, center.Y, vertexZ);

        // Check if point is inside any boundary
        bool insideBoundary = false;
        foreach (Curve boundary in boundaries) {
          if (boundary != null && boundary.IsPlanar()) {
            PointContainment containment = boundary.Contains(center, Plane.WorldXY, 0.01);
            if (containment == PointContainment.Inside ||
                containment == PointContainment.Coincident) {
              insideBoundary = true;
              break;
            }
          }
        }

        // Apply color based on height difference and boundary status
        System.Drawing.Color color;
        if (!insideBoundary || heightDiff == 0) {
          color = System.Drawing.Color.White;  // White for flat/outside boundary areas
        } else {
          // Use the specified colormap
          color = ColorMapHelper.MapValueToColor(heightDiff, minDiff, maxDiff, colorMap);
        }
        analysisMesh.VertexColors.Add(color);
      }

      // Create triangular mesh faces
      for (int i = 0; i < xCount - 1; i++) {
        for (int j = 0; j < yCount - 1; j++) {
          int v0 = i * yCount + j;
          int v1 = (i + 1) * yCount + j;
          int v2 = (i + 1) * yCount + (j + 1);
          int v3 = i * yCount + (j + 1);

          // Ensure all vertices are valid
          if (v0 < analysisMesh.Vertices.Count && v1 < analysisMesh.Vertices.Count &&
              v2 < analysisMesh.Vertices.Count && v3 < analysisMesh.Vertices.Count) {
            // Create two triangular faces instead of one quad face
            // First triangle: v0, v1, v2
            analysisMesh.Faces.AddFace(v0, v1, v2);
            // Second triangle: v0, v2, v3
            analysisMesh.Faces.AddFace(v0, v2, v3);
          }
        }
      }

      return analysisMesh;
    }
  }
}  // namespace BeingAliveTerrain
