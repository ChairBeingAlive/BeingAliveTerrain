using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace BeingAliveTerrain {
  public class CutAndFill : GH_Component {
    // Store the selected color map name
    private string selectedColorMap = "RedWhiteGreen";

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
      pManager.AddMeshParameter("CutMeshes", "CM",
                                "Mesh faces in cut areas (where material is removed).",
                                GH_ParamAccess.list);
      pManager.AddMeshParameter("FillMeshes", "FM",
                                "Mesh faces in fill areas (where material is added).",
                                GH_ParamAccess.list);
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

      // Create analysis mesh with selected colormap
      Mesh analysisMesh = CreateAnalysisMesh(gridCenters, existingHeights, heightDifferences,
                                             boundaries, gridSize, combinedBbox, selectedColorMap);

      // Extract cut and fill mesh pieces
      var (cutMeshes, fillMeshes) =
          ExtractCutFillMeshes(analysisMesh, heightDifferences, boundaries, gridSize, combinedBbox);

#region data output
      DA.SetData("AnalysisMesh", analysisMesh);
      DA.SetDataList("CutMeshes", cutMeshes);
      DA.SetDataList("FillMeshes", fillMeshes);
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

    /// <summary>
    /// Extract cut and fill mesh pieces by filtering faces based on vertex height differences
    /// </summary>
    private (List<Mesh> cutMeshes, List<Mesh> fillMeshes)
        ExtractCutFillMeshes(Mesh analysisMesh, List<double> heightDifferences,
                             List<Curve> boundaries, double gridSize, BoundingBox bbox) {
      List<Mesh> cutMeshes = new List<Mesh>();
      List<Mesh> fillMeshes = new List<Mesh>();

      // Create separate meshes for cut and fill
      Mesh cutMesh = new Mesh();
      Mesh fillMesh = new Mesh();

      // Iterate through all faces in the analysis mesh
      for (int faceIdx = 0; faceIdx < analysisMesh.Faces.Count; faceIdx++) {
        MeshFace face = analysisMesh.Faces[faceIdx];

        // Get the vertices of this face
        int v0 = face.A;
        int v1 = face.B;
        int v2 = face.C;
        bool isQuad = face.IsQuad;
        int v3 = isQuad ? face.D : -1;

        // Check if any vertex is out of bounds
        if (v0 >= heightDifferences.Count || v1 >= heightDifferences.Count ||
            v2 >= heightDifferences.Count || (isQuad && v3 >= heightDifferences.Count)) {
          continue;
        }

        // Get height differences for the face vertices
        double val0 = heightDifferences[v0];
        double val1 = heightDifferences[v1];
        double val2 = heightDifferences[v2];
        double val3 = isQuad ? heightDifferences[v3] : 0;

        // Calculate average height difference for the face
        double avgHeight = isQuad ? (val0 + val1 + val2 + val3) / 4.0 : (val0 + val1 + val2) / 3.0;

        // Check if face center is inside any boundary
        Point3d faceCenter =
            isQuad ? new Point3d((analysisMesh.Vertices[v0].X + analysisMesh.Vertices[v1].X +
                                  analysisMesh.Vertices[v2].X + analysisMesh.Vertices[v3].X) /
                                     4.0,
                                 (analysisMesh.Vertices[v0].Y + analysisMesh.Vertices[v1].Y +
                                  analysisMesh.Vertices[v2].Y + analysisMesh.Vertices[v3].Y) /
                                     4.0,
                                 (analysisMesh.Vertices[v0].Z + analysisMesh.Vertices[v1].Z +
                                  analysisMesh.Vertices[v2].Z + analysisMesh.Vertices[v3].Z) /
                                     4.0)
                   : new Point3d((analysisMesh.Vertices[v0].X + analysisMesh.Vertices[v1].X +
                                  analysisMesh.Vertices[v2].X) /
                                     3.0,
                                 (analysisMesh.Vertices[v0].Y + analysisMesh.Vertices[v1].Y +
                                  analysisMesh.Vertices[v2].Y) /
                                     3.0,
                                 (analysisMesh.Vertices[v0].Z + analysisMesh.Vertices[v1].Z +
                                  analysisMesh.Vertices[v2].Z) /
                                     3.0);

        bool insideBoundary = false;
        foreach (Curve boundary in boundaries) {
          if (boundary != null && boundary.IsPlanar()) {
            if (boundary.Contains(new Point3d(faceCenter.X, faceCenter.Y, 0), Plane.WorldXY,
                                  0.01) != PointContainment.Outside) {
              insideBoundary = true;
              break;
            }
          }
        }

        if (!insideBoundary)
          continue;

        // Classify face as cut or fill based on average height difference
        double threshold = 0.001;
        if (avgHeight < -threshold) {
          // Cut face - add to cut mesh
          int newV0 = cutMesh.Vertices.Add(analysisMesh.Vertices[v0]);
          int newV1 = cutMesh.Vertices.Add(analysisMesh.Vertices[v1]);
          int newV2 = cutMesh.Vertices.Add(analysisMesh.Vertices[v2]);

          cutMesh.VertexColors.Add(analysisMesh.VertexColors[v0]);
          cutMesh.VertexColors.Add(analysisMesh.VertexColors[v1]);
          cutMesh.VertexColors.Add(analysisMesh.VertexColors[v2]);

          if (isQuad && v3 >= 0) {
            int newV3 = cutMesh.Vertices.Add(analysisMesh.Vertices[v3]);
            cutMesh.VertexColors.Add(analysisMesh.VertexColors[v3]);
            cutMesh.Faces.AddFace(newV0, newV1, newV2, newV3);
          } else {
            cutMesh.Faces.AddFace(newV0, newV1, newV2);
          }
        } else if (avgHeight > threshold) {
          // Fill face - add to fill mesh
          int newV0 = fillMesh.Vertices.Add(analysisMesh.Vertices[v0]);
          int newV1 = fillMesh.Vertices.Add(analysisMesh.Vertices[v1]);
          int newV2 = fillMesh.Vertices.Add(analysisMesh.Vertices[v2]);

          fillMesh.VertexColors.Add(analysisMesh.VertexColors[v0]);
          fillMesh.VertexColors.Add(analysisMesh.VertexColors[v1]);
          fillMesh.VertexColors.Add(analysisMesh.VertexColors[v2]);

          if (isQuad && v3 >= 0) {
            int newV3 = fillMesh.Vertices.Add(analysisMesh.Vertices[v3]);
            fillMesh.VertexColors.Add(analysisMesh.VertexColors[v3]);
            fillMesh.Faces.AddFace(newV0, newV1, newV2, newV3);
          } else {
            fillMesh.Faces.AddFace(newV0, newV1, newV2);
          }
        }
      }

      // Add meshes to lists if they have faces
      if (cutMesh.Faces.Count > 0) {
        cutMesh.Normals.ComputeNormals();
        cutMesh.Compact();
        cutMeshes.Add(cutMesh);
      }

      if (fillMesh.Faces.Count > 0) {
        fillMesh.Normals.ComputeNormals();
        fillMesh.Compact();
        fillMeshes.Add(fillMesh);
      }

      return (cutMeshes, fillMeshes);
    }

    // Add right-click menu for color map selection using Eto
    // (cross-platform)
    public override void AppendAdditionalMenuItems(ToolStripDropDown menu) {
      base.AppendAdditionalMenuItems(menu);

      // Add a separator
      Menu_AppendSeparator(menu);

      // Add color map selection submenu
      var colorMapMenu = Menu_AppendItem(menu, "Color Map");

      foreach (var colorMap in ColorMapHelper.AvailableColorMaps) {
        Menu_AppendItem(colorMapMenu.DropDown, colorMap.Name, OnColorMapSelected, true,
                        colorMap.Name == selectedColorMap)
            .Tag = colorMap.Name;
      }
    }

    private void OnColorMapSelected(object sender, EventArgs e) {
      if (sender is ToolStripMenuItem item && item.Tag is string colorMapName) {
        selectedColorMap = colorMapName;
        ExpireSolution(true);  // Recompute with new color map
      }
    }

    // Override read/write to persist color map selection
    public override bool Write(GH_IO.Serialization.GH_IWriter writer) {
      writer.SetString("ColorMap", selectedColorMap);
      return base.Write(writer);
    }

    public override bool Read(GH_IO.Serialization.GH_IReader reader) {
      if (reader.ItemExists("ColorMap")) {
        selectedColorMap = reader.GetString("ColorMap");
      }
      return base.Read(reader);
    }
  }
}  // namespace BeingAliveTerrain
