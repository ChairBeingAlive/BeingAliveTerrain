using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace BeingAliveTerrain {
public class SlopeAnalysis : GH_Component {
  // Store the selected color map name
  private string selectedColorMap = "Viridis";

  public SlopeAnalysis()
      : base("SlopeAnalysis",
             "batSlope",
             "Analyze terrain slope and visualize with color-coded mesh.",
             "BAT",
             "Analysis") {}

  public override GH_Exposure Exposure => GH_Exposure.primary;
  public override Guid ComponentGuid => new Guid("2f4d8bc6-420c-4bf8-a094-e18a1e74981b");
  protected override System.Drawing.Bitmap Icon => Properties.Resources.slopeAnalysis;

  protected override void RegisterInputParams(GH_InputParamManager pManager) {
    pManager.AddMeshParameter("Mesh", "M", "Terrain mesh to analyze.", GH_ParamAccess.item);
    pManager.AddNumberParameter(
        "GridSize",
        "G",
        "Grid size for sampling terrain. Smaller values give more detail but  slower computation.",
        GH_ParamAccess.item,
        10.0);
    pManager.AddIntervalParameter(
        "SlopeBounds",
        "B",
        "Slope range for color mapping (degrees). Slopes outside this range  are clamped to min/max colors. " +
            "Use this to focus on specific slope ranges.",
        GH_ParamAccess.item,
        new Interval(0, 90));

    pManager[1].Optional = true;
    pManager[2].Optional = true;
  }

  protected override void RegisterOutputParams(GH_OutputParamManager pManager) {
    pManager.AddMeshParameter("AnalysisMesh",
                              "M",
                              "Color-coded mesh showing slope analysis. Colors  mapped to SlopeBounds range.",
                              GH_ParamAccess.item);
    pManager.AddPointParameter("SamplePoints", "P", "Grid points where slope was sampled.", GH_ParamAccess.list);
    pManager.AddNumberParameter("SlopeValues", "S", "Slope angle at each sample point (degrees).", GH_ParamAccess.list);
    pManager.AddVectorParameter(
        "SlopeVectors", "V", "Direction vectors showing downhill direction at each point.", GH_ParamAccess.list);
    pManager.AddBooleanParameter("InBounds",
                                 "IB",
                                 "True if slope is within SlopeBounds, False if outside. " +
                                     "Use this to filter/select points by slope range.",
                                 GH_ParamAccess.list);
  }

  protected override void SolveInstance(IGH_DataAccess DA) {
    Mesh terrain = null;
    double gridSize = 10.0;
    Interval slopeBounds = new Interval(0, 90);

    // Get input data
    if (!DA.GetData("Mesh", ref terrain))
      return;
    if (!DA.GetData("GridSize", ref gridSize))
      return;
    if (!DA.GetData("SlopeBounds", ref slopeBounds))
      return;

    // Validate inputs
    if (terrain == null) {
      AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input terrain mesh cannot be null");
      return;
    }

    if (gridSize <= 0) {
      AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "GridSize must be positive");
      gridSize = Math.Max(0.1, gridSize);
    }

    // Validate slope bounds
    if (!slopeBounds.IsValid || slopeBounds.T0 < 0 || slopeBounds.T1 > 90) {
      AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                        "SlopeBounds should be a valid interval between 0 and 90 degrees");
      slopeBounds = new Interval(Math.Max(0, slopeBounds.T0), Math.Min(90, slopeBounds.T1));
    }

    if (slopeBounds.T0 >= slopeBounds.T1) {
      AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "SlopeBounds minimum must be less than maximum");
      slopeBounds = new Interval(0, 90);
    }

    // Get terrain bounding box
    BoundingBox bbox = terrain.GetBoundingBox(false);

    // rebuild normals for accurate gradient calculations
    terrain.RebuildNormals();

    // Use half grid size for gradient sampling
    double sampleDistance = gridSize * 0.5;

    // Pre-compute gradient map for the entire terrain
    GradientMap gradientMap = TerrainAnalysisUtils.PrecomputeGradientMap(terrain, bbox, gridSize, sampleDistance);

    // Show diagnostic information about gradient calculation methods used
    if (gradientMap.MeshNormalUsedCount > 0 || gradientMap.MeshNormalFallbackCount > 0) {
      string message = $"Gradient calculation statistics: ";

      if (gradientMap.MeshNormalUsedCount > 0) {
        message += $"{gradientMap.MeshNormalUsedCount} edge points used mesh normals  (primary method), ";
      }

      if (gradientMap.MeshNormalFallbackCount > 0) {
        message += $"{gradientMap.MeshNormalFallbackCount} points used fallback  method, ";
      }

      message += $"{gradientMap.FlatTerrainCount} flat terrain points. ";

      // Only show as warning if there are actual fallbacks (calculation
      // failures)
      if (gradientMap.MeshNormalFallbackCount > 0) {
        AddRuntimeMessage(
            GH_RuntimeMessageLevel.Warning,
            message +
                "Consider using a larger GridSize or expanding the mesh  boundary if edge accuracy is critical.");
      } else {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, message);
      }
    }

    // Create lists for output data
    List<Point3d> samplePoints = new List<Point3d>();
    List<double> slopeValues = new List<double>();
    List<Vector3d> slopeVectors = new List<Vector3d>();
    List<bool> inBounds = new List<bool>();

    // Collect data from gradient map
    for (int i = 0; i < gradientMap.XCount; i++) {
      for (int j = 0; j < gradientMap.YCount; j++) {
        double x = bbox.Min.X + (i + 0.5) * gridSize;
        double y = bbox.Min.Y + (j + 0.5) * gridSize;
        double z = gradientMap.Heights[i, j];

        // Skip invalid heights
        if (double.IsNaN(z)) {
          continue;
        }

        Point3d point = new Point3d(x, y, z);
        double slope = gradientMap.Slopes[i, j];
        Vector3d gradient = gradientMap.Gradients[i, j];

        // Calculate downhill direction vector (negative gradient, normalized)
        Vector3d downhillDir = Vector3d.Zero;
        if (gradient != Vector3d.Zero) {
          downhillDir = new Vector3d(-gradient.X, -gradient.Y, 0);
          downhillDir.Unitize();
          // Scale vector by grid size for visualization
          downhillDir *= gridSize;
        }

        // Check if slope is within bounds
        bool withinBounds = !double.IsNaN(slope) && slope >= slopeBounds.T0 && slope <= slopeBounds.T1;

        samplePoints.Add(point);
        slopeValues.Add(slope);
        slopeVectors.Add(downhillDir);
        inBounds.Add(withinBounds);
      }
    }

    // Create analysis mesh with color coding
    Mesh analysisMesh = CreateSlopeAnalysisMesh(gradientMap, bbox, gridSize, slopeBounds, selectedColorMap);

    // Output results
    DA.SetData("AnalysisMesh", analysisMesh);
    DA.SetDataList("SamplePoints", samplePoints);
    DA.SetDataList("SlopeValues", slopeValues);
    DA.SetDataList("SlopeVectors", slopeVectors);
    DA.SetDataList("InBounds", inBounds);
  }

  /// <summary>
  /// Create a color-coded mesh showing slope analysis
  /// </summary>
  private Mesh CreateSlopeAnalysisMesh(
      GradientMap gradientMap, BoundingBox bbox, double gridSize, Interval slopeBounds, string colorMapName) {
    Mesh analysisMesh = new Mesh();

    // Get the specified color map
    IColorMap colorMap = ColorMapHelper.GetColorMap(colorMapName);

    int xCount = gradientMap.XCount;
    int yCount = gradientMap.YCount;

    // Use the provided slope bounds for color mapping
    double minSlope = slopeBounds.T0;
    double maxSlope = slopeBounds.T1;

    // Create a mapping from grid indices to vertex indices (for handling NaN
    // values)
    int[,] vertexIndexMap = new int[xCount, yCount];
    for (int i = 0; i < xCount; i++) {
      for (int j = 0; j < yCount; j++) {
        vertexIndexMap[i, j] = -1;  // -1 means invalid/not added
      }
    }

    // Add vertices and colors (only for valid points)
    for (int i = 0; i < xCount; i++) {
      for (int j = 0; j < yCount; j++) {
        double x = bbox.Min.X + (i + 0.5) * gridSize;
        double y = bbox.Min.Y + (j + 0.5) * gridSize;
        double z = gradientMap.Heights[i, j];
        double slope = gradientMap.Slopes[i, j];

        // Only add vertex if height is valid
        if (!double.IsNaN(z)) {
          // Store the vertex index for this grid location
          vertexIndexMap[i, j] = analysisMesh.Vertices.Count;

          analysisMesh.Vertices.Add(x, y, z);

          // Map slope value to color
          System.Drawing.Color color;
          if (double.IsNaN(slope)) {
            color = System.Drawing.Color.White;  // White for invalid slope
          } else {
            // Clamp slope to bounds - values outside range get clamped to
            // min/max colors
            double clampedSlope = Math.Max(minSlope, Math.Min(maxSlope, slope));
            color = ColorMapHelper.MapValueToColor(clampedSlope, minSlope, maxSlope, colorMap);
          }

          analysisMesh.VertexColors.Add(color);
        }
      }
    }

    // Create triangular mesh faces (only for quads where all 4 vertices are
    // valid)
    for (int i = 0; i < xCount - 1; i++) {
      for (int j = 0; j < yCount - 1; j++) {
        // Get vertex indices for the quad
        int v0 = vertexIndexMap[i, j];
        int v1 = vertexIndexMap[i + 1, j];
        int v2 = vertexIndexMap[i + 1, j + 1];
        int v3 = vertexIndexMap[i, j + 1];

        // Only create faces if all 4 vertices are valid
        if (v0 >= 0 && v1 >= 0 && v2 >= 0 && v3 >= 0) {
          // Create two triangular faces
          analysisMesh.Faces.AddFace(v0, v1, v2);
          analysisMesh.Faces.AddFace(v0, v2, v3);
        }
      }
    }

    // Compact the mesh to remove any unused vertices
    analysisMesh.Compact();
    analysisMesh.Normals.ComputeNormals();

    return analysisMesh;
  }

  // Add right-click menu for color map selection (flat structure for
  // cross-platform)
  public override void AppendAdditionalMenuItems(ToolStripDropDown menu) {
    base.AppendAdditionalMenuItems(menu);

    // Add a separator
    Menu_AppendSeparator(menu);

    // Add color map header (disabled)
    var headerItem = Menu_AppendItem(menu, "Color Map:");
    headerItem.Enabled = false;

    // Add each color map as a direct menu item
    foreach (var colorMap in ColorMapHelper.AvailableColorMaps) {
      Menu_AppendItem(menu, "  " + colorMap.Name, OnColorMapSelected, true, colorMap.Name == selectedColorMap).Tag =
          colorMap.Name;
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

public class Watershed : GH_Component {
  private string selectedColorMap = "Viridis";

  public Watershed()
      : base("Watershed",
             "batWatershed",
             "Automatically segment terrain into watershed (drainage basin) regions. " +
                 "Each basin drains to a natural low point (sink) in the terrain.",
             "BAT",
             "Analysis") {}

  public override GH_Exposure Exposure => GH_Exposure.primary;
  public override Guid ComponentGuid => new Guid("3a5e9cd7-531d-4cf9-a185-f29b2f85a92c");
  protected override System.Drawing.Bitmap Icon => null;  // TODO: Add icon

  protected override void RegisterInputParams(GH_InputParamManager pManager) {
    pManager.AddMeshParameter("Mesh", "M", "Terrain mesh to analyze.", GH_ParamAccess.item);
    pManager.AddNumberParameter(
        "Resolution",
        "R",
        "Grid resolution for watershed analysis. If not provided, automatically derived from mesh vertex density.",
        GH_ParamAccess.item);
    pManager.AddNumberParameter(
        "MinArea",
        "A",
        "Minimum watershed area. Smaller watersheds will be merged into  neighbors. Set to 0 to keep all.",
        GH_ParamAccess.item,
        0.0);
    pManager.AddNumberParameter("EdgeTolerance",
                                "T",
                                "Distance from mesh boundary to detect edge watersheds (in grid cells). " +
                                    "Watersheds with sinks within this distance are merged with interior neighbors. " +
                                    "Default is 1.5 grid cells.",
                                GH_ParamAccess.item,
                                1.5);

    pManager[1].Optional = true;
    pManager[2].Optional = true;
    pManager[3].Optional = true;
  }

  protected override void RegisterOutputParams(GH_OutputParamManager pManager) {
    pManager.AddMeshParameter("WatershedMeshes",
                              "WM",
                              "Segmented meshes for each watershed basin (original mesh faces).",
                              GH_ParamAccess.list);
    pManager.AddCurveParameter(
        "Boundaries", "B", "3D boundary curves (ridgelines) extracted from mesh edges.", GH_ParamAccess.list);
    pManager.AddPointParameter(
        "Sinks", "S", "Natural outlet points (local minima) where water accumulates.", GH_ParamAccess.list);
  }

  protected override void SolveInstance(IGH_DataAccess DA) {
    Mesh terrain = null;
    double resolution = 0;
    double minArea = 0.0;
    double edgeTolerance = 1.5;

    if (!DA.GetData("Mesh", ref terrain))
      return;
    DA.GetData("Resolution", ref resolution);
    DA.GetData("MinArea", ref minArea);
    DA.GetData("EdgeTolerance", ref edgeTolerance);

    if (terrain == null) {
      AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input terrain mesh cannot be null");
      return;
    }

    // Validate edge tolerance
    if (edgeTolerance < 0) {
      AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "EdgeTolerance must be non-negative");
      edgeTolerance = Math.Max(0, edgeTolerance);
    }

    // Get bounding box
    BoundingBox bbox = terrain.GetBoundingBox(false);

    // Auto-derive resolution from mesh if not provided
    if (resolution <= 0) {
      resolution = DeriveResolutionFromMesh(terrain, bbox);
      AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, $"Auto-derived resolution: {resolution:F2} units");
    }

    // Rebuild normals
    terrain.RebuildNormals();

    // Compute gradient map using existing terrain analysis
    double sampleDistance = resolution * 0.5;
    GradientMap gradientMap = TerrainAnalysisUtils.PrecomputeGradientMap(terrain, bbox, resolution, sampleDistance);

    // Perform watershed segmentation using gradient-based flow
    var segmentation = SegmentWatersheds(gradientMap, bbox, resolution, minArea, edgeTolerance);

    if (segmentation.BasinCount == 0) {
      AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                        "No watersheds found. Terrain may be flat or have no  clear drainage patterns.");
      return;
    }

    // Assign mesh vertices to basins
    int[] vertexBasinIds = AssignVerticesToBasins(terrain, segmentation, gradientMap, bbox, resolution);

    // Extract mesh pieces for each basin
    List<Mesh> watershedMeshes = new List<Mesh>();
    List<Curve> boundaries = new List<Curve>();
    List<Point3d> sinks = new List<Point3d>();

    for (int basinId = 0; basinId < segmentation.BasinCount; basinId++) {
      // Extract original mesh faces for this basin
      var basinMesh = ExtractBasinFromOriginalMesh(terrain, vertexBasinIds, basinId);
      if (basinMesh != null && basinMesh.Vertices.Count > 0) {
        // Clean up isolated fragments - keep only the largest connected
        // component
        basinMesh = RemoveIsolatedFragments(basinMesh);

        if (basinMesh != null && basinMesh.Vertices.Count > 0) {
          watershedMeshes.Add(basinMesh);

          // Calculate sink point directly from the cleaned watershed mesh
          // This ensures the sink is within the actual mesh after fragment
          // removal
          Point3d sinkPoint = FindLowestPointInMesh(basinMesh);
          sinks.Add(sinkPoint);
        }
      }
    }

    // Extract boundaries DIRECTLY from each watershed mesh using naked edges
    // This ensures boundaries are correct after all merging operations
    for (int i = 0; i < watershedMeshes.Count; i++) {
      var boundary = ExtractBoundaryFromWatershedMesh(watershedMeshes[i]);
      if (boundary != null) {
        boundaries.Add(boundary);
      }
    }

    // Output
    DA.SetDataList("WatershedMeshes", watershedMeshes);
    DA.SetDataList("Boundaries", boundaries);
    DA.SetDataList("Sinks", sinks);

    AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, $"Found {segmentation.BasinCount} watershed basins");
  }

  /// <summary>
  /// Auto-derive resolution from mesh vertex density
  /// </summary>
  private double DeriveResolutionFromMesh(Mesh mesh, BoundingBox bbox) {
    if (mesh.Vertices.Count == 0) {
      return 10.0;
    }

    double totalEdgeLength = 0;
    int edgeCount = 0;

    for (int i = 0; i < Math.Min(mesh.Faces.Count, 1000); i++) {
      var face = mesh.Faces[i];
      Point3f v0 = mesh.Vertices[face.A];
      Point3f v1 = mesh.Vertices[face.B];
      Point3f v2 = mesh.Vertices[face.C];

      totalEdgeLength += v0.DistanceTo(v1);
      totalEdgeLength += v1.DistanceTo(v2);
      totalEdgeLength += v2.DistanceTo(v0);
      edgeCount += 3;
    }

    double avgEdgeLength = totalEdgeLength / edgeCount;
    return avgEdgeLength * 2.5;
  }

  /// <summary>
  /// Segment terrain into watersheds using gradient-based flow simulation
  /// </summary>
  private WatershedSegmentation
  SegmentWatersheds(GradientMap gradientMap, BoundingBox bbox, double gridSize, double minArea, double edgeTolerance) {
    int xCount = gradientMap.XCount;
    int yCount = gradientMap.YCount;

    int[,] basinIds = new int[xCount, yCount];
    for (int i = 0; i < xCount; i++) {
      for (int j = 0; j < yCount; j++) {
        basinIds[i, j] = -1;
      }
    }

    // Step 1: Find drainage targets
    (int i, int j)[,] drainageTarget = new(int, int)[xCount, yCount];

    for (int i = 0; i < xCount; i++) {
      for (int j = 0; j < yCount; j++) {
        if (double.IsNaN(gradientMap.Heights[i, j])) {
          drainageTarget[i, j] = (-1, -1);
          continue;
        }

        Point3d cellCenter = new Point3d(
            bbox.Min.X + (i + 0.5) * gridSize, bbox.Min.Y + (j + 0.5) * gridSize, gradientMap.Heights[i, j]);

        Vector3d gradient = gradientMap.Gradients[i, j];

        if (gradient == Vector3d.Zero || gradient.Length < 0.001) {
          drainageTarget[i, j] = (i, j);
        } else {
          Vector3d flowDir = new Vector3d(-gradient.X, -gradient.Y, 0);
          flowDir.Unitize();

          double targetX = cellCenter.X + flowDir.X * gridSize;
          double targetY = cellCenter.Y + flowDir.Y * gridSize;

          int targetI = (int)Math.Floor((targetX - bbox.Min.X) / gridSize);
          int targetJ = (int)Math.Floor((targetY - bbox.Min.Y) / gridSize);

          targetI = Math.Max(0, Math.Min(xCount - 1, targetI));
          targetJ = Math.Max(0, Math.Min(yCount - 1, targetJ));

          if (!double.IsNaN(gradientMap.Heights[targetI, targetJ]) &&
              gradientMap.Heights[targetI, targetJ] < gradientMap.Heights[i, j]) {
            drainageTarget[i, j] = (targetI, targetJ);
          } else {
            drainageTarget[i, j] = (i, j);
          }
        }
      }
    }

    // Step 2: Trace to sinks and assign basin IDs
    Dictionary<(int, int), int> sinkToBasinId = new Dictionary<(int, int), int>();
    int currentBasinId = 0;

    for (int i = 0; i < xCount; i++) {
      for (int j = 0; j < yCount; j++) {
        if (double.IsNaN(gradientMap.Heights[i, j]) || basinIds[i, j] >= 0) {
          continue;
        }

        var sink = TraceToDrainageSink(i, j, drainageTarget, xCount, yCount);

        if (sink.i < 0) {
          continue;
        }

        if (!sinkToBasinId.ContainsKey((sink.i, sink.j))) {
          sinkToBasinId[(sink.i, sink.j)] = currentBasinId++;
        }

        basinIds[i, j] = sinkToBasinId[(sink.i, sink.j)];
      }
    }

    // Step 3: Merge small basins
    if (minArea > 0) {
      MergeSmallBasins(basinIds, gradientMap, gridSize, minArea, ref currentBasinId);
    }

    // Step 4: Merge edge basins (watersheds with sinks on boundaries)
    // These represent incomplete catchments where water flows off the mesh
    // Use cluster-based merging for efficiency and stability
    MergeEdgeBasinsClustered(basinIds, gradientMap, xCount, yCount, gridSize, edgeTolerance, ref currentBasinId);

    return new WatershedSegmentation {
      BasinIds = basinIds, BasinCount = currentBasinId, SinkToBasinId = sinkToBasinId, XCount = xCount, YCount = yCount
    };
  }

  private (int i, int j)
      TraceToDrainageSink(int startI, int startJ, (int i, int j)[,] drainageTarget, int xCount, int yCount) {
    int currentI = startI;
    int currentJ = startJ;
    HashSet<(int, int)> visited = new HashSet<(int, int)>();
    int maxSteps = xCount * yCount;

    for (int step = 0; step < maxSteps; step++) {
      if (currentI < 0 || currentI >= xCount || currentJ < 0 || currentJ >= yCount) {
        return (-1, -1);
      }

      if (visited.Contains((currentI, currentJ))) {
        return (currentI, currentJ);
      }
      visited.Add((currentI, currentJ));

      var target = drainageTarget[currentI, currentJ];

      if (target.i == currentI && target.j == currentJ) {
        return (currentI, currentJ);
      }

      currentI = target.i;
      currentJ = target.j;
    }

    return (-1, -1);
  }

  private void
  MergeSmallBasins(int[,] basinIds, GradientMap gradientMap, double gridSize, double minArea, ref int basinCount) {
    int xCount = basinIds.GetLength(0);
    int yCount = basinIds.GetLength(1);
    double cellArea = gridSize * gridSize;

    bool merged = true;
    while (merged) {
      merged = false;

      Dictionary<int, int> basinSizes = new Dictionary<int, int>();
      for (int i = 0; i < xCount; i++) {
        for (int j = 0; j < yCount; j++) {
          int basinId = basinIds[i, j];
          if (basinId >= 0) {
            basinSizes[basinId] = basinSizes.GetValueOrDefault(basinId, 0) + 1;
          }
        }
      }

      foreach (var kvp in basinSizes) {
        int basinId = kvp.Key;
        double area = kvp.Value * cellArea;

        if (area < minArea) {
          int bestNeighbor = FindBestNeighborBasin(basinIds, basinId, gradientMap, xCount, yCount);

          if (bestNeighbor >= 0) {
            for (int i = 0; i < xCount; i++) {
              for (int j = 0; j < yCount; j++) {
                if (basinIds[i, j] == basinId) {
                  basinIds[i, j] = bestNeighbor;
                }
              }
            }
            merged = true;
            break;
          }
        }
      }
    }
  }

  private int FindBestNeighborBasin(int[,] basinIds, int basinId, GradientMap gradientMap, int xCount, int yCount) {
    Dictionary<int, double> neighborMinHeights = new Dictionary<int, double>();

    for (int i = 0; i < xCount; i++) {
      for (int j = 0; j < yCount; j++) {
        if (basinIds[i, j] != basinId) {
          continue;
        }

        for (int di = -1; di <= 1; di++) {
          for (int dj = -1; dj <= 1; dj++) {
            if (di == 0 && dj == 0)
              continue;

            int ni = i + di;
            int nj = j + dj;

            if (ni >= 0 && ni < xCount && nj >= 0 && nj < yCount) {
              int neighborBasin = basinIds[ni, nj];

              if (neighborBasin >= 0 && neighborBasin != basinId) {
                double height = gradientMap.Heights[i, j];

                if (!neighborMinHeights.ContainsKey(neighborBasin) || height < neighborMinHeights[neighborBasin]) {
                  neighborMinHeights[neighborBasin] = height;
                }
              }
            }
          }
        }
      }
    }

    if (neighborMinHeights.Count == 0) {
      return -1;
    }

    return neighborMinHeights.OrderBy(kvp => kvp.Value).First().Key;
  }

  /// <summary>
  /// Merge edge basins using cluster-based approach (efficient single-pass)
  /// Edge basins are watersheds with sinks near boundaries - incomplete catchments
  /// where water flows off the mesh
  /// </summary>
  private void MergeEdgeBasinsClustered(int[,] basinIds,
                                        GradientMap gradientMap,
                                        int xCount,
                                        int yCount,
                                        double gridSize,
                                        double edgeTolerance,
                                        ref int basinCount) {
    // Step 1: Detect all edge basins (sinks near boundaries)
    HashSet<int> edgeBasinIds = DetectEdgeBasins(basinIds, gradientMap, xCount, yCount, edgeTolerance);

    if (edgeBasinIds.Count == 0) {
      return;  // No edge basins to merge
    }

    // Step 2: Build merge graph - for each edge basin, find best interior neighbor
    Dictionary<int, int> edgeBasinToInteriorNeighbor = new Dictionary<int, int>();

    foreach (int edgeBasinId in edgeBasinIds) {
      int bestNeighbor = FindBestInteriorNeighbor(basinIds, edgeBasinId, edgeBasinIds, gradientMap, xCount, yCount);

      if (bestNeighbor >= 0) {
        edgeBasinToInteriorNeighbor[edgeBasinId] = bestNeighbor;
      }
    }

    // Step 3: Build merge clusters using union-find
    // Multiple edge basins may merge into same interior basin
    // Or edge basins may form chains (edge1→edge2→interior)
    Dictionary<int, int> basinToClusterRoot = BuildMergeClusters(edgeBasinToInteriorNeighbor, edgeBasinIds);

    // Step 4: Apply merging - remap all basin IDs in single pass
    for (int i = 0; i < xCount; i++) {
      for (int j = 0; j < yCount; j++) {
        int originalBasinId = basinIds[i, j];
        if (originalBasinId >= 0 && basinToClusterRoot.ContainsKey(originalBasinId)) {
          basinIds[i, j] = basinToClusterRoot[originalBasinId];
        }
      }
    }

    // Step 5: Renumber basins to be sequential (remove gaps from merged basins)
    RenumberBasins(basinIds, xCount, yCount, ref basinCount);
  }

  /// <summary>
  /// Detect basins whose sinks are on or near the mesh boundary
  /// Returns set of edge basin IDs
  /// </summary>
  private HashSet<int>
  DetectEdgeBasins(int[,] basinIds, GradientMap gradientMap, int xCount, int yCount, double edgeTolerance) {
    HashSet<int> edgeBasinIds = new HashSet<int>();

    // Get all unique basin IDs currently in the grid
    HashSet<int> activeBasinIds = new HashSet<int>();
    for (int i = 0; i < xCount; i++) {
      for (int j = 0; j < yCount; j++) {
        if (basinIds[i, j] >= 0) {
          activeBasinIds.Add(basinIds[i, j]);
        }
      }
    }

    // For each basin, find its current sink location
    foreach (int basinId in activeBasinIds) {
      double minHeight = double.MaxValue;
      int sinkI = -1, sinkJ = -1;

      // Find the lowest point in this basin (current sink)
      for (int i = 0; i < xCount; i++) {
        for (int j = 0; j < yCount; j++) {
          if (basinIds[i, j] == basinId) {
            double height = gradientMap.Heights[i, j];
            if (!double.IsNaN(height) && height < minHeight) {
              minHeight = height;
              sinkI = i;
              sinkJ = j;
            }
          }
        }
      }

      // Check if this sink is near any boundary
      if (sinkI >= 0 && sinkJ >= 0) {
        // Calculate distance to nearest boundary (in grid cells)
        int distToLeft = sinkI;
        int distToRight = (xCount - 1) - sinkI;
        int distToBottom = sinkJ;
        int distToTop = (yCount - 1) - sinkJ;

        int minDistToBoundary = Math.Min(Math.Min(distToLeft, distToRight), Math.Min(distToBottom, distToTop));

        if (minDistToBoundary <= edgeTolerance) {
          edgeBasinIds.Add(basinId);
        }
      }
    }

    return edgeBasinIds;
  }

  /// <summary>
  /// Find the best interior (non-edge) neighbor for an edge basin
  /// Uses proximity to sink location as primary criterion
  /// </summary>
  private int FindBestInteriorNeighbor(
      int[,] basinIds, int edgeBasinId, HashSet<int> edgeBasinIds, GradientMap gradientMap, int xCount, int yCount) {
    // Find the sink location for this edge basin
    double minHeight = double.MaxValue;
    int sinkI = -1, sinkJ = -1;

    for (int i = 0; i < xCount; i++) {
      for (int j = 0; j < yCount; j++) {
        if (basinIds[i, j] == edgeBasinId) {
          double height = gradientMap.Heights[i, j];
          if (!double.IsNaN(height) && height < minHeight) {
            minHeight = height;
            sinkI = i;
            sinkJ = j;
          }
        }
      }
    }

    if (sinkI < 0 || sinkJ < 0) {
      return -1;  // No valid sink found
    }

    // Search for neighboring basins in expanding radius around sink
    int searchRadius = 1;
    int maxSearchRadius = 10;

    while (searchRadius <= maxSearchRadius) {
      // Check all cells within current search radius
      for (int di = -searchRadius; di <= searchRadius; di++) {
        for (int dj = -searchRadius; dj <= searchRadius; dj++) {
          // Only check cells on the perimeter of the search radius
          if (Math.Abs(di) != searchRadius && Math.Abs(dj) != searchRadius) {
            continue;
          }

          int checkI = sinkI + di;
          int checkJ = sinkJ + dj;

          // Check bounds
          if (checkI < 0 || checkI >= xCount || checkJ < 0 || checkJ >= yCount) {
            continue;
          }

          int neighborBasinId = basinIds[checkI, checkJ];

          // Found a different basin near the sink
          // Prefer interior basins (non-edge) but accept edge basins as fallback
          if (neighborBasinId >= 0 && neighborBasinId != edgeBasinId) {
            // Return first interior basin found (best match)
            if (!edgeBasinIds.Contains(neighborBasinId)) {
              return neighborBasinId;
            }
          }
        }
      }

      searchRadius++;
    }

    // Fallback: if no interior neighbor found, find ANY neighbor using lowest boundary
    return FindBestNeighborBasin(basinIds, edgeBasinId, gradientMap, xCount, yCount);
  }

  /// <summary>
  /// Build merge clusters using union-find approach
  /// Handles complex cases like edge1→edge2→interior or multiple edges→same interior
  /// Returns mapping from original basin ID to cluster root (final basin ID)
  /// </summary>
  private Dictionary<int, int> BuildMergeClusters(Dictionary<int, int> edgeBasinToInteriorNeighbor,
                                                  HashSet<int> edgeBasinIds) {
    // Union-Find data structure
    Dictionary<int, int> parent = new Dictionary<int, int>();

    // Initialize: each basin is its own parent
    foreach (var kvp in edgeBasinToInteriorNeighbor) {
      if (!parent.ContainsKey(kvp.Key)) {
        parent[kvp.Key] = kvp.Key;
      }
      if (!parent.ContainsKey(kvp.Value)) {
        parent[kvp.Value] = kvp.Value;
      }
    }

    // Union operation: merge edge basin with its target
    foreach (var kvp in edgeBasinToInteriorNeighbor) {
      int edgeBasin = kvp.Key;
      int targetBasin = kvp.Value;

      Union(parent, edgeBasin, targetBasin);
    }

    // Find root for all basins (with path compression)
    Dictionary<int, int> basinToClusterRoot = new Dictionary<int, int>();
    foreach (var basinId in parent.Keys) {
      basinToClusterRoot[basinId] = Find(parent, basinId);
    }

    // Ensure interior basins (non-edge) become cluster roots
    // If an interior basin merged with another interior basin, keep one as root
    foreach (var kvp in basinToClusterRoot.ToList()) {
      int basinId = kvp.Key;
      int root = kvp.Value;

      // If this is an edge basin, make sure its root is an interior basin
      if (edgeBasinIds.Contains(basinId) && edgeBasinIds.Contains(root)) {
        // Both are edge basins - need to find the ultimate interior target
        // This handles chains like edge1→edge2→interior
        int ultimateTarget = FindUltimateInteriorTarget(basinId, edgeBasinToInteriorNeighbor, edgeBasinIds);
        if (ultimateTarget >= 0) {
          basinToClusterRoot[basinId] = ultimateTarget;
        }
      }
    }

    return basinToClusterRoot;
  }

  /// <summary>
  /// Find operation for union-find (with path compression)
  /// </summary>
  private int Find(Dictionary<int, int> parent, int id) {
    if (parent[id] != id) {
      parent[id] = Find(parent, parent[id]);  // Path compression
    }
    return parent[id];
  }

  /// <summary>
  /// Union operation for union-find
  /// </summary>
  private void Union(Dictionary<int, int> parent, int id1, int id2) {
    int root1 = Find(parent, id1);
    int root2 = Find(parent, id2);

    if (root1 != root2) {
      // Always make the larger ID the parent (arbitrary choice for consistency)
      if (root1 < root2) {
        parent[root1] = root2;
      } else {
        parent[root2] = root1;
      }
    }
  }

  /// <summary>
  /// Find the ultimate interior basin target by following merge chain
  /// Handles cases like edge1→edge2→edge3→interior
  /// </summary>
  private int FindUltimateInteriorTarget(int startBasinId,
                                         Dictionary<int, int> edgeBasinToInteriorNeighbor,
                                         HashSet<int> edgeBasinIds) {
    HashSet<int> visited = new HashSet<int>();
    int currentBasin = startBasinId;
    int maxSteps = 100;  // Safety limit

    for (int step = 0; step < maxSteps; step++) {
      if (visited.Contains(currentBasin)) {
        return -1;  // Circular reference detected
      }
      visited.Add(currentBasin);

      // If this basin is interior (not edge), we found the target
      if (!edgeBasinIds.Contains(currentBasin)) {
        return currentBasin;
      }

      // Follow the merge chain
      if (edgeBasinToInteriorNeighbor.ContainsKey(currentBasin)) {
        currentBasin = edgeBasinToInteriorNeighbor[currentBasin];
      } else {
        return -1;  // No target found
      }
    }

    return -1;  // Max steps exceeded
  }

  /// <summary>
  /// Renumber basins to be sequential (remove gaps from merged basins)
  /// </summary>
  private void RenumberBasins(int[,] basinIds, int xCount, int yCount, ref int basinCount) {
    Dictionary<int, int> oldToNew = new Dictionary<int, int>();
    int newBasinId = 0;

    for (int i = 0; i < xCount; i++) {
      for (int j = 0; j < yCount; j++) {
        int oldId = basinIds[i, j];
        if (oldId >= 0) {
          if (!oldToNew.ContainsKey(oldId)) {
            oldToNew[oldId] = newBasinId++;
          }
          basinIds[i, j] = oldToNew[oldId];
        }
      }
    }

    basinCount = newBasinId;
  }

  /// <summary>
  /// Assign each mesh vertex to a watershed basin based on its location
  /// </summary>
  private int[] AssignVerticesToBasins(
      Mesh mesh, WatershedSegmentation segmentation, GradientMap gradientMap, BoundingBox bbox, double gridSize) {
    int[] vertexBasinIds = new int[mesh.Vertices.Count];

    for (int v = 0; v < mesh.Vertices.Count; v++) {
      Point3d vertex = mesh.Vertices[v];

      // Find grid cell for this vertex
      int i = (int)Math.Floor((vertex.X - bbox.Min.X) / gridSize);
      int j = (int)Math.Floor((vertex.Y - bbox.Min.Y) / gridSize);

      // Clamp to valid range
      i = Math.Max(0, Math.Min(segmentation.XCount - 1, i));
      j = Math.Max(0, Math.Min(segmentation.YCount - 1, j));

      vertexBasinIds[v] = segmentation.BasinIds[i, j];
    }

    return vertexBasinIds;
  }

  /// <summary>
  /// Extract mesh faces belonging to a specific basin (preserves original mesh
  /// geometry)
  /// </summary>
  private Mesh ExtractBasinFromOriginalMesh(Mesh originalMesh, int[] vertexBasinIds, int basinId) {
    Mesh basinMesh = new Mesh();

    // Map from original vertex indices to new basin mesh vertex indices
    Dictionary<int, int> vertexMap = new Dictionary<int, int>();

    System.Drawing.Color basinColor = GetWatershedColor(basinId);

    // Add faces where all vertices belong to this basin
    for (int f = 0; f < originalMesh.Faces.Count; f++) {
      MeshFace face = originalMesh.Faces[f];

      // Check if all vertices of this face belong to the basin
      bool allInBasin = true;
      int[] faceVertices =
          face.IsQuad ? new int[] { face.A, face.B, face.C, face.D } : new int[] { face.A, face.B, face.C };

      foreach (int vi in faceVertices) {
        if (vi >= vertexBasinIds.Length || vertexBasinIds[vi] != basinId) {
          allInBasin = false;
          break;
        }
      }

      if (allInBasin) {
        // Add vertices if not already added
        int[] newIndices = new int[faceVertices.Length];

        for (int i = 0; i < faceVertices.Length; i++) {
          int origIndex = faceVertices[i];

          if (!vertexMap.ContainsKey(origIndex)) {
            vertexMap[origIndex] = basinMesh.Vertices.Count;
            basinMesh.Vertices.Add(originalMesh.Vertices[origIndex]);
            basinMesh.VertexColors.Add(basinColor);
          }

          newIndices[i] = vertexMap[origIndex];
        }

        // Add face
        if (face.IsQuad) {
          basinMesh.Faces.AddFace(newIndices[0], newIndices[1], newIndices[2], newIndices[3]);
        } else {
          basinMesh.Faces.AddFace(newIndices[0], newIndices[1], newIndices[2]);
        }
      }
    }

    if (basinMesh.Vertices.Count > 0) {
      basinMesh.Normals.ComputeNormals();
      basinMesh.Compact();
    }

    return basinMesh;
  }

  /// <summary>
  /// Remove isolated mesh fragments, keeping only the largest connected
  /// component This cleans up small disconnected pieces that result from basin
  /// merging
  /// </summary>
  private Mesh RemoveIsolatedFragments(Mesh mesh) {
    if (mesh == null || mesh.Faces.Count == 0) {
      return mesh;
    }

    // Build face connectivity graph
    Dictionary<int, HashSet<int>> faceNeighbors = new Dictionary<int, HashSet<int>>();

    // For each face, find neighboring faces (faces that share an edge)
    for (int f = 0; f < mesh.Faces.Count; f++) {
      faceNeighbors[f] = new HashSet<int>();
    }

    // Build edge-to-faces mapping
    Dictionary<(int, int), List<int>> edgeToFaces = new Dictionary<(int, int), List<int>>();

    for (int f = 0; f < mesh.Faces.Count; f++) {
      MeshFace face = mesh.Faces[f];

      // Get face edges
      List<(int v0, int v1)> edges = new List<(int, int)>();
      if (face.IsQuad) {
        edges.Add((face.A, face.B));
        edges.Add((face.B, face.C));
        edges.Add((face.C, face.D));
        edges.Add((face.D, face.A));
      } else {
        edges.Add((face.A, face.B));
        edges.Add((face.B, face.C));
        edges.Add((face.C, face.A));
      }

      foreach (var (v0, v1) in edges) {
        // Ensure consistent edge ordering
        int minV = Math.Min(v0, v1);
        int maxV = Math.Max(v0, v1);
        var edgeKey = (minV, maxV);

        if (!edgeToFaces.ContainsKey(edgeKey)) {
          edgeToFaces[edgeKey] = new List<int>();
        }
        edgeToFaces[edgeKey].Add(f);
      }
    }

    // Build face neighbor relationships
    foreach (var facesOnEdge in edgeToFaces.Values) {
      if (facesOnEdge.Count >= 2) {
        // Faces sharing an edge are neighbors
        for (int i = 0; i < facesOnEdge.Count; i++) {
          for (int j = i + 1; j < facesOnEdge.Count; j++) {
            faceNeighbors[facesOnEdge[i]].Add(facesOnEdge[j]);
            faceNeighbors[facesOnEdge[j]].Add(facesOnEdge[i]);
          }
        }
      }
    }

    // Find connected components using BFS
    List<HashSet<int>> components = new List<HashSet<int>>();
    HashSet<int> visited = new HashSet<int>();

    for (int startFace = 0; startFace < mesh.Faces.Count; startFace++) {
      if (visited.Contains(startFace)) {
        continue;
      }

      // BFS to find all faces in this connected component
      HashSet<int> component = new HashSet<int>();
      Queue<int> queue = new Queue<int>();
      queue.Enqueue(startFace);
      visited.Add(startFace);

      while (queue.Count > 0) {
        int face = queue.Dequeue();
        component.Add(face);

        foreach (int neighbor in faceNeighbors[face]) {
          if (!visited.Contains(neighbor)) {
            visited.Add(neighbor);
            queue.Enqueue(neighbor);
          }
        }
      }

      components.Add(component);
    }

    // If there's only one component, no cleanup needed
    if (components.Count <= 1) {
      return mesh;
    }

    // Find the largest component
    HashSet<int> largestComponent = components.OrderByDescending(c => c.Count).First();

    // Create new mesh with only faces from the largest component
    Mesh cleanedMesh = new Mesh();
    Dictionary<int, int> oldToNewVertexIndex = new Dictionary<int, int>();

    foreach (int faceIndex in largestComponent) {
      MeshFace face = mesh.Faces[faceIndex];

      int[] faceVertices =
          face.IsQuad ? new int[] { face.A, face.B, face.C, face.D } : new int[] { face.A, face.B, face.C };

      int[] newIndices = new int[faceVertices.Length];

      for (int i = 0; i < faceVertices.Length; i++) {
        int oldIndex = faceVertices[i];

        if (!oldToNewVertexIndex.ContainsKey(oldIndex)) {
          oldToNewVertexIndex[oldIndex] = cleanedMesh.Vertices.Count;
          cleanedMesh.Vertices.Add(mesh.Vertices[oldIndex]);

          // Preserve vertex colors if they exist
          if (oldIndex < mesh.VertexColors.Count) {
            cleanedMesh.VertexColors.Add(mesh.VertexColors[oldIndex]);
          }
        }

        newIndices[i] = oldToNewVertexIndex[oldIndex];
      }

      // Add face
      if (face.IsQuad) {
        cleanedMesh.Faces.AddFace(newIndices[0], newIndices[1], newIndices[2], newIndices[3]);
      } else {
        cleanedMesh.Faces.AddFace(newIndices[0], newIndices[1], newIndices[2]);
      }
    }

    if (cleanedMesh.Vertices.Count > 0) {
      cleanedMesh.Normals.ComputeNormals();
      cleanedMesh.Compact();
    }

    return cleanedMesh;
  }

  /// <summary>
  /// Find the lowest point (sink) in a watershed mesh
  /// </summary>
  private Point3d FindLowestPointInMesh(Mesh mesh) {
    if (mesh == null || mesh.Vertices.Count == 0) {
      return Point3d.Unset;
    }

    Point3d lowestPoint = mesh.Vertices[0];
    double minZ = lowestPoint.Z;

    for (int i = 1; i < mesh.Vertices.Count; i++) {
      Point3d vertex = mesh.Vertices[i];
      if (vertex.Z < minZ) {
        minZ = vertex.Z;
        lowestPoint = vertex;
      }
    }

    return lowestPoint;
  }

  /// <summary>
  /// Extract boundary curve from watershed mesh by finding naked edges
  /// (boundary edges)
  /// </summary>
  private Curve ExtractBoundaryFromWatershedMesh(Mesh watershedMesh) {
    if (watershedMesh == null || watershedMesh.Vertices.Count < 3) {
      return null;
    }

    // Find naked edges (edges that belong to only one face = boundary)
    Dictionary<(int, int), int> edgeCount = new Dictionary<(int, int), int>();
    Dictionary<(int, int), Line> edgeGeometry = new Dictionary<(int, int), Line>();

    // Count each edge
    for (int f = 0; f < watershedMesh.Faces.Count; f++) {
      MeshFace face = watershedMesh.Faces[f];

      // Get face edges
      List<(int v0, int v1)> edges = new List<(int, int)>();
      if (face.IsQuad) {
        edges.Add((face.A, face.B));
        edges.Add((face.B, face.C));
        edges.Add((face.C, face.D));
        edges.Add((face.D, face.A));
      } else {
        edges.Add((face.A, face.B));
        edges.Add((face.B, face.C));
        edges.Add((face.C, face.A));
      }

      // Count each edge
      foreach (var (v0, v1) in edges) {
        // Ensure consistent edge ordering
        int minV = Math.Min(v0, v1);
        int maxV = Math.Max(v0, v1);
        var edgeKey = (minV, maxV);

        if (!edgeCount.ContainsKey(edgeKey)) {
          edgeCount[edgeKey] = 0;
          edgeGeometry[edgeKey] = new Line(watershedMesh.Vertices[v0], watershedMesh.Vertices[v1]);
        }
        edgeCount[edgeKey]++;
      }
    }

    // Collect naked edges (edges with count == 1 = boundary edges)
    List<Line> boundaryEdges = new List<Line>();
    foreach (var kvp in edgeCount) {
      if (kvp.Value == 1) {  // Naked edge = boundary
        boundaryEdges.Add(edgeGeometry[kvp.Key]);
      }
    }

    if (boundaryEdges.Count == 0) {
      return null;
    }

    // Connect boundary edges into a continuous path
    List<Point3d> boundaryPoints = ConnectBoundaryEdges(boundaryEdges);

    if (boundaryPoints.Count >= 3) {
      // Close the curve if not already closed
      if (boundaryPoints[0].DistanceTo(boundaryPoints[boundaryPoints.Count - 1]) > 0.001) {
        boundaryPoints.Add(boundaryPoints[0]);
      }

      Polyline polyline = new Polyline(boundaryPoints);
      return polyline.ToNurbsCurve();
    }

    return null;
  }

  /// <summary>
  /// Extract boundary edges from mesh where vertices belong to different basins
  /// </summary>
  private Curve ExtractBoundaryFromMeshEdges(Mesh mesh, int[] vertexBasinIds, int basinId) {
    // Find all boundary edges (edges where one vertex is in basin, other is
    // not)
    List<Line> boundaryEdges = new List<Line>();
    HashSet<(int, int)> processedEdges = new HashSet<(int, int)>();

    for (int f = 0; f < mesh.Faces.Count; f++) {
      MeshFace face = mesh.Faces[f];

      // Get face edges
      List<(int v0, int v1)> edges = new List<(int, int)>();
      if (face.IsQuad) {
        edges.Add((face.A, face.B));
        edges.Add((face.B, face.C));
        edges.Add((face.C, face.D));
        edges.Add((face.D, face.A));
      } else {
        edges.Add((face.A, face.B));
        edges.Add((face.B, face.C));
        edges.Add((face.C, face.A));
      }

      // Check each edge
      foreach (var (v0, v1) in edges) {
        // Ensure consistent edge ordering
        int minV = Math.Min(v0, v1);
        int maxV = Math.Max(v0, v1);

        if (processedEdges.Contains((minV, maxV))) {
          continue;
        }

        int basin0 = vertexBasinIds[v0];
        int basin1 = vertexBasinIds[v1];

        // Boundary edge: one vertex in basin, other not
        if ((basin0 == basinId && basin1 != basinId) || (basin1 == basinId && basin0 != basinId)) {
          boundaryEdges.Add(new Line(mesh.Vertices[v0], mesh.Vertices[v1]));
          processedEdges.Add((minV, maxV));
        }
      }
    }

    if (boundaryEdges.Count == 0) {
      return null;
    }

    // Connect boundary edges into a polyline
    List<Point3d> boundaryPoints = ConnectBoundaryEdges(boundaryEdges);

    if (boundaryPoints.Count >= 3) {
      // Close the curve
      if (boundaryPoints[0].DistanceTo(boundaryPoints[boundaryPoints.Count - 1]) > 0.001) {
        boundaryPoints.Add(boundaryPoints[0]);
      }

      Polyline polyline = new Polyline(boundaryPoints);
      return polyline.ToNurbsCurve();
    }

    return null;
  }

  /// <summary>
  /// Connect boundary edges into a continuous path
  /// </summary>
  private List<Point3d> ConnectBoundaryEdges(List<Line> edges) {
    if (edges.Count == 0) {
      return new List<Point3d>();
    }

    List<Point3d> path = new List<Point3d>();
    HashSet<int> usedEdges = new HashSet<int>();

    // Start with first edge
    path.Add(edges[0].From);
    path.Add(edges[0].To);
    usedEdges.Add(0);

    // Greedily connect edges
    while (usedEdges.Count < edges.Count) {
      Point3d lastPoint = path[path.Count - 1];
      double minDist = double.MaxValue;
      int bestEdge = -1;
      bool reverseEdge = false;

      // Find closest unused edge
      for (int i = 0; i < edges.Count; i++) {
        if (usedEdges.Contains(i)) {
          continue;
        }

        double distFrom = lastPoint.DistanceTo(edges[i].From);
        double distTo = lastPoint.DistanceTo(edges[i].To);

        if (distFrom < minDist) {
          minDist = distFrom;
          bestEdge = i;
          reverseEdge = false;
        }

        if (distTo < minDist) {
          minDist = distTo;
          bestEdge = i;
          reverseEdge = true;
        }
      }

      if (bestEdge < 0 || minDist > 0.1) {
        break;  // Can't connect more edges
      }

      // Add next edge
      if (reverseEdge) {
        path.Add(edges[bestEdge].From);
      } else {
        path.Add(edges[bestEdge].To);
      }

      usedEdges.Add(bestEdge);
    }

    return path;
  }

  private System.Drawing.Color GetWatershedColor(int basinId) {
    double hue = (basinId * 137.5) % 360.0;
    double saturation = 0.7;
    double value = 0.9;

    double c = value * saturation;
    double x = c * (1 - Math.Abs((hue / 60.0) % 2 - 1));
    double m = value - c;

    double r = 0, g = 0, b = 0;
    if (hue < 60) {
      r = c;
      g = x;
    } else if (hue < 120) {
      r = x;
      g = c;
    } else if (hue < 180) {
      g = c;
      b = x;
    } else if (hue < 240) {
      g = x;
      b = c;
    } else if (hue < 300) {
      r = x;
      b = c;
    } else {
      r = c;
      b = x;
    }

    return System.Drawing.Color.FromArgb(
        (int)Math.Round((r + m) * 255), (int)Math.Round((g + m) * 255), (int)Math.Round((b + m) * 255));
  }

  public override void AppendAdditionalMenuItems(ToolStripDropDown menu) {
    base.AppendAdditionalMenuItems(menu);
  }

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

/// <summary>
/// Watershed segmentation result
/// </summary>
public class WatershedSegmentation {
  public int[,] BasinIds { get; set; }
  public int BasinCount { get; set; }
  public Dictionary<(int, int), int> SinkToBasinId { get; set; }
  public int XCount { get; set; }
  public int YCount { get; set; }

  public int GetBasinCellCount(int basinId) {
    int count = 0;
    for (int i = 0; i < XCount; i++) {
      for (int j = 0; j < YCount; j++) {
        if (BasinIds[i, j] == basinId) {
          count++;
        }
      }
    }
    return count;
  }

  public Point3d GetSinkPoint(int basinId, GradientMap gradientMap, BoundingBox bbox, double gridSize) {
    double minHeight = double.MaxValue;
    int sinkI = -1, sinkJ = -1;

    for (int i = 0; i < XCount; i++) {
      for (int j = 0; j < YCount; j++) {
        if (BasinIds[i, j] == basinId) {
          double height = gradientMap.Heights[i, j];
          if (!double.IsNaN(height) && height < minHeight) {
            minHeight = height;
            sinkI = i;
            sinkJ = j;
          }
        }
      }
    }

    if (sinkI >= 0) {
      return new Point3d(bbox.Min.X + (sinkI + 0.5) * gridSize, bbox.Min.Y + (sinkJ + 0.5) * gridSize, minHeight);
    }

    return Point3d.Unset;
  }
}
}  // namespace BeingAliveTerrain
