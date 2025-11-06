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
  public class RainSim : GH_Component {
    public RainSim()
        : base("RainSimulation", "batRainSim", "Waterflow simulation on the target terrain.", "BAT",
               "Analysis") {}

    public override GH_Exposure Exposure => GH_Exposure.primary;
    public override Guid ComponentGuid => new Guid("143c0ae4-319b-4ac7-9983-d3890d63870a");
    protected override System.Drawing.Bitmap Icon => Properties.Resources.rainSim;

    protected override void RegisterInputParams(GH_InputParamManager pManager) {
      pManager.AddMeshParameter("Mesh", "M", "Terrain to run the rain simulation.",
                                GH_ParamAccess.item);
      pManager.AddPointParameter(
          "Start", "P", "Start points for the rain simulation on or above the target terrain.",
          GH_ParamAccess.list);
      pManager.AddIntegerParameter("Steps", "S",
                                   "Maximum number of steps to run for the simulation.",
                                   GH_ParamAccess.item, 100);
      pManager.AddNumberParameter("StepSize", "D", "Distance to move in each simulation step.",
                                  GH_ParamAccess.item, 1.0);

      pManager[2].Optional = true;
      pManager[3].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager) {
      pManager.AddCurveParameter(
          "Flowline", "C",
          "Polyline curves that represent the simulation results of the water flow.",
          GH_ParamAccess.list);
    }

    protected override void SolveInstance(IGH_DataAccess DA) {
      Mesh terrain = null;
      List<Point3d> startPoints = new List<Point3d>();
      int maxSteps = 100;
      double stepSize = 1.0;

      // Get input data
      if (!DA.GetData("Mesh", ref terrain))
        return;
      if (!DA.GetDataList("Start", startPoints))
        return;
      if (!DA.GetData("Steps", ref maxSteps))
        return;
      if (!DA.GetData("StepSize", ref stepSize))
        return;

      // Validate inputs
      if (terrain == null) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input terrain mesh cannot be null");
        return;
      }

      if (startPoints == null || startPoints.Count == 0) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "At least one start point is required");
        return;
      }

      if (maxSteps <= 0) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Steps must be positive");
        maxSteps = Math.Max(1, maxSteps);
      }

      if (stepSize <= 0) {
        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "StepSize must be positive");
        stepSize = Math.Max(0.1, stepSize);
      }

      // Run simulation for each start point
      List<Curve> flowlines = new List<Curve>();

      foreach (Point3d startPoint in startPoints) {
        List<Point3d> flowPath = SimulateWaterFlowGradient(terrain, startPoint, maxSteps, stepSize);

        // Create polyline from flow path if we have at least 2 points
        if (flowPath.Count >= 2) {
          Polyline polyline = new Polyline(flowPath);
          flowlines.Add(polyline.ToNurbsCurve());
        }
      }

      // Output flowlines
      DA.SetDataList("Flowline", flowlines);
    }

    /// <summary>
    /// Simulate water flow using gradient-based continuous direction
    /// This produces smooth, natural-looking flow lines
    /// </summary>
    private List<Point3d> SimulateWaterFlowGradient(Mesh terrain, Point3d startPoint, int maxSteps,
                                                    double stepSize) {
      List<Point3d> flowPath = new List<Point3d>();

      // Project start point onto terrain
      Point3d currentPoint = TerrainAnalysisUtils.ProjectPointOntoTerrain(terrain, startPoint);

      // If projection failed, return empty path
      if (double.IsNaN(currentPoint.Z)) {
        return flowPath;
      }

      // Add start point
      flowPath.Add(currentPoint);

      // Use half the step size for gradient sampling (for more accurate gradient calculation)
      double sampleDistance = stepSize * 0.5;

      // Run simulation
      for (int step = 0; step < maxSteps; step++) {
        // Calculate the steepest descent direction using terrain gradient
        Vector3d descentDirection =
            TerrainAnalysisUtils.GetSteepestDescentDirection(terrain, currentPoint, sampleDistance);

        // Stop if no valid descent direction (flat area, depression, or edge)
        if (descentDirection == Vector3d.Zero || descentDirection.Length < 0.001) {
          break;
        }

        // Move in the descent direction by stepSize
        Point3d nextPointXY = currentPoint + descentDirection * stepSize;

        // Get the actual height at the new XY position
        double nextHeight = TerrainAnalysisUtils.GetMeshHeightAtPoint(terrain, nextPointXY);

        // Stop if we moved outside the mesh bounds
        if (double.IsNaN(nextHeight)) {
          break;
        }

        // Create the next point with proper height
        Point3d nextPoint = new Point3d(nextPointXY.X, nextPointXY.Y, nextHeight);

        // Optional: Check if we're actually going downhill (numerical stability)
        // Stop if height increased (shouldn't happen with gradient descent, but terrain sampling
        // can be noisy)
        if (nextPoint.Z > currentPoint.Z + 0.001) {
          break;
        }

        // Add point to path and continue
        flowPath.Add(nextPoint);
        currentPoint = nextPoint;
      }

      return flowPath;
    }
  }

  public class SlopeAnalysis : GH_Component {
    // Store the selected color map name
    private string selectedColorMap = "Viridis";

    public SlopeAnalysis()
        : base("SlopeAnalysis", "batSlope",
               "Analyze terrain slope and visualize with color-coded mesh.", "BAT", "Analysis") {}

    public override GH_Exposure Exposure => GH_Exposure.primary;
    public override Guid ComponentGuid => new Guid("2f4d8bc6-420c-4bf8-a094-e18a1e74981b");
    protected override System.Drawing.Bitmap Icon => Properties.Resources.slopeAnalysis;

    protected override void RegisterInputParams(GH_InputParamManager pManager) {
      pManager.AddMeshParameter("Mesh", "M", "Terrain mesh to analyze.", GH_ParamAccess.item);
      pManager.AddNumberParameter(
          "GridSize", "G",
          "Grid size for sampling terrain. Smaller values give more detail but slower computation.",
          GH_ParamAccess.item, 10.0);
      pManager.AddIntervalParameter(
          "SlopeBounds", "B",
          "Slope range for color mapping (degrees). Slopes outside this range are clamped to min/max colors. " +
              "Use this to focus on specific slope ranges.",
          GH_ParamAccess.item, new Interval(0, 90));

      pManager[1].Optional = true;
      pManager[2].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager) {
      pManager.AddMeshParameter(
          "AnalysisMesh", "M",
          "Color-coded mesh showing slope analysis. Colors mapped to SlopeBounds range.",
          GH_ParamAccess.item);
      pManager.AddPointParameter("SamplePoints", "P", "Grid points where slope was sampled.",
                                 GH_ParamAccess.list);
      pManager.AddNumberParameter("SlopeValues", "S", "Slope angle at each sample point (degrees).",
                                  GH_ParamAccess.list);
      pManager.AddVectorParameter("SlopeVectors", "V",
                                  "Direction vectors showing downhill direction at each point.",
                                  GH_ParamAccess.list);
      pManager.AddBooleanParameter("InBounds", "IB",
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
        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                          "SlopeBounds minimum must be less than maximum");
        slopeBounds = new Interval(0, 90);
      }

      // Get terrain bounding box
      BoundingBox bbox = terrain.GetBoundingBox(false);

      // Use half grid size for gradient sampling
      double sampleDistance = gridSize * 0.5;

      // Pre-compute gradient map for the entire terrain
      GradientMap gradientMap =
          TerrainAnalysisUtils.PrecomputeGradientMap(terrain, bbox, gridSize, sampleDistance);

      // Show warning if mesh normal fallback was used (extreme edge cases)
      if (gradientMap.MeshNormalFallbackCount > 0) {
        AddRuntimeMessage(
            GH_RuntimeMessageLevel.Warning,
            $"Mesh normal fallback used for {gradientMap.MeshNormalFallbackCount} points at mesh edges. " +
                "Consider using a smaller GridSize for better edge accuracy.");
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
          bool withinBounds =
              !double.IsNaN(slope) && slope >= slopeBounds.T0 && slope <= slopeBounds.T1;

          samplePoints.Add(point);
          slopeValues.Add(slope);
          slopeVectors.Add(downhillDir);
          inBounds.Add(withinBounds);
        }
      }

      // Create analysis mesh with color coding
      Mesh analysisMesh =
          CreateSlopeAnalysisMesh(gradientMap, bbox, gridSize, slopeBounds, selectedColorMap);

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
    private Mesh CreateSlopeAnalysisMesh(GradientMap gradientMap, BoundingBox bbox, double gridSize,
                                         Interval slopeBounds, string colorMapName) {
      Mesh analysisMesh = new Mesh();

      // Get the specified color map
      IColorMap colorMap = ColorMapHelper.GetColorMap(colorMapName);

      int xCount = gradientMap.XCount;
      int yCount = gradientMap.YCount;

      // Use the provided slope bounds for color mapping
      double minSlope = slopeBounds.T0;
      double maxSlope = slopeBounds.T1;

      // Create a mapping from grid indices to vertex indices (for handling NaN values)
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
              // Clamp slope to bounds - values outside range get clamped to min/max colors
              double clampedSlope = Math.Max(minSlope, Math.Min(maxSlope, slope));
              color = ColorMapHelper.MapValueToColor(clampedSlope, minSlope, maxSlope, colorMap);
            }

            analysisMesh.VertexColors.Add(color);
          }
        }
      }

      // Create triangular mesh faces (only for quads where all 4 vertices are valid)
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

    // Add right-click menu for color map selection (flat structure for cross-platform)
    public override void AppendAdditionalMenuItems(ToolStripDropDown menu) {
      base.AppendAdditionalMenuItems(menu);

      // Add a separator
      Menu_AppendSeparator(menu);

      // Add color map header (disabled)
      var headerItem = Menu_AppendItem(menu, "Color Map:");
      headerItem.Enabled = false;

      // Add each color map as a direct menu item
      foreach (var colorMap in ColorMapHelper.AvailableColorMaps) {
        Menu_AppendItem(menu, "  " + colorMap.Name, OnColorMapSelected, true,
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
}
