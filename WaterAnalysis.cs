using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
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
}
