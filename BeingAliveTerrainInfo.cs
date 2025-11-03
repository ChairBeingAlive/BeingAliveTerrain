using System;
using System.Drawing;
using Grasshopper;
using Grasshopper.Kernel;

namespace BeingAliveTerrain {
  public class BeingAliveTerrainInfo : GH_AssemblyInfo {
    public override string Name => "BeingAliveTerrain";

    // Return a 24x24 pixel bitmap to represent this GHA library.
    public override Bitmap Icon => Properties.Resources.appIcon;

    // Return a short string describing the purpose of this GHA library.
    public override string Description =>
        "Terrain Plugin developed at the Chair of Being Alive @ ETH Zurich.";

    public override Guid Id => new Guid("d3b8edbb-54d2-4d80-ac59-a04b9f77f1dd");

    // Return a string identifying you or your company.
    public override string AuthorName => "Dr.Zhao.MA";

    // Return a string representing your preferred contact details.
    public override string AuthorContact => "zhma@ethz.ch";

    // public override string AssemblyVersion => GetType().Assembly.GetName().Version.ToString();
    public override string AssemblyVersion => "0.1.3";

    // this is currently the variable used by McNeel for plugin system
    public override string Version => AssemblyVersion;
  }
}
