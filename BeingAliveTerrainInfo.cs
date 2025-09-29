using System;
using System.Drawing;
using Grasshopper;
using Grasshopper.Kernel;

namespace BeingAliveTerrain {
  public class BeingAliveTerrainInfo : GH_AssemblyInfo {
    public override string Name => "BAT";

    // Return a 24x24 pixel bitmap to represent this GHA library.
    public override Bitmap Icon => null;

    // Return a short string describing the purpose of this GHA library.
    public override string Description => "";

    public override Guid Id => new Guid("d3b8edbb-54d2-4d80-ac59-a04b9f77f1dd");

    // Return a string identifying you or your company.
    public override string AuthorName => "Dr.Zhao.MA";

    // Return a string representing your preferred contact details.
    public override string AuthorContact => "zhma@ethz.ch";

    // Return a string representing the version.  This returns the same version as the assembly.
    public override string AssemblyVersion => GetType().Assembly.GetName().Version.ToString();
  }
}
