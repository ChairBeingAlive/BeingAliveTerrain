using System;
using System.Drawing;
using Grasshopper;
using Grasshopper.Kernel;

namespace BeingAliveTerrain {
public class batInfo : GH_AssemblyInfo {
  public override string Name => "BeingAliveTerrain";

  // Return a 24x24 pixel bitmap to represent this GHA library.
  public override Bitmap Icon => Properties.Resources.appIcon;

  // Return a short string describing the purpose of this GHA library.
  public override string Description =>
      "A mesh-based terrain modelling and analysis plugin developed at the Chair of Being Alive @ ETH Zurich.";

  public override Guid Id => new Guid("d3b8edbb-54d2-4d80-ac59-a04b9f77f1dd");

  // Return a string identifying you or your company.
  public override string AuthorName => "Dr.Zhao.MA";

  // Return a string representing your preferred contact details.
  public override string AuthorContact => "https://github.com/ChairBeingAlive/BeingAliveTerrain";

  // public override string AssemblyVersion => GetType().Assembly.GetName().Version.ToString();
  public override string AssemblyVersion => "0.2.1";

  // this is currently the variable used by McNeel for plugin system
  public override string Version => AssemblyVersion;
}

// update plugin icons and category in the tab
public class BAT_CategoryIcon : GH_AssemblyPriority {
  public override GH_LoadingInstruction PriorityLoad() {
    Grasshopper.Instances.ComponentServer.AddCategoryIcon("BAT", Properties.Resources.appIcon);
    Grasshopper.Instances.ComponentServer.AddCategorySymbolName("BAT", 'B');
    return GH_LoadingInstruction.Proceed;
  }
}

}
