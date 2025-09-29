using System;
using System.Drawing;

namespace BeingAliveTerrain {
  /// <summary>
  /// Base interface for color mapping functions
  /// </summary>
  public interface IColorMap {
    /// <summary>
    /// Maps a normalized value [0,1] to a color
    /// </summary>
    /// <param name="value">Normalized value between 0 and 1</param>
    /// <returns>Color corresponding to the value</returns>
    Color GetColor(double value);

    /// <summary>
    /// Name of the color map
    /// </summary>
    string Name { get; }
  }

  /// <summary>
  /// Parula colormap implementation (MATLAB-inspired)
  /// </summary>
  public class ParulaColorMap : IColorMap {
    public string Name => "Parula";

    public Color GetColor(double value) {
      // Clamp value to [0, 1]
      value = Math.Max(0, Math.Min(1, value));

      double r, g, b;

      if (value < 0.25) {
        double t = value / 0.25;
        r = 0.2078 + t * (0.1647 - 0.2078);
        g = 0.1663 + t * (0.5213 - 0.1663);
        b = 0.5292 + t * (0.9786 - 0.5292);
      } else if (value < 0.5) {
        double t = (value - 0.25) / 0.25;
        r = 0.1647 + t * (0.2365 - 0.1647);
        g = 0.5213 + t * (0.7720 - 0.5213);
        b = 0.9786 + t * (0.7741 - 0.9786);
      } else if (value < 0.75) {
        double t = (value - 0.5) / 0.25;
        r = 0.2365 + t * (0.7854 - 0.2365);
        g = 0.7720 + t * (0.9424 - 0.7720);
        b = 0.7741 + t * (0.3602 - 0.7741);
      } else {
        double t = (value - 0.75) / 0.25;
        r = 0.7854 + t * (0.9763 - 0.7854);
        g = 0.9424 + t * (0.9831 - 0.9424);
        b = 0.3602 + t * (0.0538 - 0.3602);
      }

      return Color.FromArgb((int)(r * 255), (int)(g * 255), (int)(b * 255));
    }
  }

  /// <summary>
  /// Jet colormap implementation (rainbow colors)
  /// </summary>
  public class JetColorMap : IColorMap {
    public string Name => "Jet";

    public Color GetColor(double value) {
      value = Math.Max(0, Math.Min(1, value));

      double r, g, b;

      if (value < 0.125) {
        r = 0;
        g = 0;
        b = 0.5 + value * 4;
      } else if (value < 0.375) {
        r = 0;
        g = (value - 0.125) * 4;
        b = 1;
      } else if (value < 0.625) {
        r = (value - 0.375) * 4;
        g = 1;
        b = 1 - (value - 0.375) * 4;
      } else if (value < 0.875) {
        r = 1;
        g = 1 - (value - 0.625) * 4;
        b = 0;
      } else {
        r = 1 - (value - 0.875) * 4;
        g = 0;
        b = 0;
      }

      return Color.FromArgb((int)(r * 255), (int)(g * 255), (int)(b * 255));
    }
  }

  /// <summary>
  /// Plasma colormap implementation (perceptually uniform)
  /// </summary>
  public class PlasmaColorMap : IColorMap {
    public string Name => "Plasma";

    public Color GetColor(double value) {
      value = Math.Max(0, Math.Min(1, value));

      double r, g, b;

      if (value < 0.25) {
        double t = value / 0.25;
        r = 0.050 + t * (0.308 - 0.050);
        g = 0.030 + t * (0.061 - 0.030);
        b = 0.528 + t * (0.727 - 0.528);
      } else if (value < 0.5) {
        double t = (value - 0.25) / 0.25;
        r = 0.308 + t * (0.624 - 0.308);
        g = 0.061 + t * (0.185 - 0.061);
        b = 0.727 + t * (0.616 - 0.727);
      } else if (value < 0.75) {
        double t = (value - 0.5) / 0.25;
        r = 0.624 + t * (0.897 - 0.624);
        g = 0.185 + t * (0.465 - 0.185);
        b = 0.616 + t * (0.281 - 0.616);
      } else {
        double t = (value - 0.75) / 0.25;
        r = 0.897 + t * (0.940 - 0.897);
        g = 0.465 + t * (0.797 - 0.465);
        b = 0.281 + t * (0.165 - 0.281);
      }

      return Color.FromArgb((int)(r * 255), (int)(g * 255), (int)(b * 255));
    }
  }

  /// <summary>
  /// Viridis colormap implementation (perceptually uniform, colorblind-friendly)
  /// </summary>
  public class ViridisColorMap : IColorMap {
    public string Name => "Viridis";

    public Color GetColor(double value) {
      value = Math.Max(0, Math.Min(1, value));

      double r, g, b;

      if (value < 0.25) {
        double t = value / 0.25;
        r = 0.267 + t * (0.229 - 0.267);
        g = 0.005 + t * (0.322 - 0.005);
        b = 0.329 + t * (0.545 - 0.329);
      } else if (value < 0.5) {
        double t = (value - 0.25) / 0.25;
        r = 0.229 + t * (0.127 - 0.229);
        g = 0.322 + t * (0.566 - 0.322);
        b = 0.545 + t * (0.550 - 0.545);
      } else if (value < 0.75) {
        double t = (value - 0.5) / 0.25;
        r = 0.127 + t * (0.369 - 0.127);
        g = 0.566 + t * (0.788 - 0.566);
        b = 0.550 + t * (0.382 - 0.550);
      } else {
        double t = (value - 0.75) / 0.25;
        r = 0.369 + t * (0.993 - 0.369);
        g = 0.788 + t * (0.906 - 0.788);
        b = 0.382 + t * (0.144 - 0.382);
      }

      return Color.FromArgb((int)(r * 255), (int)(g * 255), (int)(b * 255));
    }
  }

  /// <summary>
  /// Red-White-Blue colormap implementation (white center for cut/fill analysis)
  /// </summary>
  public class RedWhiteBlueColorMap : IColorMap {
    public string Name => "RedWhiteBlue";

    public Color GetColor(double value) {
      value = Math.Max(0, Math.Min(1, value));

      double r, g, b;

      if (value < 0.5) {
        // Interpolate from red (0) to white (0.5)
        double t = value * 2; // Scale to [0,1]
        r = 1.0;
        g = t;
        b = t;
      } else {
        // Interpolate from white (0.5) to blue (1)
        double t = (value - 0.5) * 2; // Scale to [0,1]
        r = 1.0 - t;
        g = 1.0 - t;
        b = 1.0;
      }

      return Color.FromArgb((int)(r * 255), (int)(g * 255), (int)(b * 255));
    }
  }

  /// <summary>
  /// Red-White-Green colormap implementation (white center with intensity fading)
  /// </summary>
  public class RedWhiteGreenColorMap : IColorMap {
    public string Name => "RedWhiteGreen";

    public Color GetColor(double value) {
      value = Math.Max(0, Math.Min(1, value));

      double r, g, b;

      if (value < 0.5) {
        // Interpolate from red (0) to white (0.5)
        // As value approaches 0.5, colors fade to white
        double t = value * 2; // Scale to [0,1]
        double intensity = 1.0 - t; // Higher intensity when further from center
        
        r = 1.0; // Always red component
        g = t; // Fade to white
        b = t; // Fade to white
      } else {
        // Interpolate from white (0.5) to green (1)
        // As value moves away from 0.5, colors get more intense
        double t = (value - 0.5) * 2; // Scale to [0,1]
        double intensity = t; // Higher intensity when further from center
        
        r = 1.0 - t; // Fade from white
        g = 1.0; // Always green component
        b = 1.0 - t; // Fade from white
      }

      return Color.FromArgb((int)(r * 255), (int)(g * 255), (int)(b * 255));
    }
  }

  /// <summary>
  /// Utility class for color mapping operations
  /// </summary>
  public static class ColorMapHelper {
    /// <summary>
    /// Available color maps
    /// </summary>
    public static readonly IColorMap[] AvailableColorMaps = {
      new ParulaColorMap(), new JetColorMap(), new PlasmaColorMap(), new ViridisColorMap(), new RedWhiteBlueColorMap(), new RedWhiteGreenColorMap()
    };

    /// <summary>
    /// Get a color map by name
    /// </summary>
    /// <param name="name">Name of the color map</param>
    /// <returns>Color map instance or Parula if not found</returns>
    public static IColorMap GetColorMap(string name) {
      foreach (var colorMap in AvailableColorMaps) {
        if (colorMap.Name.Equals(name, StringComparison.OrdinalIgnoreCase)) {
          return colorMap;
        }
      }
      return new ParulaColorMap();  // Default fallback
    }

    /// <summary>
    /// Maps a value from input range to normalized [0,1] range and applies color map
    /// </summary>
    /// <param name="value">Input value</param>
    /// <param name="minValue">Minimum of input range</param>
    /// <param name="maxValue">Maximum of input range</param>
    /// <param name="colorMap">Color map to use</param>
    /// <returns>Mapped color</returns>
    public static Color MapValueToColor(double value, double minValue, double maxValue,
                                        IColorMap colorMap) {
      if (maxValue <= minValue) {
        return colorMap.GetColor(0.5);  // Neutral color for invalid range
      }

      double normalized = (value - minValue) / (maxValue - minValue);
      return colorMap.GetColor(normalized);
    }
  }
}
