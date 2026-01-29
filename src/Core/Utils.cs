using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BeingAliveTerrain {
  internal class SysUtils {
    public static System.Drawing.Bitmap cvtByteBitmap(byte[] byteBitmap) {
      using (var stream = new System.IO.MemoryStream(byteBitmap)) {
        return new System.Drawing.Bitmap(stream: stream);
      }
    }
  }
}
