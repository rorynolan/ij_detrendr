import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;

import java.util.zip.DataFormatException;

class MyImg extends ImagePlus {
  MyImg() {
    super();
  }
  MyImg(ImagePlus imPlus) {
    super("new_MyImg", imPlus.getStack());
    int[] dim = imPlus.getDimensions();
    if (dim[3] > 1 && dim[4] > 1) {
      Exception e = new DataFormatException(
              "Your image has volume (z) and time (frame) axes, " +
                      "i.e. it has both slices and frames. " +
                      "detrendr is made to deal with one or the other, " +
                      "but it can't deal with both.");
      IJ.handleException(e);
    }
    if (dim[3] == 1 && dim[4] == 1) {
      Exception e = new DataFormatException(
              "Your image only has one frame. To be detrendable, " +
                      "an image must have more than one frame."
      );
      IJ.handleException(e);
    }
    ImageStack imStack = imPlus.getStack();
    if (imStack.isRGB() && dim[2] == 1) {
      ImageStack newStack = new ImageStack(width, height, dim[2] * dim[3]);
      for (int i = 1; i <= dim[3]; ++i) {
        ImageProcessor iProc = newStack.getProcessor(i);
      }
    }
    setDimensions(dim[2], Math.max(dim[3], dim[4]), Math.min(dim[3], dim[4]));
  }
}
