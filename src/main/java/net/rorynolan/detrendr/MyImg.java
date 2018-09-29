package net.rorynolan.detrendr;

import Jama.Matrix;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.ChannelSplitter;
import ij.plugin.CompositeConverter;
import ij.process.ImageProcessor;

import java.util.zip.DataFormatException;

class MyImg {

  static ImagePlus makeMine(ImagePlus imPlus)
          throws DataFormatException {
    int[] dim = imPlus.getDimensions();
    if (dim[3] > 1 && dim[4] > 1) {
      DataFormatException e = new DataFormatException(
              "Your image has volume (z) and time (frame) axes, " +
                      "i.e. it has both slices and frames. " +
                      "detrendr is made to deal with one or the other, " +
                      "but it can't deal with both.");
      throw e;
    }
    if (dim[3] == 1 && dim[4] == 1) {
      DataFormatException e = new DataFormatException(
              "Your image only has one frame. To be detrendable, " +
                      "an image must have more than one frame."
      );
      throw e;
    }
    ImagePlus out;
    if (imPlus.getType() == ImagePlus.COLOR_RGB) {
      out = CompositeConverter.makeComposite(imPlus);
      dim = out.getDimensions();
    } else {
      out = imPlus;
    }
    out.setDimensions(dim[2], Math.max(dim[3], dim[4]), Math.min(dim[3], dim[4]));
    return out;
  }

  static ImagePlus makeMyOneCh(ImagePlus imPlus,
                               int channel)
          throws IllegalArgumentException {  // channel is 1-based
    if (channel < 1) {
      IllegalArgumentException e = new IllegalArgumentException(
              "\n  * The channel number must be greater than or equal to 1."
      );
      throw e;
    }
    ImagePlus[] channelArr = ChannelSplitter.split(imPlus);
    int nCh = imPlus.getDimensions()[2];
    if (channel > nCh) {
      IllegalArgumentException e = new IllegalArgumentException(
              "\n  * You have requested channel " + channel +
                      " of the image, but the image has only " +
                      nCh + "channels in total."
      );
      throw e;
    }
    return channelArr[channel - 1];
  }

  static Matrix convertToMatrix(ImagePlus oneChImgPlus) {
    int[] imDim = oneChImgPlus.getDimensions();
    Matrix out = new Matrix(imDim[3], imDim[0] * imDim[1]);
    for (int slice = 0; slice != imDim[3]; ++slice) {
      oneChImgPlus.setSlice(slice);
      ImageProcessor ip = oneChImgPlus.getProcessor();
      for (int x = 0; x != imDim[0]; ++x) {
        for (int y = 0; y != imDim[1]; ++y) {
          out.set(slice, x + y * imDim[0], ip.get(x, y));
        }
      }
    }
    return out;
  }

  static ImagePlus convertToImagePlus(Matrix mat, int width)
          throws IllegalArgumentException {
    int matNCol = mat.getColumnDimension(), matNRow = mat.getRowDimension();
    if ((matNCol % width) != 0) {
      IllegalArgumentException e = new IllegalArgumentException(
              "width must divide evenly into the number of columns in mat. \n" +
                      "  * Your mat has " + matNCol + " columns and your width is " +
                      width + "."
      );
      throw e;
    }
    int height = matNCol / width;
    ImageStack imStack = ImageStack.create(width, height, matNRow, 16);
    for (int slice = 0; slice != matNRow; ++slice) {
      for (int x = 0; x != width; ++x) {
        for (int y = 0; y != height; ++y) {
          imStack.setVoxel(x, y, slice, mat.get(slice, x + y * width));
        }
      }
    }
    return new ImagePlus("created_from_jama_matrix", imStack);
  }

}
