package net.rorynolan.detrendr;

import Jama.Matrix;
import ij.ImagePlus;
import ij.process.AutoThresholder;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import java.util.zip.DataFormatException;

/**
 * Class containing static function for tresholding stacks.
 */
class Thresh {
  /**
   * Threshold a stack based on the mean intensity profile. For a given xy pixel position in the stack, either all of
   * the pixels in all of the frames at that position are thresholded away to zero or none are.
   * @param oneChImPlus A one channel, multi-frame ImagePlus.
   * @param method For choosing the threshold automatically. Must be one of "Default", "Huang", "Intermodes", "IsoData",
   *               "Li", "MaxEntropy", "Mean", "MinError","Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy",
   *               "Shanbhag", "Triangle", or "Yen".
   * @param manualThresh A positive number. Pixels with a mean intensity of less than manualThresh will be thresholded
   *                     away to zero. If manualThresh is set (greater than zero), it overrides method, i.e. it does not
   *                     matter what method is set to. To be clear, if you want to use manualThresh, set
   *                     method = "None".
   * @return The stack-thresholded image.
   * @throws DataFormatException if oneChImPlus has more than one channel or only one frame.
   */
  static ImagePlus stackThresh(ImagePlus oneChImPlus,
                                       String method,
                                       double manualThresh)
          throws DataFormatException {
    oneChImPlus = MyImg.assertOneChManyFrames(oneChImPlus, "stackThresh");
    int width = oneChImPlus.getDimensions()[0];
    Matrix mat = MyImg.convertToMatrix(oneChImPlus);
    int nRow = mat.getRowDimension();
    if (manualThresh > 0) {
      double[] colMeans = MyStats.colsMeans(mat);
      for (int i = 0; i != colMeans.length; ++i) {
        if (colMeans[i] < manualThresh) {
          for (int row = 0; row != nRow; ++row) {
            mat.set(row, i, 0);
          }
        }
      }
      return MyImg.convertToImagePlus(mat, width);
    }
    double[] colSums = MyStats.colsSums(mat);
    int nCol = colSums.length;
    SummaryStatistics sumStat = new SummaryStatistics();
    for (int i = 0; i != nCol; ++i) {
      sumStat.addValue(colSums[i]);
    }
    int min = (int) sumStat.getMin();
    int max = (int) sumStat.getMax();
    int prelimHistLength = max - min + 1;
    int histLength = (int) Math.pow(2, 8);
    int[] hist = new int[histLength];
    double adjustFactor = 1;
    if (prelimHistLength > histLength) {
      adjustFactor = (double) histLength / prelimHistLength;
      for (int col = 0; col != nCol; ++col) {
        int index = (int) (colSums[col] * adjustFactor);
        ++hist[index - min];
      }
    } else {
      for (int col = 0; col != nCol; ++col) {
        ++hist[(int) colSums[col] - min];
      }
    }
    AutoThresholder autoT = new AutoThresholder();
    double thresh = -1;
    switch (method) {
      case "Default":
        thresh = autoT.getThreshold("Default", hist);
        break;
      case "Huang":
        thresh = autoT.getThreshold("Huang", hist);
        break;
      case "Intermodes":
        thresh = autoT.getThreshold("Intermodes", hist);
        break;
      case "IsoData":
        thresh = autoT.getThreshold("IJ_IsoData", hist);
        break;
      case "Li":
        thresh = autoT.getThreshold("Li", hist);
        break;
      case "MaxEntropy":
        thresh = autoT.getThreshold("MaxEntropy", hist);
        break;
      case "Mean":
        thresh = autoT.getThreshold("Mean", hist);
        break;
      case "MinError":
        thresh = autoT.getThreshold("MinError", hist);
        break;
      case "Minimum":
        thresh = autoT.getThreshold("Minimum", hist);
        break;
      case "Moments":
        thresh = autoT.getThreshold("Moments", hist);
        break;
      case "Otsu":
        thresh = autoT.getThreshold("Otsu", hist);
        break;
      case "Percentile":
        thresh = autoT.getThreshold("Percentile", hist);
        break;
      case "RenyiEntropy":
        thresh = autoT.getThreshold("RenyiEntropy", hist);
        break;
      case "Shanbhag":
        thresh = autoT.getThreshold("Shanbhag", hist);
        break;
      case "Triangle":
        thresh = autoT.getThreshold("Triangle", hist);
        break;
      case "Yen":
        thresh = autoT.getThreshold("Yen", hist);
        break;
    }
    if (thresh < 0) {
      throw new IllegalArgumentException(
              "Threshold was calculated to be less than zero. " +
                      "This is an unexpected error and should never happen.");
    }
    if (adjustFactor != 1) {
      thresh /= adjustFactor;
    }
    thresh += min;
    for (int i = 0; i != colSums.length; ++i) {
      if (colSums[i] < thresh) {
        for (int row = 0; row != nRow; ++row) {
          mat.set(row, i, 0);
        }
      }
    }
    return MyImg.convertToImagePlus(mat, width);
  }
}

