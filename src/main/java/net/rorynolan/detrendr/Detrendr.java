package net.rorynolan.detrendr;

import Jama.Matrix;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.plugin.ChannelSplitter;
import ij.plugin.CompositeConverter;
import ij.plugin.RGBStackMerge;
import ij.process.AutoThresholder;
import ij.process.ImageProcessor;
import net.imagej.ImageJ;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import java.util.NoSuchElementException;
import java.util.stream.IntStream;
import java.util.zip.DataFormatException;


/**
 * The Detrendr ImageJ plugin class.
 */
@Plugin(type = Command.class, menuPath = "Plugins>Detrendr")
public class Detrendr implements Command {

  /**
   * This is just for the developers of the package to test how it's going.
   *
   * @param args Should be left blank.
   * @throws NoSuchElementException
   */
  public static void main(String[] args) throws NoSuchElementException {
    Class<?> cls = Detrendr.class;
    String url = cls.getResource("/" + cls.getName().replace('.', '/') + ".class").toString();
    int lastIdx = url.lastIndexOf('/');
    String pluginsDir = url.substring(5, lastIdx);
    System.out.println("Plugins directory: " + pluginsDir);
    System.setProperty("plugins.dir", pluginsDir);
    final ImageJ ij = net.imagej.Main.launch(args);
  }

  @Parameter
  private LogService LogService;

  @Parameter(
          label = "Automatic thresholding method",
          choices = {"None", "Default", "Huang", "Intermodes",
                  "IsoData", "Li", "MaxEntropy", "Mean", "MinError",
                  "Minimum", "Moments", "Otsu", "Percentile",
                  "RenyiEntropy", "Shanbhag", "Triangle", "Yen"},
          description = "Select 'None' to forego thresholding. " +
                  "To apply thresholding to your image as a preprocessing step, " +
                  "select one of the methods. Browse these methods at " +
                  "https://imagej.net/Auto_Threshold."
  )
  private String autoThreshMethod;

  @Parameter(
          label = "Manual Threshold", initializer = "myZero", min = "0",
          description = "Manually threshold your image. " +
                  "If this value is greater than zero, it overrides " +
                  "automatic thresholding.",
          persist = false
  )
  private double manualThresh;

  @Parameter(
          label = "Seed",
          description = "Seed for random number generation " +
                  "(if you don't know what this is for, " +
                  "leave at 1, 1 is fine).",
          initializer = "myOne",
          min = "1", max = "9999"
  )
  private int seed;

  private double myZero() {
    return 0;
  }

  private int myOne() {
    return 1;
  }

  /**
   * This is the method that is actually called when the plugin is opened in ImageJ.
   */
  @Override
  public void run() {
    int nOpenImgs = WindowManager.getWindowCount();
    if (nOpenImgs == 0) {
      Exception e = new NoSuchElementException(
              "\n  * The Detrendr plugin detrends " +
                      "the active image. \n" +
                      "    - You have no open images, " +
                      "so the operation has failed."
      );
      LogService.error(e);
      return;
    }
    ImagePlus img = WindowManager.getCurrentImage();
    String title = img.getTitle();
    try {
      ImagePlus preProcessed;
      ImagePlus out;
      if (manualThresh > 0 || (!autoThreshMethod.equals("None"))) {
        preProcessed = stackThresh(img, autoThreshMethod, manualThresh);
        out = detrend(preProcessed, seed);
      } else {
        out = detrend(img.duplicate(), seed);
      }
      out.setTitle(insertIntoFileName(title, "detrended"));
      out.show();
    } catch (DataFormatException e) {
      LogService.error(e);
    }
  }

  // Image Manipulation -------------------------------------------------------------------------

  private ImagePlus makeMine(ImagePlus imPlus)
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

  private ImagePlus makeMyOneCh(ImagePlus imPlus,
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

  ImagePlus assertOneChManyFrames(ImagePlus imPlus, String currentFun)
          throws DataFormatException {
    int nCh = imPlus.getDimensions()[2];
    if (nCh != 1) {
      throw new IllegalArgumentException("The function " + currentFun +
              "() expects an ImagePlus with one channel.\n" +
              "  * You have passed an ImagePlus with " + nCh + "channels.");
    }
    return makeMine(imPlus);
  }

  Matrix convertToMatrix(ImagePlus oneChImgPlus)
          throws DataFormatException {
    oneChImgPlus = assertOneChManyFrames(oneChImgPlus, "convertToMatrix");
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

  ImagePlus convertToImagePlus(Matrix mat, int width)
          throws IllegalArgumentException {
    int matNCol = mat.getColumnDimension(), matNRow = mat.getRowDimension();
    if (width < 1) {
      throw new IllegalArgumentException(
              "width must be positive.\n" +
                      "  * You have specified a width of " + width + "."
      );
    }
    if ((matNCol % width) != 0) {
      throw new IllegalArgumentException(
              "width must divide evenly into the number of columns in mat. \n" +
                      "  * Your mat has " + matNCol + " columns and your width is " +
                      width + "."
      );
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

  // --------------------------------------------------------------------------------------------

  // Matrix Manipulation ------------------------------------------------------------------------

  private double[] getRow(Matrix mat, int row) {
    // include a check on value of row if get time
    int nCol = mat.getColumnDimension();
    double[] out = new double[nCol];
    for (int col = 0; col != nCol; ++col) {
      out[col] = mat.get(row, col);
    }
    return out;
  }

  private double[] getCol(Matrix mat, int col) {
    // include a check on value of col if get time
    int nRow = mat.getRowDimension();
    double[] out = new double[nRow];
    for (int row = 0; row != nRow; ++row) {
      out[row] = mat.get(row, col);
    }
    return out;
  }

  // --------------------------------------------------------------------------------------------

  // String Manipulation

  private String last(String[] s) {
    int l = s.length;
    return s[l - 1];
  }

  private String insertIntoFileName(String fileName, String s) {
    int l = fileName.length();
    if (fileName.equals(".")) {
      return s;
    }
    if (fileName.matches(".*\\..*")) {
      if (l > 1 &&
              fileName.substring(0, 1).equals(".") &&
              (!fileName.substring(1, l).matches("\\."))) {
        return (fileName + "_" + s);
      }
      if (fileName.substring(l - 1, l).equals(".")) {
        return fileName.substring(0, l - 1) + "_" + s;
      }
      String[] dotSplit = fileName.split("\\.");
      StringBuilder out = new StringBuilder();
      int dsl = dotSplit.length;
      for (int i = 0; i != dsl - 1; ++i) {
        out.append(dotSplit[i]);
        if (i != dsl - 2) {
          out.append(".");
        }
      }
      out.append("_").append(s).append(".").append(last(dotSplit));
      return out.toString();
    }
    return fileName + s;
  }

  // --------------------------------------------------------------------------------------------

  // Summary Statistics -------------------------------------------------------------------------

  private double sum(double[] data) {
    double sum = 0.0;
    for (double a : data)
      sum += a;
    return sum;
  }

  double mean(double[] data) {
    return sum(data) / data.length;
  }

  double var(double[] data) {
    double mean = mean(data);
    double temp = 0;
    for (double a : data)
      temp += (a - mean) * (a - mean);
    return temp / (data.length - 1);
  }

  private double brightnessB(SummaryStatistics sumStat) {
    return sumStat.getVariance() / sumStat.getMean();
  }

  private double[] colsBrightnessB(Matrix mat) {
    int nCol = mat.getColumnDimension();
    int nRow = mat.getRowDimension();
    double[] out = new double[nCol];
    SummaryStatistics sumStat = new SummaryStatistics();
    for (int col = 0; col != nCol; ++col) {
      for (int row = 0; row != nRow; ++row) {
        sumStat.addValue(mat.get(row, col));
      }
      out[col] = brightnessB(sumStat);
      sumStat.clear();
    }
    return out;
  }

  double colsMeanBrightnessB(Matrix mat) {
    double[] cBB = colsBrightnessB(mat);
    SummaryStatistics sumStat = new SummaryStatistics();
    for (int col = 0; col != cBB.length; ++col) {
      if (!Double.isNaN(cBB[col])) {
        sumStat.addValue(cBB[col]);
      }
    }
    long sumStatN = sumStat.getN();
    if (sumStatN == 0) {
      return Double.NaN;
    } else {
      return sumStat.getMean();
    }
  }

  private double[] rowsSums(Matrix mat) {
    int nRow = mat.getRowDimension();
    int nCol = mat.getColumnDimension();
    double[] out = new double[nRow];
    SummaryStatistics sumStat = new SummaryStatistics();
    for (int row = 0; row != nRow; ++row) {
      for (int col = 0; col != nCol; ++col) {
        sumStat.addValue(mat.get(row, col));
      }
      out[row] = sumStat.getSum();
      sumStat.clear();
    }
    return out;
  }

  private double[] colsSums(Matrix mat) {
    int nRow = mat.getRowDimension();
    int nCol = mat.getColumnDimension();
    double[] out = new double[nCol];
    SummaryStatistics sumStat = new SummaryStatistics();
    for (int col = 0; col != nCol; ++col) {
      for (int row = 0; row != nRow; ++row) {
        sumStat.addValue(mat.get(row, col));
      }
      out[col] = sumStat.getSum();
      sumStat.clear();
    }
    return out;
  }

  double[] rowsMeans(Matrix mat) {
    int nRow = mat.getRowDimension();
    int nCol = mat.getColumnDimension();
    double[] out = rowsSums(mat);
    for (int row = 0; row != nRow; ++row) {
      out[row] /= nCol;
    }
    return out;
  }

  double[] colsMeans(Matrix mat) {
    int nRow = mat.getRowDimension();
    int nCol = mat.getColumnDimension();
    double[] out = colsSums(mat);
    for (int col = 0; col != nCol; ++col) {
      out[col] /= nRow;
    }
    return out;
  }

  private boolean allAreNum(double[] data, double num) {
    for (double d : data) {
      if (d != num) {
        return false;
      }
    }
    return true;
  }

  private double absDiff(double x, double y) {
    return Math.abs(x - y);
  }

  private int whichMin(double[] x) {
    int minIndex = 0;
    double min = x[minIndex];
    int l = x.length;
    for (int i = 0; i != l; ++i) {
      if (x[i] < min) {
        minIndex = i;
        min = x[minIndex];
      }
    }
    return minIndex;
  }

  // --------------------------------------------------------------------------------------------

  // Thresholding -------------------------------------------------------------------------------

  /**
   * Threshold a stack based on the mean intensity profile. For a given xy pixel position in the stack, either all of
   * the pixels in all of the frames at that position are thresholded away to zero or none are.
   *
   * @param oneChImPlus  A one channel, multi-frame ImagePlus.
   * @param method       For choosing the threshold automatically. Must be one of "Default", "Huang", "Intermodes", "IsoData",
   *                     "Li", "MaxEntropy", "Mean", "MinError","Minimum", "Moments", "Otsu", "Percentile", "RenyiEntropy",
   *                     "Shanbhag", "Triangle", or "Yen".
   * @param manualThresh A positive number. Pixels with a mean intensity of less than manualThresh will be thresholded
   *                     away to zero. If manualThresh is set (greater than zero), it overrides method, i.e. it does not
   *                     matter what method is set to. To be clear, if you want to use manualThresh, set
   *                     method = "None".
   * @return The stack-thresholded image.
   * @throws DataFormatException if oneChImPlus has more than one channel or only one frame.
   */
  ImagePlus stackThresh(ImagePlus oneChImPlus,
                               String method,
                               double manualThresh)
          throws DataFormatException {
    oneChImPlus = assertOneChManyFrames(oneChImPlus, "stackThresh");
    int width = oneChImPlus.getDimensions()[0];
    Matrix mat = convertToMatrix(oneChImPlus);
    int nRow = mat.getRowDimension();
    if (manualThresh > 0) {
      double[] colMeans = colsMeans(mat);
      for (int i = 0; i != colMeans.length; ++i) {
        if (colMeans[i] < manualThresh) {
          for (int row = 0; row != nRow; ++row) {
            mat.set(row, i, 0);
          }
        }
      }
      return convertToImagePlus(mat, width);
    }
    double[] colSums = colsSums(mat);
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
    return convertToImagePlus(mat, width);
  }

  // --------------------------------------------------------------------------------------------

  // Simulation ---------------------------------------------------------------------------------

  private Matrix poisMat(int nPx, double[] frameMeans,
                         int seed) {
    int nFrames = frameMeans.length;
    Matrix out = new Matrix(nFrames, nPx);
    for (int frame = 0; frame != nFrames; ++frame) {
      PoissonDistribution pd = new PoissonDistribution(frameMeans[frame]);
      pd.reseedRandomGenerator(seed++);
      for (int px = 0; px != nPx; ++px) {
        out.set(frame, px, pd.sample());
      }
    }
    return out;
  }

  private Matrix simMat(ImagePlus oneChImPlus, int seed)
          throws DataFormatException {
    oneChImPlus = assertOneChManyFrames(oneChImPlus, "simMat");
    Matrix oneChMat = convertToMatrix(oneChImPlus);
    int oneChMatNCol = oneChMat.getColumnDimension();
    int oneChMatNRow = oneChMat.getRowDimension();
    int nNonZeroCols = 0;
    for (int col = 0; col != oneChMatNCol; ++col) {
      if (!allAreNum(getCol(oneChMat, col), 0)) {
        ++nNonZeroCols;
      }
    }
    double[] frameMeans = rowsSums(oneChMat);
    for (int frame = 0; frame != oneChMatNRow; ++frame) {
      frameMeans[frame] /= nNonZeroCols;
    }
    return poisMat(nNonZeroCols, frameMeans, seed);
  }

  // --------------------------------------------------------------------------------------------

  // Detrending ---------------------------------------------------------------------------------

  private long getMaxSwaps(Matrix mat) {
    int nRow = mat.getRowDimension();
    double[] rowSums = rowsSums(mat);
    double meanRowSum = mean(rowSums);
    double[] framesCanLose = new double[nRow];
    double[] framesCanGet = framesCanLose.clone();
    for (int frame = 0; frame != nRow; ++frame) {
      framesCanLose[frame] = Math.max(0, rowSums[frame] - meanRowSum);
      framesCanGet[frame] = Math.max(0, meanRowSum - rowSums[frame]);
    }
    return (long) Math.min(sum(framesCanGet),
            sum(framesCanLose));
  }

  private long getMaxSwaps(double[] framesCanLose, double[] framesCanGet) {
    return (long) Math.min(sum(framesCanLose),
            sum(framesCanGet));
  }

  private void performSwaps(Matrix mat, Matrix origMat,
                            long nSwaps, long seed) {
    if (nSwaps == 0) {
      return;
    }
    int nRow = mat.getRowDimension(), nCol = mat.getColumnDimension();
    double[] rowSums = rowsSums(mat);
    double meanRowSum = mean(rowSums);
    double origMean = mean(origMat.getRowPackedCopy());
    double[] frameWeights = rowsMeans(origMat);
    for (int w = 0; w != nRow; ++w) {
      frameWeights[w] -= origMean;
    }
    double[] frameLosingWeights = frameWeights.clone();
    double[] frameGettingWeights = frameWeights.clone();
    for (int w = 0; w != nRow; ++w) {
      if (frameWeights[w] > 0) {
        frameGettingWeights[w] = 0;
      } else {
        frameLosingWeights[w] = 0;
        frameGettingWeights[w] *= -1;
      }
    }
    int[] frameIndex = IntStream.range(0, nRow).toArray();
    EnumeratedIntegerDistribution frameLosingIntDist =
            new EnumeratedIntegerDistribution(frameIndex, frameLosingWeights);
    frameLosingIntDist.reseedRandomGenerator(seed++);
    EnumeratedIntegerDistribution frameGettingIntDist =
            new EnumeratedIntegerDistribution(frameIndex, frameGettingWeights);
    frameGettingIntDist.reseedRandomGenerator(seed++);
    EnumeratedIntegerDistribution[] pxGiveWeights =
            new EnumeratedIntegerDistribution[nRow];
    int[] pxIndex = IntStream.range(0, nCol).toArray();
    Matrix weightMat = origMat.copy();
    for (int row = 0; row != nRow; ++row) {
      double[] weightRow = getRow(weightMat, row);
      double[] matRow = getRow(mat, row);
      for (int col = 0; col != nCol; ++col) {
        if (matRow[col] == 0) {
          weightRow[col] = 0;
        }
      }
      pxGiveWeights[row] = new EnumeratedIntegerDistribution(pxIndex, weightRow);
      pxGiveWeights[row].reseedRandomGenerator(seed++);
    }
    double[] framesCanLose = new double[nRow];
    double[] framesCanGet = framesCanLose.clone();
    for (int frame = 0; frame != nRow; ++frame) {
      framesCanLose[frame] = Math.max(0, rowSums[frame] - meanRowSum);
      framesCanGet[frame] = Math.max(0, meanRowSum - rowSums[frame]);
    }
    long maxSwaps = getMaxSwaps(framesCanLose, framesCanGet);
    long performingSwaps = Math.min(nSwaps, maxSwaps);
    for (long swapNum = 0; swapNum != performingSwaps; ++swapNum) {
      int frameLosing = frameLosingIntDist.sample();
      int pxLosing = pxGiveWeights[frameLosing].sample();
      mat.set(frameLosing, pxLosing, mat.get(frameLosing, pxLosing) - 1);
      if (mat.get(frameLosing, pxLosing) == 0) {
        weightMat.set(frameLosing, pxLosing, 0);
        double[] weightRow = getRow(weightMat, frameLosing);
        pxGiveWeights[frameLosing] = new EnumeratedIntegerDistribution(pxIndex, weightRow);
        pxGiveWeights[frameLosing].reseedRandomGenerator(seed++);
      }
      int frameGetting = frameGettingIntDist.sample();
      mat.set(frameGetting, pxLosing, mat.get(frameGetting, pxLosing) + 1);
      if (swapNum != performingSwaps - 1) {
        if ((--framesCanLose[frameLosing]) <= 0) {
          frameLosingWeights[frameLosing] = 0;
          frameLosingIntDist = new EnumeratedIntegerDistribution(frameIndex, frameLosingWeights);
        }
        if ((--framesCanGet[frameGetting]) <= 0) {
          frameGettingWeights[frameGetting] = 0;
          frameGettingIntDist = new EnumeratedIntegerDistribution(frameIndex, frameGettingWeights);
        }
      }
    }
  }

  /**
   * Perform the swaps for Robin Hood detrending.
   *
   * @param oneChImPlus A one-channel, multi-frame ImagePlus.
   * @param nSwaps      The number of swaps to perform, probably calculated with calcIdealSwaps.
   * @param seed        For random number generation, to select which counts to swap.
   * @return The detrended image.
   * @throws DataFormatException if oneChImPlus has more than one channel or only one frame.
   */
  ImagePlus performSwaps(ImagePlus oneChImPlus, long nSwaps, int seed)
          throws DataFormatException {
    oneChImPlus = assertOneChManyFrames(oneChImPlus, "performSwaps");
    if (nSwaps < 0) {
      throw new IllegalArgumentException(
              "nSwaps must be positive.\n" +
                      "  * You have specified nSwaps = " + nSwaps + "."
      );
    }
    int width = oneChImPlus.getDimensions()[0];
    Matrix mat = convertToMatrix(oneChImPlus);
    performSwaps(mat, mat.copy(), nSwaps, seed);
    return convertToImagePlus(mat, width);
  }

  /**
   * Calculate the ideal number of swaps to make on an image during Robin Hood detrending.
   *
   * @param oneChImPlus A one-channel, multi-frame ImagePlus.
   * @param seed        For random number generation during simulations.
   * @return The ideal number of swaps.
   * @throws DataFormatException if oneChImPlus has more than one channel or only one frame.
   */
  long calcIdealSwaps(ImagePlus oneChImPlus, int seed)
          throws DataFormatException {
    oneChImPlus = assertOneChManyFrames(oneChImPlus, "calcIdealSwaps");
    Matrix simMat = simMat(oneChImPlus, seed++);
    double lowerMeanB = colsMeanBrightnessB(simMat);
    if (lowerMeanB <= 1) {
      return 0;
    }
    Matrix leastMat = simMat.copy();
    long maxSwaps = getMaxSwaps(simMat);
    long mostSwaps = 1, leastSwaps = 0;
    Matrix mostMat = leastMat.copy();
    performSwaps(mostMat, simMat, mostSwaps, seed++);
    while (colsMeanBrightnessB(mostMat) > 1) {
      leastMat = mostMat.copy();
      leastSwaps = mostSwaps;
      long moreSwaps = Math.min(mostSwaps, maxSwaps - mostSwaps);
      performSwaps(mostMat, simMat, moreSwaps, seed++);
      mostSwaps += moreSwaps;
      if (mostSwaps == maxSwaps) {
        break;
      }
      if (LogService != null) {  // prevents log attempts in wrong context
        LogService.info("Trying " + mostSwaps + " swaps . . .");
      }
    }
    if (colsMeanBrightnessB(mostMat) > 1) {
      return maxSwaps;
    }
    long middleSwaps = (leastSwaps + mostSwaps) / 2;
    Matrix middleMat = leastMat.copy();
    performSwaps(middleMat, simMat, middleSwaps - leastSwaps, seed++);
    while ((mostSwaps - middleSwaps > 0) && (middleSwaps - leastSwaps > 0)) {
      double middleMatMeanB = colsMeanBrightnessB(middleMat);
      if (middleMatMeanB > 1) {
        mostSwaps = middleSwaps;
        mostMat = middleMat;
        middleSwaps = (leastSwaps + middleSwaps) / 2;
      } else if (middleMatMeanB < 1) {
        leastSwaps = middleSwaps;
        leastMat = middleMat;
        middleSwaps = (middleSwaps + mostSwaps) / 2;
      } else {
        return middleSwaps;  // this line will probably never be run
      }
      if (LogService != null) {  // prevents log attempts in wrong context
        LogService.info("Trying " + middleSwaps + " swaps . . .");
      }
      middleMat = leastMat.copy();
      performSwaps(middleMat, simMat, middleSwaps - leastSwaps, seed++);
    }
    long[] swaps = {leastSwaps, middleSwaps, mostSwaps};
    Matrix[] mats = {leastMat, middleMat, mostMat};
    double diffs1[] = new double[3];
    for (int i = 0; i != 3; ++i) {
      diffs1[i] = absDiff(colsMeanBrightnessB(mats[i]), 1);
    }
    long out = swaps[whichMin(diffs1)];
    if (LogService != null) {  // prevents log attempts in wrong context
      LogService.info("Settling on " + out + " swaps :-)");
    }
    return out;
  }

  /**
   * Detrend an ImagePlus with the Robin Hood algorithm, automatically choosing the number of swaps.
   *
   * @param imPlus The multi-frame image to detrend.
   * @param seed   For random number generation for simulations and swapping.
   * @return The detrended image.
   * @throws DataFormatException if imPlus has time and z axes.
   */
  ImagePlus detrend(ImagePlus imPlus, int seed)
          throws DataFormatException {
    ImagePlus myImPlus = makeMine(imPlus);
    int nCh = myImPlus.getDimensions()[2];
    // check for nCh <= 7
    ImagePlus[] channelArr = new ImagePlus[nCh];
    ImagePlus[] outChannelArr = new ImagePlus[nCh];
    for (int i = 0; i != nCh; ++i) {
      channelArr[i] = makeMyOneCh(myImPlus, i + 1);
      if (LogService != null) {  // prevents log attempts in wrong context
        LogService.info("Detrending channel " + (i + 1) + " . . .");
      }
      long nSwaps = calcIdealSwaps(channelArr[i], seed++);
      outChannelArr[i] = performSwaps(channelArr[i], nSwaps, seed++);
      if (LogService != null) {  // prevents log attempts in wrong context
        LogService.info("Finished detrending channel " + (i + 1) + ".\n");
      }
    }
    if (LogService != null) {  // prevents log attempts in wrong context
      LogService.info("Done.");
    }
    if (outChannelArr.length == 1) {
      return outChannelArr[0];
    }
    return RGBStackMerge.mergeChannels(outChannelArr, false);
  }

}