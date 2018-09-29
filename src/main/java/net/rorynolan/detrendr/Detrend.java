package net.rorynolan.detrendr;

import Jama.Matrix;
import ij.ImagePlus;
import ij.plugin.RGBStackMerge;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;

import java.util.stream.IntStream;
import java.util.zip.DataFormatException;

class Detrend {

  private static long getMaxSwaps(Matrix mat) {
    int nRow = mat.getRowDimension();
    double[] rowSums = MyStats.rowsSums(mat);
    double meanRowSum = MyStats.mean(rowSums);
    double[] framesCanLose = new double[nRow];
    double[] framesCanGet = framesCanLose.clone();
    for (int frame = 0; frame != nRow; ++frame) {
      framesCanLose[frame] = Math.max(0, rowSums[frame] - meanRowSum);
      framesCanGet[frame] = Math.max(0, meanRowSum - rowSums[frame]);
    }
    return (long) Math.min(MyStats.sum(framesCanGet),
                           MyStats.sum(framesCanLose));
  }
  private static long getMaxSwaps(double[] framesCanLose, double[] framesCanGet) {
    return (long) Math.min(MyStats.sum(framesCanLose),
                           MyStats.sum(framesCanGet));
  }

  private static void performSwaps(Matrix mat, Matrix origMat,
                                   long nSwaps, long seed) {
    if (nSwaps == 0) {
      return;
    }
    int nRow = mat.getRowDimension(), nCol = mat.getColumnDimension();
    double[] rowSums = MyStats.rowsSums(mat);
    double meanRowSum = MyStats.mean(rowSums);
    double origMean = MyStats.mean(origMat.getRowPackedCopy());
    double[] frameWeights = MyStats.rowsMeans(origMat);
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
      double[] weightRow = MyMat.getRow(weightMat, row);
      double[] matRow = MyMat.getRow(mat, row);
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
        double[] weightRow = MyMat.getRow(weightMat, frameLosing);
        pxGiveWeights[frameLosing] = new EnumeratedIntegerDistribution(pxIndex, weightRow);
        pxGiveWeights[frameLosing].reseedRandomGenerator(seed++);
      }
      int frameGetting = frameGettingIntDist.sample();
      mat.set(frameGetting, pxLosing, mat.get(frameGetting, pxLosing) + 1);
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
  private static ImagePlus performSwaps(ImagePlus oneChImPlus, long nSwaps, int seed) {
    int width = oneChImPlus.getDimensions()[0];
    Matrix mat = MyImg.convertToMatrix(oneChImPlus);
    performSwaps(mat, mat.copy(), nSwaps, seed);
    return MyImg.convertToImagePlus(mat, width);
  }

  private static long calcIdealSwaps(ImagePlus oneChImPlus, int seed) {
    Matrix simMat = Sim.simMat(oneChImPlus, seed++);
    double lowerMeanB = MyStats.colsMeanBrightnessB(simMat);
    if (lowerMeanB <= 1) {
      return 0;
    }
    Matrix leastMat = simMat.copy();
    long maxSwaps = getMaxSwaps(simMat);
    long mostSwaps = 1, leastSwaps = 0;
    Matrix mostMat = leastMat.copy();
    performSwaps(mostMat, simMat, mostSwaps, seed++);
    while (MyStats.colsMeanBrightnessB(mostMat) > 1) {
      leastMat = mostMat.copy();
      leastSwaps = mostSwaps;
      long moreSwaps = Math.min(mostSwaps, maxSwaps - mostSwaps);
      performSwaps(mostMat, simMat, moreSwaps, seed++);
      mostSwaps += moreSwaps;
      if (mostSwaps == maxSwaps) {
        break;
      }
    }
    if (MyStats.colsMeanBrightnessB(mostMat) > 1) {
      return maxSwaps;
    }
    long middleSwaps = (leastSwaps + mostSwaps) / 2;
    Matrix middleMat = leastMat.copy();
    performSwaps(middleMat, simMat, middleSwaps - leastSwaps, seed++);
    while ((mostSwaps - middleSwaps > 0) && (middleSwaps - leastSwaps > 0)) {
      double middleMatMeanB = MyStats.colsMeanBrightnessB(middleMat);
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
      middleMat = leastMat.copy();
      performSwaps(middleMat, simMat, middleSwaps - leastSwaps, seed++);
    }
    long[] swaps = {leastSwaps, middleSwaps, mostSwaps};
    Matrix[] mats = {leastMat, middleMat, mostMat};
    double diffs1[] = new double[3];
    for (int i = 0; i != 3; ++i) {
      diffs1[i] = MyStats.absDiff(MyStats.colsMeanBrightnessB(mats[i]), 1);
    }
    return swaps[MyStats.whichMin(diffs1)];
  }

  static ImagePlus detrend(ImagePlus imPlus, int seed)
  throws DataFormatException {
    ImagePlus myImPlus = MyImg.makeMine(imPlus);
    int nCh = myImPlus.getDimensions()[2];
    // check for nCh <= 7
    ImagePlus[] channelArr = new ImagePlus[nCh];
    ImagePlus[] outChannelArr = new ImagePlus[nCh];
    for (int i = 0; i != nCh; ++i) {
      channelArr[i] = MyImg.makeMyOneCh(myImPlus, i + 1);
      long nSwaps = calcIdealSwaps(channelArr[i], seed++);
      outChannelArr[i] = performSwaps(channelArr[i], nSwaps, seed++);
    }
    if (outChannelArr.length == 1) {
      return outChannelArr[0];
    }
    return RGBStackMerge.mergeChannels(outChannelArr, false);
  }

}
