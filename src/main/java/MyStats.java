import Jama.Matrix;
import ij.ImagePlus;
import ij.ImageStack;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import java.util.zip.DataFormatException;

class MyStats {

  static long sum(int[] data) {
    long sum = 0;
    for(int a : data)
      sum += a;
    return sum;
  }
  static double sum(double[] data) {
    double sum = 0.0;
    for(double a : data)
      sum += a;
    return sum;
  }

  static double mean(int[] data) {
    return sum(data) / data.length;
  }
  static double mean(double[] data) {
    return sum(data) / data.length;
  }

  static double var(int[] data) {
    double mean = mean(data);
    double temp = 0;
    for(double a : data)
      temp += (a - mean) * (a - mean);
    return temp / (data.length - 1);
  }

  static double brightnessB(int[] data) {
    double mean = mean(data);
    double temp = 0;
    for(double a : data)
      temp += (a - mean) * (a - mean);
    double var = temp / (data.length - 1);
    return var / mean;
  }
  static double brightnessB(SummaryStatistics sumStat) {
    return sumStat.getVariance() / sumStat.getMean();
  }
  static Matrix brightnessB(ImagePlus imPlus)
          throws DataFormatException {
    ImagePlus myImPlus = MyImg.makeMine(imPlus);
    int[] d = myImPlus.getDimensions();
    if (d[2] != 1) {
      DataFormatException e = new DataFormatException(
              "Can only calculate the brightness of a 1-channel ImagePlus. " +
                      "\n  * You have attempted to calculate the brightness " +
                      " of a " + d[2] + "-channel ImagePlus."
      );
      throw e;
    }
    ImageStack stack = myImPlus.getStack();
    Matrix out = new Matrix(d[1], d[0]);
    SummaryStatistics sumStat = new SummaryStatistics();
    for (int x = 0; x != d[0]; ++x) {
      for (int y = 0; y != d[1]; ++y) {
        for (int z = 0; z != d[3]; ++z) {
          sumStat.addValue(stack.getVoxel(x, y, z));
        }
        out.set(y, x, brightnessB(sumStat));
        sumStat.clear();
      }
    }
    return out;
  }

  static double[] colsBrightnessB(Matrix mat) {
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

  static double colsMeanBrightnessB(Matrix mat) {
    return mean(colsBrightnessB(mat));
  }

  static double[] rowsSums(Matrix mat) {
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
  static double[] colsSums(Matrix mat) {
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
  static double[] rowsMeans(Matrix mat) {
    int nRow = mat.getRowDimension();
    int nCol = mat.getColumnDimension();
    double[] out = rowsSums(mat);
    for (int row = 0; row != nRow; ++row) {
      out[row] /= nCol;
    }
    return out;
  }
  static double[] colsMeans(Matrix mat) {
    int nRow = mat.getRowDimension();
    int nCol = mat.getColumnDimension();
    double[] out = colsSums(mat);
    for (int col = 0; col != nCol; ++col) {
      out[col] /= nRow;
    }
    return out;
  }

  static boolean allAreNum(double[] data, double num) {
    for (double d : data) {
      if (d != num) {
        return false;
      }
    }
    return true;
  }

  static double absDiff(double x, double y) {
    return Math.abs(x - y);
  }
  static int whichMin(double[] x) {
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

}