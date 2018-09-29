package net.rorynolan.detrendr;

import Jama.Matrix;
import ij.ImagePlus;
import org.apache.commons.math3.distribution.PoissonDistribution;

import java.util.zip.DataFormatException;

class Sim {

  private static Matrix poisMat(int nPx, double[] frameMeans,
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

  static Matrix simMat(ImagePlus oneChImPlus, int seed)
          throws DataFormatException {
    oneChImPlus = MyImg.assertOneChManyFrames(oneChImPlus, "simMat");
    Matrix oneChMat = MyImg.convertToMatrix(oneChImPlus);
    int oneChMatNCol = oneChMat.getColumnDimension();
    int oneChMatNRow = oneChMat.getRowDimension();
    int nNonZeroCols = 0;
    for (int col = 0; col != oneChMatNCol; ++col) {
      if (!MyStats.allAreNum(MyMat.getCol(oneChMat, col), 0)) {
        ++nNonZeroCols;
      }
    }
    double[] frameMeans = MyStats.rowsSums(oneChMat);
    for (int frame = 0; frame != oneChMatNRow; ++frame) {
      frameMeans[frame] /= nNonZeroCols;
    }
    return poisMat(nNonZeroCols, frameMeans, seed);
  }

}
