package net.rorynolan.detrendr;

import Jama.Matrix;

class MyMat {
  static double[] getRow(Matrix mat, int row) {
    // include a check on value of row if get time
    int nCol = mat.getColumnDimension();
    double[] out = new double[nCol];
    for (int col = 0; col != nCol; ++col) {
      out[col] = mat.get(row, col);
    }
    return out;
  }
  static double[] getCol(Matrix mat, int col) {
    // include a check on value of col if get time
    int nRow = mat.getRowDimension();
    double[] out = new double[nRow];
    for (int row = 0; row != nRow; ++row) {
      out[row] = mat.get(row, col);
    }
    return out;
  }
}
