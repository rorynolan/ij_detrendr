package net.rorynolan.detrendr;

import Jama.Matrix;
import ij.ImagePlus;
import org.junit.jupiter.api.Test;

import java.util.zip.DataFormatException;

import static org.junit.jupiter.api.Assertions.*;

class DetrendTest {

  @Test
  void testDetrendingAndBrightness() {
    try {
      String nandbImgUrl = "https://github.com/rorynolan/nandb/raw/master/inst/extdata/50.tif";
      ImagePlus img = new ImagePlus(nandbImgUrl);
      Matrix imgMat = MyImg.convertToMatrix(img);
      double[] frameMeans = MyStats.rowsMeans(imgMat);
      double frameMeansVar = MyStats.var(frameMeans);
      double origMeanB = MyStats.colsMeanBrightnessB(imgMat);
      assertEquals(1.045, origMeanB, 0.001);
      for (int seed = 1; seed != 50; ++seed) {
        ImagePlus detrended = Detrend.detrend(img, seed);
        Matrix detrendedMat = MyImg.convertToMatrix(detrended);
        double[] detrendedFrameMeans = MyStats.rowsMeans(detrendedMat);
        double detrendedFrameMeansVar = MyStats.var(detrendedFrameMeans);
        assertTrue(frameMeansVar >= detrendedFrameMeansVar);
        double detrendedMeanB = MyStats.colsMeanBrightnessB(detrendedMat);
        assertEquals(1.044, detrendedMeanB, 0.002);
      }
    } catch (DataFormatException e) {
      System.out.println(e);
      assertEquals("this", "that");
    }
  }
}