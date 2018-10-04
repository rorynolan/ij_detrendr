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
      Detrendr detrendr = new Detrendr();
      Matrix imgMat = detrendr.convertToMatrix(img);
      double[] frameMeans = detrendr.rowsMeans(imgMat);
      double frameMeansVar = detrendr.var(frameMeans);
      double origMeanB = detrendr.colsMeanBrightnessB(imgMat);
      assertEquals(1.045, origMeanB, 0.001);
      for (int seed = 1; seed != 50; ++seed) {
        ImagePlus detrended = detrendr.detrend(img, seed);
        Matrix detrendedMat = detrendr.convertToMatrix(detrended);
        double[] detrendedFrameMeans = detrendr.rowsMeans(detrendedMat);
        double detrendedFrameMeansVar = detrendr.var(detrendedFrameMeans);
        assertTrue(frameMeansVar >= detrendedFrameMeansVar);
        double detrendedMeanB = detrendr.colsMeanBrightnessB(detrendedMat);
        assertEquals(1.044, detrendedMeanB, 0.002);
      }
    } catch (DataFormatException e) {
      System.out.println(e);
      assertEquals("this", "that");
    }
  }

  @Test
  void testFileNameStuff() {
    Detrendr detrendr = new Detrendr();
    assertEquals("abc_d", detrendr.insertIntoFileName("abc", "d"));
    assertEquals("d", detrendr.insertIntoFileName(".", "d"));
    assertEquals("abc_d.tif", detrendr.insertIntoFileName("abc.tif", "d"));
    assertEquals("abc_d", detrendr.insertIntoFileName("abc.", "d"));
    assertEquals("abc_d", detrendr.insertIntoFileName(".abc", "d"));
    assertEquals("abc2.3_d.tif",
            detrendr.insertIntoFileName("abc2.3.tif", "d"));
  }
}