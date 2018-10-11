package net.rorynolan.detrendr;

import Jama.Matrix;
import ij.ImagePlus;
import org.junit.Test;

import java.util.zip.DataFormatException;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;


public class DetrendTest {

  @Test
  public void testDetrendingAndBrightness() {
    try {
      String nandbImgUrl = "https://github.com/rorynolan/nandb/raw/master/inst/extdata/two_ch.tif";
      ImagePlus img = new ImagePlus(nandbImgUrl);
      Detrendr detrendr = new Detrendr();
      ImagePlus imgHuanged = detrendr.stackThresh(img, "Huang", 0);
      int nCh = img.getDimensions()[2];
      Matrix[] imgMats = new Matrix[nCh];
      double[] origMeanB = new double[nCh];
      double[] expectedMeanBs = {1.035, 1.259};
      for (int i = 0; i != nCh; ++i) {
        imgMats[i] = detrendr.convertToMatrix(detrendr.makeMyOneCh(imgHuanged, i + 1));
        origMeanB[i] = detrendr.colsMeanBrightnessB(imgMats[i]);
        assertEquals(expectedMeanBs[i], origMeanB[i], 0.001);
      }
      for (int seed = 1; seed != 20; ++seed) {
        ImagePlus detrended = detrendr.detrend(img, seed);
        Matrix[] detrendedMats = new Matrix[nCh];
        double[] detrendedMeanBs = new double[nCh];
        double[] expectedDetrendedMeanBs = {1.032, 1.241};
        for (int i = 0; i != nCh; ++i) {
          detrendedMats[i] = detrendr.convertToMatrix(detrendr.makeMyOneCh(detrended, i + 1));
          detrendedMeanBs[i] = detrendr.colsMeanBrightnessB(detrendedMats[i]);
          assertEquals(expectedDetrendedMeanBs[i], detrendedMeanBs[i], 0.02);
        }
      }
    } catch (DataFormatException e) {
      System.out.println(e);
      assertEquals("this", "that");
    }
  }

  @Test
  public void testFileNameStuff() {
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