import ij.ImagePlus;
import ij.WindowManager;
import net.imagej.ImageJ;

import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.command.Command;

import java.util.NoSuchElementException;
import java.util.zip.DataFormatException;

@Plugin(type = Command.class, menuPath = "Plugins>Detrendr")
public class Detrendr implements Command {
  public static void main(String[] args) throws NoSuchElementException {
    Class<?> cls = Detrendr.class;
    String url = cls.getResource("/" + cls.getName().replace('.', '/') + ".class").toString();
    int lastIdx = url.lastIndexOf('/');
    String pluginsDir = url.substring(5, lastIdx);
    System.out.println(pluginsDir);
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
        preProcessed = Thresh.stackThresh(img, autoThreshMethod, manualThresh);
        out = Detrend.detrend(preProcessed, seed);
      } else {
        out = Detrend.detrend(img.duplicate(), seed);
      }
      out.setTitle(MyString.insertIntoFileName(title, "detrended"));
      out.show();
    } catch (DataFormatException e) {
      LogService.error(e);
    }
  }

}
