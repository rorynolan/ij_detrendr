import ij.ImagePlus;
import ij.WindowManager;
import net.imagej.ImageJ;

import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.command.Command;

@Plugin(type = Command.class, menuPath = "Plugins>Detrendr")
public class Detrendr implements Command {
  public static void main(String[] args) {
    Class<?> cls = Detrendr.class;
    String url = cls.getResource("/" + cls.getName().replace('.', '/') + ".class").toString();
    int lastIdx = url.lastIndexOf('/');
    String pluginsDir = url.substring(5, lastIdx);
    System.out.println(pluginsDir);
    System.setProperty("plugins.dir", pluginsDir);
    final ImageJ ij = net.imagej.Main.launch(args);
  }

  @Parameter(label = "Thresholding",
             persist = true,
             choices = {"None", "IJDefault", "Huang", "Huang2", "Intermodes", "IsoData",
                        "Li", "MaxEntropy", "Mean", "MinErrorI", "Minimum", "Moments", "Otsu", "Percentile",
                        "RenyiEntropy", "Shanbhag", "Triangle"})
  private String tip_radius;

  @Parameter
  LogService logService;

  @Override
  public void run() {
    ImagePlus img = WindowManager.getCurrentImage();
    MyOneChImg oneCh = new MyOneChImg(img, 2);
    oneCh.show();
  }
  
}
