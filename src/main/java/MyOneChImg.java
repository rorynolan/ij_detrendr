import ij.ImagePlus;
import ij.ImageStack;

class MyOneChImg extends MyImg {
  MyOneChImg() {
    super();
  }
  MyOneChImg(ImagePlus inImg, int channel) {
    super(inImg);
    ImageStack inStack = inImg.getStack();
    int[] dim = getDimensions();
    ImageStack outStack = new ImageStack(width, height, dim[3]);
    for (int i = 1; i <= dim[3]; ++i) {
      int inSlice = inImg.getStackIndex(channel, i, 1);
      int[] arrSlice = (int[]) inStack.getPixels(inSlice);
      outStack.setPixels(arrSlice, i);
    }
    setStack(outStack, 1, dim[3], 1);
  }
}
