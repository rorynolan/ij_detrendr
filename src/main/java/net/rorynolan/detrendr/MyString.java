package net.rorynolan.detrendr;

class MyString {

  private static String last(String[] s) {
    int l = s.length;
    return s[l - 1];
  }

  static String insertIntoFileName(String fileName, String s) {
    int l = fileName.length();
    if (fileName.equals(".")) {
      return s;
    }
    if (fileName.matches(".*\\..*")) {
      if (l > 1 &&
              fileName.substring(0, 1).equals(".") &&
              (!fileName.substring(1, l).matches("\\."))) {
        return(fileName + "_" + s);
      }
      System.out.println(fileName.substring(l - 1, l));
      if (fileName.substring(l - 1, l).equals(".")) {
        return fileName.substring(0, l - 1) + "_" + s;
      }
      String[] dotSplit = fileName.split("\\.");
      StringBuilder out = new StringBuilder();
      int dsl = dotSplit.length;
      for (int i = 0; i != dsl - 1; ++i) {
        out.append(dotSplit[i]);
        if (i != dsl - 2) {
          out.append(".");
        }
      }
      out.append("_").append(s).append(".").append(last(dotSplit));
      return out.toString();
    }
    return fileName + s;
  }

}
