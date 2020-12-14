import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;
import java.util.TreeSet;
import java.util.*;
import java.lang.*;
import java.util.stream.Collectors;
import java.util.Collection;
import java.nio.charset.StandardCharsets;

// GUAVA
import com.google.common.hash.*;
import com.google.common.hash.Hashing;
import com.google.common.hash.HashFunction;

public class SaveScreen {

  String[] genomeNames;
  ArrayList<HashSet<Integer>> sketch;
  int numGenomes;

  SaveScreen(ScreenGenerator sg) throws Exception {

    this.genomeNames = Settings.GENOMES;
    this.sketch = sg.sketch_hash;
    this.numGenomes = genomeNames.length;

    // Save screens to individual files
    if (Settings.INDIVIDUAL_SCREENS) {
      for (int a = 0; a < numGenomes; a++)
      {
        String filename = Settings.SCREEN_LOCATION + genomeNames[a] + ".bin";
        saveToFile(sketch.get(a), filename);

      }
    // Save screen to a single file
    } else {
        String filename = Settings.SCREEN_LOCATION + "screen" + ".bin";
        saveToFileWhole(sketch, filename);
    }

    System.out.println("Screens generated and saved!");


    // Testing loading
    // String filename = Settings.SCREEN_LOCATION + genomeNames[1] + ".bin";
    // String filename = Settings.SCREEN_LOCATION + "screen" + ".bin";
    // HashSet<Integer> s = loadScreen(filename);
    // System.out.println(s.size());
    // System.out.println(s.getClass().getName());

  }

  // Save individual genome sketch to file
  void saveToFile(HashSet<Integer> sketch, String f) throws Exception {
    FileOutputStream fos = new FileOutputStream(f);
    ObjectOutputStream oos = new ObjectOutputStream(fos);
    oos.writeObject(sketch);
    oos.close();
    fos.close();
  }

  // TODO - write this.
  void saveToFileWhole(ArrayList<HashSet<Integer>> sketch, String f) throws Exception {
    FileOutputStream fos = new FileOutputStream(f);
    ObjectOutputStream oos = new ObjectOutputStream(fos);
    oos.writeObject(sketch);
    oos.close();
    fos.close();
  }

  // // TEST - loading the saved file
  // HashSet<Integer> loadScreen(String f) throws Exception {
  //   FileInputStream fis = new FileInputStream(f);
  //   ObjectInputStream ois = new ObjectInputStream(fis);
  //   HashSet<Integer> s = (HashSet<Integer>) ois.readObject();
  //   return s;
  // }


}
