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
  int numGenomes;

  SaveScreen(Screen sg) throws Exception {

    this.genomeNames = Settings.GENOMES;
    this.numGenomes = genomeNames.length;

    // Save screens to individual files
    for (int a = 0; a < numGenomes; a++)
    {
      String filename = Settings.SCREEN_LOCATION + genomeNames[a] + ".bin";
      saveToFile(sg.sketch_hash.get(a), filename);

    }
    System.out.println("Screens saved!");

    // Save weights
    if (Settings.WEIGHTED) {
      String weights_filename = Settings.SCREEN_LOCATION + "weights.bin";
      saveWeightsToFile(sg.weights, weights_filename);
      System.out.println("Weights saved!");
    }

    // Save Order - TODO

    // Clusters - TODO

    // Experiment parameters
    String paramsFile = Settings.SCREEN_LOCATION + "params.txt";
    saveParams(sg, Settings.SCREEN_TYPE, Settings.HASH_TYPE, paramsFile);

  }

  // Save individual genome sketch to file
  void saveToFile(HashSet<Integer> sketch, String f) throws Exception {
    FileOutputStream fos = new FileOutputStream(f);
    ObjectOutputStream oos = new ObjectOutputStream(fos);
    oos.writeObject(sketch);
    oos.close();
    fos.close();
  }

  // Save weights to file
  void saveWeightsToFile(Map<Integer, Integer> weights, String f) throws Exception {
    FileOutputStream fos = new FileOutputStream(f);
    ObjectOutputStream oos = new ObjectOutputStream(fos);
    oos.writeObject(weights);
    oos.close();
    fos.close();
  }

  void saveParams(Screen sg, String screen_type, String hash_type, String f) throws Exception {
    PrintWriter out = new PrintWriter(new File(f));
    out.println("Screen Details:");
    out.println("Screen Type: " + screen_type);
    out.println("Hash Function used: " + hash_type);
    out.println("K: " + sg.k);
    out.println("Target Matches: " + sg.targetMatches);
    out.println("Read Length: " + sg.readLen);
    out.println("Read Error: " + sg.readErr);
    out.println("Window Size: " + sg.window);
    out.close();
  }

  // TODO - delete when ready
  // void saveToFileWhole(ArrayList<HashSet<Integer>> sketch, String f) throws Exception {
  //   FileOutputStream fos = new FileOutputStream(f);
  //   ObjectOutputStream oos = new ObjectOutputStream(fos);
  //   oos.writeObject(sketch);
  //   oos.close();
  //   fos.close();
  // }


}
