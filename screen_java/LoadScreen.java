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

public class LoadScreen {

  String[] genomeNames;
  ArrayList<HashSet<Integer>> sketch;
  Map<Integer,Integer> weights;

  LoadScreen() throws Exception {

    System.out.println("Loading screens...");

    // Load all .bin files in the directory, except for weights
    List<String> sketches = new ArrayList<String>();
    File screenDir = new File(Settings.SCREEN_LOCATION);
    for (File f : screenDir.listFiles()) {
      if ((f.getName().endsWith(("bin"))) && !(f.getName().startsWith(("weights")))) {
        sketches.add(f.getName());
      }
    }
    Collections.sort(sketches);
    this.genomeNames = new String[sketches.size()];
    genomeNames = sketches.toArray(genomeNames);

    // READ IN SKETCHES FROM THESE FILES
    int num_sketches = genomeNames.length;

    sketch = new ArrayList<HashSet<Integer>>();

    for (int i = 0; i < num_sketches; i++){
      String filename = Settings.SCREEN_LOCATION + genomeNames[i];
      sketch.add(loadSketch(filename));
    }

    System.out.println("Screens loaded!");

    weights = new HashMap<Integer, Integer>();
    if (Settings.WEIGHTED) {
      String filename = Settings.SCREEN_LOCATION + "weights.bin";
      weights = loadWeights(filename);
      System.out.println("Weights loaded!");
      // System.out.println(weights.size());
    }

  }

  // Loads a single sketch
  HashSet<Integer> loadSketch(String f) throws Exception {
    FileInputStream fis = new FileInputStream(f);
    ObjectInputStream ois = new ObjectInputStream(fis);
    HashSet<Integer> s = (HashSet<Integer>) ois.readObject();
    return s;
  }

  // Loads the weights for the generated screen
  Map<Integer,Integer> loadWeights(String f) throws Exception {
    FileInputStream fis = new FileInputStream(f);
    ObjectInputStream ois = new ObjectInputStream(fis);
    HashMap<Integer,Integer> w = (HashMap<Integer,Integer>) ois.readObject();
    return w;
  }

}
