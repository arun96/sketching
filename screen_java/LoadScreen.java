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

  LoadScreen() throws Exception {

    System.out.println("Loading screens...");

    List<String> sketches = new ArrayList<String>();
    File screenDir = new File(Settings.SCREEN_LOCATION);
    for (File readFile : screenDir.listFiles()) {
      if (readFile.getName().endsWith(("bin"))) {
        sketches.add(readFile.getName());
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
      HashSet<Integer> s = loadSketch(filename);
      this.sketch.add(s);
    }

    System.out.println("Screens loaded!");

  }

  HashSet<Integer> loadSketch(String f) throws Exception {
    FileInputStream fis = new FileInputStream(f);
    // System.out.println("File Stream Started");
    ObjectInputStream ois = new ObjectInputStream(fis);
    // System.out.println("Object Stream Started");
    HashSet<Integer> s = (HashSet<Integer>) ois.readObject();
    // System.out.println("Object Loaded");
    return s;
  }

}
