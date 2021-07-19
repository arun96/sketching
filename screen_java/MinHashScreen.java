import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;
import java.util.TreeSet;
import java.util.*;
import java.lang.*;
import java.util.Collection;

public class MinHashScreen extends Screen {

  public static void main(String[] args) throws Exception
  {
  }

  // Screen Generator for MinHash-based screen
  MinHashScreen() throws Exception
  {
    if (Settings.FIXED) {
      System.out.println("Generating MinHash-based screen with fixed size = " + Settings.FIXED_SIZE + ", using " + getHashName(Settings.HASH_TYPE) + " hash function.");
    } else {
      System.out.println("Generating MinHash-based Screen, using " + getHashName(Settings.HASH_TYPE) + " hash function.");
    }

    // Store variables
    this.targetMatches = Settings.TARGET_MATCHES;
    this.readLen = Settings.READ_LENGTH;
    this.readErr = Settings.READ_ERROR;
    this.genomeFolder = Settings.GENOME_FOLDER;
    this.k = Settings.K;
    this.genomeNames = Settings.GENOMES;
    this.numGenomes = genomeNames.length;
    this.window = 0;

    // Get genome lengths
    int[] genomeLengths = new int[numGenomes];

    // Array for the genomes
    String[] genomes = new String[numGenomes];

    // Read in the genomes, save length and the genome
    for (int i = 0; i < numGenomes; i++)
    {
      String gnm = getGenome(genomeFolder + genomeNames[i]);
      genomes[i] = gnm;
      genomeLengths[i] = gnm.length();
    }

    // Compute sketch size/ set it if fixed
    sketch_size = new int[numGenomes];
    for (int j = 0; j < numGenomes; j++)
    {
      if (Settings.FIXED) {
        sketch_size[j] = Settings.FIXED_SIZE;
      } else {
        sketch_size[j] = getSketchSize(genomeLengths[j], readLen, readErr, targetMatches, k);

        // Max sketch size
        if (sketch_size[j] > (genomeLengths[j] - k)){
          sketch_size[j] = genomeLengths[j] - k;
        }
      }
    }

    // Get sketch using genomes and sketch sizes
    // Store the sketch
    sketch_hash = new ArrayList<HashSet<Integer>>();

    // For each genome
    for (int x = 0; x < numGenomes; x++)
    {
      sketch_hash.add(getMinHashes(genomes[x], sketch_size[x], k));
    }

    if (Settings.WEIGHTED) {

      // Create weights dictionary
      System.out.println("Calculating Weights...");
      weights = new HashMap<Integer, Integer>();

      HashSet<Integer> total_sketch = new HashSet<Integer>();
      total_sketch = combineSketch(sketch_hash);

      for (int x = 0; x < numGenomes; x++)
      {
        weights = getSketchWeights(genomes[x], total_sketch, k, weights);
      }

      // Debugging
      // max_weight = Collections.max(weights.values());
      // System.out.println(max_weight);
      // System.out.println(total_sketch.size());
      // System.out.println(weights.size());
    }
  }

  // Option for pre-generated screens
  MinHashScreen(String[] genomes, ArrayList<HashSet<Integer>> sketch, Map<Integer, Integer> loaded_weights) throws Exception {
    System.out.println("Creating MinHash-based Screen...");

    weights = new HashMap<Integer, Integer>();
    sketch_hash = new ArrayList<HashSet<Integer>>();

    this.targetMatches = Settings.TARGET_MATCHES;
    this.readLen = Settings.READ_LENGTH;
    this.readErr = Settings.READ_ERROR;
    this.genomeFolder = Settings.SCREEN_LOCATION;
    this.k = Settings.K;
    this.genomeNames = genomes;
    this.numGenomes = genomeNames.length;
    this.window = 0;
    sketch_hash = sketch;
    weights = loaded_weights;
  }

}
