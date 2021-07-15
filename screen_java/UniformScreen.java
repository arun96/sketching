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

public class UniformScreen extends Screen {

  public static void main(String[] args) throws Exception
  {
  }

  //Screen Generator for Uniform Sampling
  UniformScreen() throws Exception
  {
    if (Settings.FIXED) {
      System.out.println("Generating Uniformly Sampled Screen with fixed size = " + Settings.FIXED_SIZE + ", using " + getHashName(Settings.HASH_TYPE) + " hash function.");
    } else {
      System.out.println("Generating Uniformly Sampled Screen, using " + getHashName(Settings.HASH_TYPE) + " hash function.");
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
      // System.out.println(genomeLengths[i]);
    }

    // Compute sketch size/set it if fixed
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
      // Row corresponding to this genome
      sketch_hash.add(new HashSet<Integer>());

      // distance between sketched k-mers
      int spacing = (int) genomeLengths[x]/sketch_size[x];

      // Jump through the genome and store the k-mers at those positions
      int p = 0;
      while (p < (genomeLengths[x]-k)){
        // Get k-mer
        int start = p;
        int end = start + k;
        String mer = genomes[x].substring(start, end);

        String selected_mer = getCanonical(mer);
        int hash_val = getHash(selected_mer);

        sketch_hash.get(x).add(hash_val);
        p = p + spacing;
      }
    }
  }

  // Option for pre-generated screens
  UniformScreen(String[] genomes, ArrayList<HashSet<Integer>> sketch) throws Exception
  {
    System.out.println("Creating Uniformly Sampled Screen...");
    this.targetMatches = Settings.TARGET_MATCHES;
    this.readLen = Settings.READ_LENGTH;
    this.readErr = Settings.READ_ERROR;
    this.genomeFolder = Settings.SCREEN_LOCATION;
    this.k = Settings.K;
    this.genomeNames = genomes;
    this.numGenomes = genomeNames.length;
    this.window = 0;
    sketch_hash = sketch;
  }

}
