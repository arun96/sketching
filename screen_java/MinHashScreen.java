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

public class MinHashScreen extends ScreenGenerator {

  public static void main(String[] args) throws Exception
  {

    // Constructor

  }

  // Screen Generator for MinHash-based screen
  MinHashScreen(String gf, String[] g, int readLength, double readError, int tm, int kmer, String hashType) throws Exception
  {
    System.out.println("Generating MinHash-based Screen, using " + getHashName(hashType) + " hash function.");
    this.screenType = "v";

    // Store variables
    this.targetMatches = tm;
    this.readLen = readLength;
    this.readErr = readError;
    this.genomeFolder = gf;
    this.k = kmer;
    this.numGenomes = g.length;
    this.genomeNames = g;
    this.window = 0;

    // Get genome lengths
    int[] genomeLengths = new int[numGenomes];

    // Array for the genomes
    String[] genomes = new String[numGenomes];

    // Read in the genomes, save length and the genome
    for (int i = 0; i < numGenomes; i++)
    {
      String gnm = getGenome(genomeFolder + g[i]);
      genomes[i] = gnm;
      genomeLengths[i] = gnm.length();
    }

    // Compute sketch size
    sketch_size = new int[numGenomes];
    for (int j = 0; j < numGenomes; j++)
    {
      sketch_size[j] = getSketchSize(genomeLengths[j], readLen, readErr, targetMatches, k);
    }

    // Get sketch using genomes and sketch sizes

    // Store the sketch
    sketch = new ArrayList<HashSet<String>>();
    sketch_hash = new ArrayList<HashSet<Integer>>();

    // For each genome
    for (int x = 0; x < numGenomes; x++)
    {
      // Row corresponding to this genome
      sketch.add(new HashSet<String>());
      sketch_hash.add(getMinHashes(genomes[x], sketch_size[x], k, hashType));
    }
  }

  // Screen Generator for fixed-size MinHash-based sketches
  MinHashScreen(String gf, String[] g, int readLength, double readError, int kmer, String fixed, int fixedSize, String hashType) throws Exception
  {
    System.out.println("Generating MinHash-based screen with fixed size = " + fixedSize + ", using " + getHashName(hashType) + " hash function.");
    this.screenType = "f";

    // Store variables
    this.targetMatches = 0;
    this.readLen = readLength;
    this.readErr = readError;
    this.genomeFolder = gf;
    this.k = kmer;
    this.numGenomes = g.length;
    this.genomeNames = g;
    this.window = 0;

    // Get genome lengths
    int[] genomeLengths = new int[numGenomes];

    // Array for the genomes
    String[] genomes = new String[numGenomes];

    // Read in the genomes, save length and the genome
    for (int i = 0; i < numGenomes; i++)
    {
      String gnm = getGenome(genomeFolder + g[i]);
      genomes[i] = gnm;
      genomeLengths[i] = gnm.length();
    }

    // Fixed sketch size for all genomes!
    sketch_size = new int[numGenomes];
    for (int j = 0; j < numGenomes; j++)
    {
      sketch_size[j] = fixedSize;
    }

    // Get sketch using genomes and the fixed size

    // Store the sketch
    sketch = new ArrayList<HashSet<String>>();
    sketch_hash = new ArrayList<HashSet<Integer>>();

    // For each genome
    for (int x = 0; x < numGenomes; x++)
    {
      // Row corresponding to this genome
      sketch.add(new HashSet<String>());
      // sketch_hash.add(new HashSet<Integer>());
      sketch_hash.add(getMinHashes(genomes[x], sketch_size[x], k, hashType));
    }
  }
}
