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

public class MinimizerScreen extends ScreenGenerator{

  public static void main(String[] args) throws Exception
  {

    // Constructor

  }

  // Screen Generator for minimzer-based approach, with calculated windowsize
  MinimizerScreen(String gf, String[] g, int readLength, double readError, int tm, int kmer, String mini, boolean track, String hashType) throws Exception
  {
    // TRUE if we want to save minimizer positions
    boolean MINIMIZER_TRACKING = track;

    // Multiplier for adjusting window size - finalized on 2.0
    double multiplier = 2.0;

    System.out.println("Generating Minimizer-Based Screen with calculated Window Size, using " + getHashName(hashType) + " hash function.");
    this.screenType = "v";

    // Store variables
    this.targetMatches = tm;
    this.readLen = readLength;
    this.readErr = readError;
    this.genomeFolder = gf;
    this.k = kmer;
    this.numGenomes = g.length;
    this.genomeNames = g;

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
    int[] window_sizes = new int[numGenomes];
    for (int j = 0; j < numGenomes; j++)
    {
      sketch_size[j] = getSketchSize(genomeLengths[j], readLen, readErr, targetMatches, k);
      window_sizes[j] = (int) ((genomeLengths[j]/sketch_size[j])*multiplier);
    }

    // set the window size - agnostic of genomes
    this.window = window_sizes[0];

    // Get sketch using genomes and window sizes

    sketch = new ArrayList<HashSet<String>>();
    sketch_hash = new ArrayList<HashSet<Integer>>();

    for (int x = 0; x < numGenomes; x++)
    {
      // Row corresponding to this genome
      sketch.add(new HashSet<String>());
      HashSet<Integer> minimizers = getAllMinimizers(genomes[x], window_sizes[x], k, hashType);
      sketch_hash.add(minimizers);
    }
  }

  // Screen Generator for minimzer-based approach, with specified window size
  MinimizerScreen(String gf, String[] g, int readLength, double readError, int kmer, String mini, int windowSize, boolean track, String hashType) throws Exception
  {
    // TRUE if we want to save minimizer positions
    boolean MINIMIZER_TRACKING = track;

    System.out.println("Generating Minimizer-Based Screen with specified Window Size = " + windowSize + ", using " + getHashName(hashType) + " hash function.");

    // Store variables
    this.targetMatches = 0;
    this.readLen = readLength;
    this.readErr = readError;
    this.genomeFolder = gf;
    this.k = kmer;
    this.numGenomes = g.length;
    this.genomeNames = g;
    this.window = windowSize;

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
    int[] window_sizes = new int[numGenomes];
    for (int j = 0; j < numGenomes; j++)
    {
      window_sizes[j] = windowSize;
    }
    // set the window size - agnostic of genomes
    this.window = windowSize;

    // Get sketch using genomes and window sizes

    sketch = new ArrayList<HashSet<String>>();
    sketch_hash = new ArrayList<HashSet<Integer>>();

    for (int x = 0; x < numGenomes; x++)
    {
      // Row corresponding to this genome
      sketch.add(new HashSet<String>());
      HashSet<Integer> minimizers = getAllMinimizers(genomes[x], window_sizes[x], k, hashType);
      sketch_hash.add(minimizers);
    }
  }

}
