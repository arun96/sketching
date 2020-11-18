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

public class UniformScreen extends ScreenGenerator {

  public static void main(String[] args) throws Exception
  {
  }

  //Screen Generator for Uniform Sampling
  UniformScreen(String gf, String[] g, int readLength, double readError, int tm, int kmer, String uniform, String hashType) throws Exception
  {
    System.out.println("Generating Uniformly Sampled Screen, using " + getHashName(hashType) + " hash function.");

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
      // System.out.println(genomeLengths[i]);
    }

    // Compute sketch size
    sketch_size = new int[numGenomes];
    for (int j = 0; j < numGenomes; j++)
    {
      sketch_size[j] = getSketchSize(genomeLengths[j], readLen, readErr, targetMatches, k);
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
        int hash_val = getHash(selected_mer, hashType);

        sketch_hash.get(x).add(hash_val);
        p = p + spacing;
      }
    }
  }

}
