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

    // Constructor

  }

  //Screen Generator for Uniform Sampling
  UniformScreen(String gf, String[] g, int readLength, double readError, int tm, int kmer, String uniform) throws Exception
  {
    System.out.println("Generating Uniformly Sampled Screen...");
    this.screenType = "u";


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
    // TODO - Turn this into a hash function. For now, it's just random start positions.

    // Store the sketch
    sketch = new ArrayList<ArrayList<String>>();

    // For each genome
    for (int x = 0; x < numGenomes; x++)
    {
      // Row corresponding to this genome
      sketch.add(new ArrayList<String>());

      // Set for tracking which kmers are in sketch
      Set<String> sketch_set  = new HashSet<String>();

      // distance between sketched k-mers
      int spacing = (int) genomeLengths[x]/sketch_size[x];

      // For each sketch in a given genome (there are a total of sketch_size[x] of them)
      int p = 0;
      while (p < sketch_size[x]){
        // Get k-mer
        int start = p * spacing;
        int end = start + k;
        String mer = genomes[x].substring(start, end);

        // Get canonical kmer
        String reversed_mer = reverseComplement(mer);
        String selected_mer = "";
        if (mer.compareTo(reversed_mer) > 0){
          selected_mer = reversed_mer;
        } else {
          selected_mer = mer;
        }

        // TODO - find good solution for duplicate k-mers.

        // If it is already in the sketch
        if (sketch_set.contains(selected_mer)){
          p++;
          // continue;
        // Add to sketch, increment count
        } else {
          sketch_set.add(selected_mer);
          sketch.get(x).add(selected_mer);
          p++;
        }
      }
    }
  }

}
