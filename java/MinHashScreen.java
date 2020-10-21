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
  MinHashScreen(String gf, String[] g, int readLength, double readError, int tm, int kmer) throws Exception
  {
    System.out.println("Generating Hash-based Screen...");
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
      // System.out.println(genomeLengths[i]);
    }


    // Compute sketch size
    sketch_size = new int[numGenomes];
    for (int j = 0; j < numGenomes; j++)
    {
      sketch_size[j] = getSketchSize(genomeLengths[j], readLen, readErr, targetMatches, k);
    }

    // Get sketch using genomes and sketch sizes
    Random r = new Random();

    // Store the sketch
    sketch = new ArrayList<ArrayList<String>>();

    // For each genome
    for (int x = 0; x < numGenomes; x++)
    {

      // TODO - Finish the hashcode based method
      // // Row corresponding to this genome
      // sketch.add(new ArrayList<String>());
      //
      // // Set for tracking which kmers are in sketch
      // Set<String> sketch_set  = new HashSet<String>();
      //
      // // Store the hash values for the whole genome
      // List<Integer> hash_values = new new ArrayList<>();
      //
      // int p = 0;
      //
      // // Build the list of
      // while (p < genomeLengths[x] - kmer)
      // {
      //   String mer = genomes[x].substring(p, p+kmer);
      //   String reversed_mer = reverseComplement(mer);
      //   String selected_mer = "";
      //   if (mer.compareTo(reversed_mer) > 0){
      //     selected_mer = reversed_mer;
      //   } else {
      //     selected_mer = mer;
      //   }
      //
      //   if (sketch_set.contains(selected_mer)){
      //     p++;
      //     continue;
      //   // Add to sketch, increment count
      //   } else {
      //     sketch_set.add(selected_mer);
      //     hash_values.add(selected_mer.hashCode())
      //     p++;
      // }

      // TODO - Turn this into a hash function. For now, it's just random start positions.
      // Row corresponding to this genome
      sketch.add(new ArrayList<String>());

      // Set for tracking which kmers are in sketch
      Set<String> sketch_set  = new HashSet<String>();


      // For each sketch in a given genome (there are a total of sketch_size[x] of them)
      int p = 0;

      while (p < sketch_size[x]){
        // Get k-mer
        int start = r.nextInt(genomeLengths[x]-k+1);
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

        // If it is already in the sketch
        if (sketch_set.contains(selected_mer)){
          continue;
        // Add to sketch, increment count
        } else {
          sketch_set.add(selected_mer);
          sketch.get(x).add(selected_mer);
          p++;
        }
      }
    }
  }

  // Screen Generator for fixed size sketches
  MinHashScreen(String gf, String[] g, int readLength, double readError, int kmer, String fixed, int fixedSize) throws Exception
  {
    System.out.println("Generating Fixed Size Screen...");
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
      // System.out.println(genomeLengths[i]);
    }

    // Fixed sketch size for all genomes!
    sketch_size = new int[numGenomes];
    for (int j = 0; j < numGenomes; j++)
    {
      sketch_size[j] = fixedSize;
    }

    // Get sketch using genomes and the fixed size
    // TODO - Turn this into a hash function. For now, it's just random start positions.
    Random r = new Random();

    // Store the sketch
    sketch = new ArrayList<ArrayList<String>>();

    // For each genome
    for (int x = 0; x < numGenomes; x++)
    {
      // Row corresponding to this genome
      sketch.add(new ArrayList<String>());

      // Set for tracking which kmers are in sketch
      Set<String> sketch_set  = new HashSet<String>();

      // For each sketch in a given genome (there are a total of fixedSize of them)
      int p = 0;
      while (p < sketch_size[x]){
        // Get k-mer
        int start = r.nextInt(genomeLengths[x]-k+1);
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

        // If it is already in the sketch
        if (sketch_set.contains(selected_mer)){
          continue;
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
