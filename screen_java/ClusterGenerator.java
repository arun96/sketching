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
import java.util.stream.Collectors;
import java.util.Collection;
import java.nio.charset.StandardCharsets;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;

// GUAVA
import com.google.common.hash.*;
import com.google.common.hash.Hashing;
import com.google.common.hash.HashFunction;

// CLUSTERING
import com.clust4j.algo.*;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;

public class ClusterGenerator{


  ArrayList<ArrayList<Integer>> clusters;
  double[][] similarity_matrix;

  ClusterGenerator() throws Exception{

    System.out.println("Clustering Genomes...");

    int numGenomes = Settings.GENOMES.length;
    ArrayList<HashSet<Integer>> sketch = new ArrayList<HashSet<Integer>>();

    // Read in the genomes
    for (int i = 0; i < numGenomes; i++)
    {
      String gnm = getGenome(Settings.GENOME_FOLDER + Settings.GENOMES[i]);
      sketch.add(getMinHashes(gnm, Settings.CLUSTER_SKETCH_SIZE, Settings.K));
    }

    // Similarity matrix
    similarity_matrix = new double[numGenomes][numGenomes];
    for (int x = 0; x < numGenomes; x++) {
      for (int y = 0; y < numGenomes; y++) {
        HashSet<Integer> intersection = new HashSet<Integer>(sketch.get(x));
        intersection.retainAll(sketch.get(y));
        similarity_matrix[x][y] = (double) intersection.size() / (double) Settings.CLUSTER_SKETCH_SIZE;
      }
    }

    // See the matrix
    // printMatrix(similarity_matrix);

    // Clustering - clust4j
    System.out.println("Clustering genomes...");
    final Array2DRowRealMatrix mat = new Array2DRowRealMatrix(similarity_matrix);
    HierarchicalAgglomerativeParameters cluster_params = new HierarchicalAgglomerativeParameters();
    cluster_params.setNumClusters(Settings.NUM_CLUSTERS);
    HierarchicalAgglomerative a = cluster_params.fitNewModel(mat);
    int[] results = a.getLabels();
    System.out.println(Arrays.toString(results));
    System.out.println("Assigned clusters - combining sketches...");


  }

  // ----- DISPLAY FUNCTION ------

  void printMatrix(double[][] m) throws Exception {
    for (double[] x : m) {
      for (double y : x) {
        System.out.print(y + " ");
      }
      System.out.println();
    }
  }
  // ----- UTILITY FUNCTIONS ------

  // Get the lexicographically smaller of a k-mer and its reverse complement
  String getCanonical(String seq){
    String reversed_mer = reverseComplement(seq);
    String selected_mer = "";
    if (seq.compareTo(reversed_mer) > 0){
      selected_mer = reversed_mer;
    } else {
      selected_mer = seq;
    }
    return selected_mer;
  }

  // Generates the reverse complement of a DNA sequence
  String reverseComplement(String sequence){
    String reversed_tmp = sequence.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").toUpperCase();
    String reversed = new StringBuffer(reversed_tmp).reverse().toString();
    return reversed;
  }

  // Given a string and the hash function to be used, returns the hashed sequence
  int getHash(String seq){
    if (Settings.HASH_TYPE.equals("h")){
      return seq.hashCode();
    } else if (Settings.HASH_TYPE.equals("mmh3")) {
      int hashVal = Hashing.murmur3_32().hashString(seq, StandardCharsets.UTF_8).asInt();
      return hashVal;
    } else if (Settings.HASH_TYPE.equals("mmh3_128")) {
      int hashVal = Hashing.murmur3_128().hashString(seq, StandardCharsets.UTF_8).asInt();
      return hashVal;
    } else {
      // Default - use the hashCode version
      return seq.hashCode();
    }
  }

  // Gets n minimal hashes from a given string
  HashSet<Integer> getMinHashes(String g, int sketch_size, int k){
    int gl = g.length();
    int boundary = gl - k;
    HashSet<Integer> hashmers_set = new HashSet<Integer>();

    // Build list of hashed mers
    for (int p = 0; p < boundary; p++) {
      int start = p;
      int end = p + k;
      String curr = g.substring(start, end);
      hashmers_set.add(getHash(getCanonical(curr)));
    }

    // Sort this list
    ArrayList<Integer> hashmers = new ArrayList<Integer>(hashmers_set);
    Collections.sort(hashmers);

    HashSet<Integer> minhashvals = new HashSet<Integer>();

    // Return the n minimal ones
    for (int q = 0; q < sketch_size; q++) {
      minhashvals.add(hashmers.get(q));
    }
    return minhashvals;
  }


  // ----- I/O FUNCTIONS ------

  // Load genome from a given file
  String getGenome(String fn) throws Exception {
    Scanner input = new Scanner(new FileInputStream(new File(fn)));
    StringBuilder sb = new StringBuilder("");
    while(input.hasNext())
    {
      String line = input.nextLine();
      if(line.length() == 0 || line.startsWith(">"))
      {
        continue;
      }
      sb.append(line.toUpperCase());
    }
    String s = sb.toString();
    return s;
  }

}
