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
import java.util.HashMap;

// GUAVA
import com.google.common.hash.*;
import com.google.common.hash.Hashing;
import com.google.common.hash.HashFunction;

// CLUSTERING
// Clust4j
import com.clust4j.algo.*;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;

// Heirachical Clustering Java
import com.apporiented.algorithm.clustering.*;
import com.apporiented.algorithm.clustering.visualization.*;


public class ClusterGenerator{

  // int [] cluster_assignments;
  // ArrayList<ArrayList<Integer>> clusters;

  Cluster cluster;
  HashMap<String, ArrayList<String>> cluster_sketch_map;
  HashMap<String, Integer> genome_sketch_map;
  HashMap<String, Integer> cluster_height_map;

  ClusterGenerator() throws Exception{

    int numGenomes = Settings.GENOMES.length;
    ArrayList<HashSet<Integer>> sketch = new ArrayList<HashSet<Integer>>();

    // Mapping between the genome names and the sketch number
    genome_sketch_map = new HashMap<String, Integer>();

    // Read in the genomes and sketch them
    for (int i = 0; i < numGenomes; i++)
    {
      String gnm = getGenome(Settings.GENOME_FOLDER + Settings.GENOMES[i]);
      genome_sketch_map.put(Settings.GENOMES[i], i);
      sketch.add(getMinHashes(gnm, Settings.CLUSTER_SKETCH_SIZE, Settings.K));
    }

    // Similarity matrix
    double[][] similarity_matrix = new double[numGenomes][numGenomes];
    for (int x = 0; x < numGenomes; x++) {
      for (int y = 0; y < numGenomes; y++) {
        HashSet<Integer> intersection = new HashSet<Integer>(sketch.get(x));
        intersection.retainAll(sketch.get(y));
        similarity_matrix[x][y] = (double) intersection.size() / (double) Settings.CLUSTER_SKETCH_SIZE;
      }
    }
    // See the matrix
    // printMatrix(similarity_matrix);

    // Clustering - clust4j <----- NOT IN USE
    // final Array2DRowRealMatrix mat = new Array2DRowRealMatrix(similarity_matrix);
    // HierarchicalAgglomerativeParameters cluster_params = new HierarchicalAgglomerativeParameters();
    // cluster_params.setNumClusters(Settings.NUM_CLUSTERS);
    // HierarchicalAgglomerative a = cluster_params.fitNewModel(mat);
    // cluster_assignments = a.getLabels();
    // System.out.println(Arrays.toString(cluster_assignments));

    // Clustering - HCJ
    ClusteringAlgorithm alg = new DefaultClusteringAlgorithm();
    cluster = alg.performClustering(similarity_matrix, Settings.GENOMES, new SingleLinkageStrategy());
    // Prints the hierarchy
    // cluster.toConsole(2);

    cluster_sketch_map = new HashMap<String, ArrayList<String>>();

    // TODO - map between clusters and how much the sketches in them need to shrunk
    cluster_height_map = new HashMap<String, Integer>();

    // Gets the sketches in each cluster
    ArrayList<String> top_cluster_sketches = getClusterSketches(cluster, cluster_sketch_map, cluster_height_map);
    //System.out.println(top_cluster_sketches.toString());
    //System.out.println(cluster_sketch_map.toString());
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
    String forward_mer = seq.toUpperCase();
    String reversed_mer = reverseComplement(forward_mer);
    String selected_mer = "";
    if (forward_mer.compareTo(reversed_mer) > 0){
      selected_mer = reversed_mer;
    } else {
      selected_mer = forward_mer;
    }
    return selected_mer;
  }

  // Get the reverse complement of a sequence
  String reverseComplement(String sequence)
  {
    String reversed_tmp = sequence.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").toUpperCase();
    // String reversed_tmp = sequence.toUpperCase();
    // reversed_tmp = reversed_tmp.replace("A", "t").replace("T", "t").replace("C", "g").replace("G", "c");
    // reversed_tmp = reversed_tmp.toUpperCase();
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

  ArrayList<String> getClusterSketches(Cluster cluster, HashMap<String, ArrayList<String>> map, HashMap<String, Integer> height_map) {
    ArrayList<String> sketches = new ArrayList<String>();
    if (cluster.isLeaf()) {
      sketches.add(cluster.getName());
      // map.put(cluster.getName(), sketches);

      //TODO
      // Leaves have height 1
      height_map.put(cluster.getName(), 1);

    }
    else {
      List<Cluster> children = cluster.getChildren();
      int[] child_heights = new int[children.size()];
      for (int i = 0; i < children.size(); i++) {
        ArrayList<String> child_sketches = getClusterSketches(children.get(i), map, height_map);
        sketches.addAll(child_sketches);
        child_heights[i] = height_map.get(children.get(i).getName());
      }
      map.put(cluster.getName(), sketches);

      // Get the height of the cluster
      int height = Arrays.stream(child_heights).max().getAsInt() + 1;
      height_map.put(cluster.getName(), height);
    }
    return sketches;
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
