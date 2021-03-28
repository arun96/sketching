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

// Heirachical Clustering Java
import com.apporiented.algorithm.clustering.*;
import com.apporiented.algorithm.clustering.visualization.*;


public class ReadClassifierCluster {

  // List of sketches
  ArrayList<HashSet<Integer>> sketch_hash;

  // Read to be processed
  String read;

  // Matches read shares with each sketch
  int[] scores;

  // Data structure to store hashed read
  ArrayList<Integer> read_hashes;

  // Parameters for deciding hashing method
  int window;

  int k;

  // True Source of the read
  int source;

  // Min number of matches for read to be classified
  int threshold;

  // Read's number
  int read_number;

  // TODO - we need a way to store the clusters
  // int[] cluster_assignments;
  // ArrayList<ArrayList<Integer>> clusters;
  Cluster cluster;
  HashMap<String, ArrayList<String>> cluster_sketch_map;
  HashMap<String, Integer> genome_sketch_map;

  // Keep track of the status of the read
  int predicted;
  int score;
  int correct;
  int incorrect;
  int insufficient;
  int tied;

  ReadClassifierCluster(ArrayList<HashSet<Integer>> sketch_hash, String read, int window, int source, int read_number, Cluster cluster, HashMap<String, ArrayList<String>> cluster_sketch_map, HashMap<String, Integer> genome_sketch_map){

    scores = new int[sketch_hash.size()];

    this.sketch_hash = sketch_hash;

    this.window = window;

    this.read = read;

    this.k = Settings.K;

    this.source = source;

    this.threshold = Settings.THRESHOLD;

    this.read_number = read_number;

    // Cluster Information
    this.cluster = cluster;
    this.cluster_sketch_map = cluster_sketch_map;
    this.genome_sketch_map = genome_sketch_map;

    // read status
    predicted = 0;
    score = 0;
    correct = 0;
    incorrect = 0;
    insufficient = 0;
    tied = 0;
  }

  // Misclassification matrix can be pieced together from saved read logs
  void classifyRead() {

    // System.out.println(source + " " + read_number);

    // Get read k-mers
    if (window > 0) {
      read_hashes = getAllMinimizers(read, window, k);
    } else {
      read_hashes = getReadKmersHash(read, k);
    }

    boolean tie = false;
    ArrayList<Integer> scores = new ArrayList<Integer>();

    // TODO - this has to be updated to compare to the clustered hashes
    // int[] read_scores = screenReadHash(sketch_hash, read_hashes);
    String [] path = screenReadCluster(sketch_hash, cluster, cluster_sketch_map, genome_sketch_map, read_hashes, scores, tie);

    predicted = genome_sketch_map.get(path[path.length - 1]);

    // TODO - figure out some way of storing the score
    score = (int) scores.get(scores.size() - 1);
    System.out.println(source + " " + predicted + " " + score);

    // Update counts
    if (predicted == source && score > threshold) {
      // Correctly classified
      correct++;
      if (tie) {
        // Tie was broken correctly
        tied++;
      }
    } else {
      if (score == 0){
        // Not enough to classify - should not happen
        insufficient++;
      } else{
        // Misclassified - update counts
        incorrect++;
        if (Settings.TRACK_MISCLASSIFIED){
          System.out.println(source + " " + read_number);
        }
        if (tie){
          // Tie, but this time broken incorrectly
          tied++;
        }
      }
    }

    if (Settings.READ_LOGGING) {
      saveReadResultsCluster(Settings.READ_LOCATION, source, read_number, path, score, predicted, source);
      // saveReadResults(Settings.READ_LOCATION, source, read_number, path, source);
    }
  }

  // ----- HELPER FUNCTIONS FOR SCREENING ------
  // Break Read into Kmers - HASH VERSION
  ArrayList<Integer> getReadKmersHash(String r, int k)
  {
    ArrayList<Integer> readMers = new ArrayList<Integer>();

    if (r == null){
      return readMers;
    }

    int readLen = r.length();

    for (int i = 0; i < readLen-k; i++)
    {
      // Slice the kmer, and chose the canonical ones
      String readMer = r.substring(i, i+k);
      String selectedMer = getCanonical(readMer);
      readMers.add(getHash(selectedMer));
    }
    return readMers;
  }

  // Screen read kmers against screen - HASH VERSION
  int[] screenReadHash(ArrayList<HashSet<Integer>> sketch, ArrayList<Integer> readMers){

    // Number of genomes in screen
    int num = sketch.size();
    // Track number of matches
    int[] scores = new int[num];
    // For each sketch in the screen
    for (int i = 0; i < num; i++)
    {
      // Get and store overlap
      int score = getOverlapHash(sketch.get(i), readMers);
      scores[i] = score;
    }
    return scores;
  }

  // Screen read kmers against clustered approach
  String[] screenReadCluster(ArrayList<HashSet<Integer>> sketch, Cluster cluster, HashMap<String, ArrayList<String>> cluster_sketch_map, HashMap<String, Integer> genome_sketch_map, ArrayList<Integer> readMers, ArrayList<Integer> scores, boolean tie){

    int pred = 0;
    int selected = 0;
    int latest_score = 0;
    boolean leaf = false;
    ArrayList<String> path_list = new ArrayList<String>();

    while(true) {

      // if (children.size() == 0) {
      //if (!cluster_sketch_map.containsKey(cluster.getName())) {
      if (cluster.isLeaf()) {
        // path_list.add(cluster.getName());
        break;
      }

      List<Cluster> children = cluster.getChildren();


      int[] matches = new int[children.size()];

      // For each child
      for (int i = 0; i < children.size(); i++) {

        String n = children.get(i).getName();
        System.out.println("Current Cluster:");
        System.out.println(n);

        if (children.get(i).isLeaf()) {

          int si = genome_sketch_map.get(n);
          int s = getOverlapHash(sketch.get(si), readMers);

          matches[i] = s;

        } else {

          // Gets the sketches under this cluster
          ArrayList<String> child_sketches = cluster_sketch_map.get(n);
          System.out.println(child_sketches.toString());

          // Get the sketch indices
          ArrayList<Integer> child_sketches_indices = new ArrayList<Integer>();
          for (int j = 0; j < child_sketches.size(); j++) {
            String sn = child_sketches.get(j);
            child_sketches_indices.add(genome_sketch_map.get(sn));
          }

          // Create combined sketch
          HashSet<Integer> cluster_sketches = new HashSet<Integer>();
          for (int k = 0; k < child_sketches_indices.size(); k++) {
            cluster_sketches.addAll(sketch.get(child_sketches_indices.get(k)));
          }

          // Compute overlap
          int s = getOverlapHash(cluster_sketches, readMers);
          // Store the overlap
          matches[i] = s;
        }
      }

      // Get the cluster to descend into
      selected = getMaxIndex(matches);

      // Checking for tie
      if (sameCounts(matches, selected) > 0) {
        tie = true;
      }

      // Descend into cluster and update
      cluster = children.get(selected);
      path_list.add(cluster.getName());
      scores.add(matches[selected]);
    }

    // Return the path that was taken
    //int[] path_list = path.stream().mapToInt(i -> i).toArray();
    String[] path = path_list.toArray(new String[0]);
    return path;
  }

  // Finds number of shared elements between a sketch and a list of query kmers - HASH VERSION
  int getOverlapHash(HashSet<Integer> sketch_set, ArrayList<Integer> readMers_list)
  {
    int overlap = 0;

    // DEBUGGING - tracking number of duplicates
    int dupe = 0;
    // DEBUGGING - track the number of duplicates k-mers
    Set<Integer> already_counted  = new HashSet<Integer>();

    // Iterate through all the k-mers in the read, compare them to the sketch set
    for (int i = 0; i < readMers_list.size(); i++) {
      if (sketch_set.contains(readMers_list.get(i))) {
        if (already_counted.contains(readMers_list.get(i))) {
          // DEBUGGING - keeping track of duplicate k-mers
          dupe++;
        } else {
          already_counted.add(readMers_list.get(i));
          overlap++;
        }
      }
    }
    return overlap;
  }

  // Check if there are any ties for top value in the list
  int sameCounts(int[] matches, int topMatchIndex){
    int m = 0;
    int topMatchVal = matches[topMatchIndex];
    for (int r = 0; r < matches.length; r++){
      if (matches[r] == topMatchVal && r != topMatchIndex){
        m++;
      }
    }
    return m;
  }

  // ----- UTILITY FUNCTIONS ------
  // Essentially argmax, but with default value 0
  // Not using argmax, so that I can add a tie-breaking mechanism in there
  int getMaxIndex(int[] numArray){
    int max = 0;
    int max_pos = 0;
    for (int i = 0; i < numArray.length; i++)
    {
      if (numArray[i] > max) {
        max = numArray[i];
        max_pos = i;
      }
    }
    return max_pos;
  }

  // Get the reverse complement of a sequence
  String reverseComplement(String sequence)
  {
    String reversed_tmp = sequence.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").toUpperCase();
    String reversed = new StringBuffer(reversed_tmp).reverse().toString();
    return reversed;
  }

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

  // ----- MINIMIZER HELPER FUNCTIONS ------
  // Gets all minimizers from a given string with specified window size and k
  ArrayList<Integer> getAllMinimizers(String g, int window_size, int k){
    // Get number of kmers and windows in this string
    int num_mers = g.length() - k + 1;
    int num_windows = g.length() - window_size + 1;

    // Initialize data structures to store kmers and minimizers
    int[] string_hashes = new int[num_mers];
    HashSet<Integer> minimizers = new HashSet<Integer>();

    // Get all the k-mer hashes from the input string
    for (int i = 0; i < num_mers; i++){
      string_hashes[i] = getHash(getCanonical(g.substring(i, i+k)));
    }

    //iterate through the windows
    for (int j = 0; j < num_windows; j++){
      // Find the minimal kmer hash for the window that starts at this point
      int min_val = string_hashes[j];
      // Look through all kmers in this window
      for (int x = j + 1; x < (j + window_size - k); x++){
        int curr = string_hashes[x];
        if (curr < min_val) {
          min_val = curr;
        }
      }
      // Add this window's minimizer
      minimizers.add(min_val);
    }
    ArrayList<Integer> minimizer_list = new ArrayList<Integer>(minimizers);
    return minimizer_list;
  }

  // ----- READ LOGGING -----

  void saveReadResults(String location, int readset, int readnumber, int[] readscores, int s) {
    String filename = location + readset + "_" + readnumber + ".log";
    try {
      PrintWriter out = new PrintWriter(new File(filename));
      out.println(Arrays.toString(readscores));
      out.println(s);
      out.close();
    } catch (Exception e) {
      // Debugging
      // System.out.println(readnumber);
      e.printStackTrace();
    }
  }

  // TODO - update this when we have decided on clusters
  void saveReadResultsCluster(String location, int readset, int readnumber, String[] path, int score, int prediction, int source) {
    String filename = location + readset + "_" + readnumber + ".log";
    try {
      PrintWriter out = new PrintWriter(new File(filename));
      out.println(Arrays.toString(path));
      out.println(score);
      out.println(prediction);
      out.println(source);
      out.close();
    } catch (Exception e) {
      // Debugging
      // System.out.println(readnumber);
      e.printStackTrace();
    }
  }

}
