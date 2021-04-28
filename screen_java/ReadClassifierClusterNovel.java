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


public class ReadClassifierClusterNovel {

  // List of sketches
  // ArrayList<HashSet<Integer>> sketch_hash;

  // Read to be processed
  String read;

  // Matches read shares with each sketch
  int[] scores;

  // Data structure to store hashed read
  ArrayList<Integer> read_hashes;

  // Parameters for deciding hashing method
  int window;

  int k;

  // readset number
  int readSet;

  // Min number of matches for read to be classified
  int threshold;

  // Read's number
  int read_number;


  Cluster cluster;
  HashMap<String, ArrayList<String>> cluster_sketch_map;
  HashMap<String, Integer> genome_sketch_map;
  HashMap<String, Integer> cluster_height_map;

  // Keep track of the status of the read
  int predicted;
  int score;
  int insufficient;
  int tied;

  ReadClassifierClusterNovel(String read, int window, int readSet, int read_number, Cluster cluster, HashMap<String, ArrayList<String>> cluster_sketch_map, HashMap<String, Integer> genome_sketch_map, HashMap<String, Integer> cluster_height_map, int sketch_size){

    scores = new int[sketch_size];

    // this.sketch_hash = sketch_hash;

    this.window = window;

    this.read = read;

    this.k = Settings.K;

    this.readSet = readSet;

    this.threshold = Settings.THRESHOLD;

    this.read_number = read_number;

    // Cluster Information
    this.cluster = cluster;
    this.cluster_sketch_map = cluster_sketch_map;
    this.genome_sketch_map = genome_sketch_map;
    this.cluster_height_map = cluster_height_map;

    // read status
    predicted = 0;
    score = 0;
    insufficient = 0;
    tied = 0;
  }

  // Misclassification matrix can be pieced together from saved read logs
  void classifyRead(ArrayList<HashSet<Integer>> sketch_hash, HashMap<String, HashSet<Integer>> cluster_map) {

    // Get read k-mers
    if (window > 0) {
      read_hashes = getAllMinimizers(read, window, k);
    } else {
      read_hashes = getReadKmersHash(read, k);
    }

    boolean tie = false;
    ArrayList<Integer> scores = new ArrayList<Integer>();

    String [] path = screenReadCluster(sketch_hash, cluster_map, cluster, cluster_sketch_map, genome_sketch_map, read_hashes, scores, tie);

    predicted = genome_sketch_map.get(path[path.length - 1]);

    score = (int) scores.get(scores.size() - 1);

    // System.out.println(predicted + " " + score + " " + readSet + " " + read_number);

    // Check if insufficient or tied
    if (score == 0 || score < threshold) {
      // Not enough to classify - should not happen often
      insufficient++;
    }
    if (tie) {
      tied++;
    }

    if (Settings.READ_LOGGING) {
      saveReadResultsCluster(Settings.READ_LOCATION, readSet, read_number, path, score, predicted, sketch_hash.size());
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
  String[] screenReadCluster(ArrayList<HashSet<Integer>> sketch, HashMap<String, HashSet<Integer>> cluster_map, Cluster cluster, HashMap<String, ArrayList<String>> cluster_sketch_map, HashMap<String, Integer> genome_sketch_map, ArrayList<Integer> readMers, ArrayList<Integer> scores, boolean tie){

    int pred = 0;
    int selected = 0;
    int latest_score = 0;
    boolean leaf = false;
    ArrayList<String> path_list = new ArrayList<String>();

    while(true) {

      if (cluster.isLeaf()) {
        // path_list.add(cluster.getName());
        break;
      }

      List<Cluster> children = cluster.getChildren();


      int[] matches = new int[children.size()];

      // For each child
      for (int i = 0; i < children.size(); i++) {

        String n = children.get(i).getName();
        // System.out.println("Current Cluster:");
        // System.out.println(n);

        if (children.get(i).isLeaf()) {

          int si = genome_sketch_map.get(n);
          int s = getOverlapHash(sketch.get(si), readMers);

          matches[i] = s;

        } else {

          int ss = getOverlapHash(cluster_map.get(n), readMers);
          matches[i] = ss;
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

  // TODO - finalize
  void saveReadResultsCluster(String location, int readSet, int readnumber, String[] path, int score, int prediction, int screen_size) {
    String filename = location + readSet + "_" + readnumber + ".log";
    try {
      PrintWriter out = new PrintWriter(new File(filename));
      out.println(screen_size);
      out.println(Arrays.toString(path));
      out.println(prediction);
      out.println(score);
      out.close();
    } catch (Exception e) {
      // Debugging
      // System.out.println(readnumber);
      e.printStackTrace();
    }
  }

}
