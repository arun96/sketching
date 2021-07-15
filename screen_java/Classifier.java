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

// Class that contains all basic functions needed for read classification 
public class Classifier {

  // ---- Parameters for all types of classifiers (matched, novel, clustering) ----

  // Read to be processed
  public String read;

  // Matches read shares with each sketch
  public int[] scores;

  // Data structure to store hashed read
  public ArrayList<Integer> read_hashes;

  // Parameter for deciding hashing method
  public int window;

  // K
  public int k;

  // Min number of matches for read to be classified
  public int threshold;

  // Read's number
  public int read_number;

  // Keeping track of status of the read
  public int predicted;
  public int score;
  public int insufficient;
  public int tied;

  // Whether read filtering is active or not
  public boolean filtered_out;

  // ---- Parameters for classifiers where source of read is known ----

  // True Source of the read
  public int source;
  // Track the classification of the read
  public int correct;
  public int incorrect;


  // ---- Parameters for classification without known source ----
  public int read_set;

  /// Parameters for clustering
  public Cluster cluster;
  public HashMap<String, ArrayList<String>> cluster_sketch_map;
  public HashMap<String, Integer> genome_sketch_map;
  public HashMap<String, Integer> cluster_height_map;


  public static void main(String[] args) throws Exception
  {
  }

  // ----- HELPER FUNCTIONS FOR SCREENING ------

  // Breaks read into kmer hashes
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

  // Screens read kmer hashes against sketch
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

  // Finds number of shared elements between a sketch and a list of query kmers
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
    // String reversed_tmp = sequence.toUpperCase();
    // reversed_tmp = reversed_tmp.replace("A", "t").replace("T", "t").replace("C", "g").replace("G", "c");
    // reversed_tmp = reversed_tmp.toUpperCase();
    String reversed = new StringBuffer(reversed_tmp).reverse().toString();
    return reversed;
  }

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

  // Check if there are any ties for max value in a list
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

  // Function for filtering "bad" reads out
  // TODO: Update for read qualities and % of bad lengths, and to fix NPE
  boolean read_filtering(String read, int read_length) {
    boolean filter_out = false;
    String pattern = "[acgtACGT]+";
    if (read.matches(pattern)) {
      filter_out = false;
    } else {
      filter_out = true;
    }
    return filter_out;
  }

  // Count occurrences of a base in a read
  int count_bases(String read, String base) {
    int count = read.length() - read.replace(base, "").length();
    return count;
  }

}