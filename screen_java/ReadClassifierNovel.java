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


public class ReadClassifierNovel {

  // Read to be processed
  String read;

  // Matches read shares with each sketch
  int[] scores;

  // Data structure to store hashed read
  ArrayList<Integer> read_hashes;

  // Parameters for deciding hashing method
  int window;

  int k;

  // Min number of matches for read to be classified
  int threshold;

  // Read's number and Read Set
  int read_number;
  int read_set;

  // Keep track of the status of the read
  int predicted;
  int score;
  int insufficient;
  int tied;

  boolean filtered_out;

  ReadClassifierNovel(int readSet, String read, int window, int read_number, int sketch_size) {

    scores = new int[sketch_size];

    // this.sketch_hash = sketch_hash;

    this.window = window;

    this.read = read;

    this.k = Settings.K;

    this.threshold = Settings.THRESHOLD;

    this.read_number = read_number;

    this.read_set = readSet;

    this.filtered_out = false;

    // read status
    predicted = 0;
    score = 0;
    insufficient = 0;
    tied = 0;
  }

  // Misclassification matrix can be pieced together from saved read logs
  void classifyRead(ArrayList<HashSet<Integer>> sketch_hash) {

    // System.out.println(source + " " + read_number);

    // TODO - generalize filtering
    // Read Filtering
    if (Settings.FILTER_READS) {
      filtered_out = read_filtering(read, Settings.READ_LENGTH);
      if (filtered_out) {
        return;
      }
    }

    // Get read k-mers
    if (window > 0) {
      read_hashes = getAllMinimizers(read, window, k);
    } else {
      read_hashes = getReadKmersHash(read, k);
    }

    int[] read_scores = screenReadHash(sketch_hash, read_hashes);

    predicted = getMaxIndex(read_scores);

    score = read_scores[predicted];
    // System.out.println(score);

    // Check if tied or insufficient
    if (read_scores[predicted] == 0 || read_scores[predicted] < threshold){
      insufficient++;
    } else {
      if (sameCounts(read_scores, predicted) > 0) {
        tied++;
      }
    }

    // saveReadResults(Settings.READ_LOCATION, read_set, read_number, predicted, score, tied);
    saveReadResults(Settings.READ_LOCATION, read_set, read_number, read_scores);
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

  // Function for filtering "bad" reads out
  // TODO: Update for read qualities and % of bad lengths
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
  // Save read details to a file
  void saveReadResults(String location, int readset, int readnumber, int[] readscores) {
    String filename = location + readset + "_" + readnumber + ".log";
    try {
      PrintWriter out = new PrintWriter(new File(filename));
      out.println(Arrays.toString(readscores));
      out.close();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
