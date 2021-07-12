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
import java.nio.charset.StandardCharsets;

// GUAVA
import com.google.common.hash.*;
import com.google.common.hash.Hashing;
import com.google.common.hash.HashFunction;

public class ScreenGenerator {
  // Key Variables
  public int targetMatches;
  public int readLen;
  public double readErr;
  public String genomeFolder;
  public int k;
  public int numGenomes;
  public int[] sketch_size;
  public String[] genomeNames;
  public ArrayList<HashSet<Integer>> sketch_hash;
  public int window;

  // TODO - add option to store weight/order for MinHash
  public Map<Integer, Integer> weights;

  // Weight - need some connection between k-mer and weights
  // Regular sketch, but a dictionary matching k-mer hashes to weights
  // Map<Integer, Integer> map = new HashMap<Integer, Integer>();
  // Could create it from the sketch, instead of from the genomes themselves?

  // Order - need to store a list of tuples
  // Normal - ArrayList<HashSet<Integer>>
  // We need ArrayList of ArrayList of Tuples

  public static void main(String[] args) throws Exception
  {
  }

  // ----- UTILITY FUNCTIONS ------
  // Compute the sketch size, given a genome, read lengths, error rates and # of target matches
  int getSketchSize(int genomeLength, int readLength, double readError, int targetMatches, int k){
    // double unaffectedChance = 1.0 - readError;
    // double multiplier = Math.pow(unaffectedChance, k);
    //
    // double genomeLengthD = (double) genomeLength;
    // double readLengthD = (double) readLength;
    // double targetMatchesD = (double) targetMatches;
    //
    // double numerator = (targetMatchesD * genomeLengthD);
    // double denominator = (readLengthD * multiplier);
    //
    // double num_sketches = numerator/denominator;
    // return (int) num_sketches;

    // NEW
    // Compute number of necessary matches per read length given error rate, extend to whole genome
    double kmer_error_rate = Math.pow((1-readError), k);
    double necessary_matches = ((double) targetMatches)/kmer_error_rate;
    double genome_segements = ((double) genomeLength) / ((double) readLength);
    double total_num_sketches = necessary_matches * genome_segements;
    return (int) total_num_sketches;
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

  // Returns a user-interpretable name for the chosen hash function
  String getHashName(String hashType){
    String ret = "";
    if (hashType.equals("h")) {
      ret = "default";
    } else if (hashType.equals("mmh3")) {
      ret = "32-bit MurmurHash3";
    } else if (hashType.equals("mmh3_128")) {
      ret = "128-bit MurmurHash3";
    } else {
      ret = "UNKNOWN";
    }
    return ret;
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

  // ----- I/O HELPER FUNCTIONS ------

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

  // Loads all .fasta/.fna/.fa files in the genomes directory
  List<String> getGenomeFiles(String directory) {
    List<String> textFiles = new ArrayList<String>();
    File dir = new File(directory);
    for (File file : dir.listFiles()) {
      if (file.getName().endsWith((".fasta")) || file.getName().endsWith((".fna")) || file.getName().endsWith((".fa"))) {
        textFiles.add(file.getName());
      }
    }
    return textFiles;
  }

  // ----- MINHASH HELPER FUNCTIONS -----

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

  // TODO - function to get weights

  //TODO - an ordered minhash version too

  // ----- EXHAUSTIVE HASH HELPER FUNCTIONS -----

  // Gets n minimal hashes from a given string
  HashSet<Integer> getAllHashes(String g, int k){
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

    return hashmers_set;
  }

  // ----- MINIMIZER HELPER FUNCTIONS ------
  // Gets all minimizers from a given string
  HashSet<Integer> getAllMinimizers(String g, int window_size, int k){
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
    //Iterate through the windows
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
    return minimizers;
  }

  // ----- WEIGHTED MINHASH HELPER FUNCTIONS ------
  // TODO
  // Map<Integer, Integer> getWeights

  // ----- ORDER MINHASH HELPER FUNCTIONS ------
  // TODO
}
