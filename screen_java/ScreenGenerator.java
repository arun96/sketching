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

  public static void main(String[] args) throws Exception
  {
  }

  // ----- UTILITY FUNCTIONS ------
  // Compute the sketch size, given a genome, read lengths, error rates and # of target matches
  int getSketchSize(int genomeLength, int readLength, double readError, int targetMatches, int k){
    double unaffectedChance = 1.0 - readError;
    double multiplier = Math.pow(unaffectedChance, k);

    double genomeLengthD = (double) genomeLength;
    double readLengthD = (double) readLength;
    double targetMatchesD = (double) targetMatches;

    double numerator = (targetMatchesD * genomeLengthD);
    double denominator = (readLengthD * multiplier);

    double num_sketches = numerator/denominator;
    return (int) num_sketches;
  }

  // Generates the reverse complement of a DNA sequence
  String reverseComplement(String sequence){
    String reversed_tmp = sequence.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").toUpperCase();
    String reversed = new StringBuffer(reversed_tmp).reverse().toString();
    return reversed;
  }

  // Given a string and the hash function to be used, returns the hashed sequence
  int getHash(String seq, String hashType){
    if (hashType.equals("h")){
      return seq.hashCode();
    } else if (hashType.equals("mmh3")) {
      // TODO - fix this
      int hashVal = Hashing.murmur3_32().hashString(seq, StandardCharsets.UTF_8).asInt();
      return hashVal;
    } else if (hashType.equals("mmh3_128")) {
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
      ret = "MurmurHash3";
    } else {
      ret = "UNKNOWN";
    }
    return ret;
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
  HashSet<Integer> getMinHashes(String g, int sketch_size, int k, String hashType){
    int gl = g.length();
    int boundary = gl - k;
    HashSet<Integer> hashmers_set = new HashSet<Integer>();

    // Build list of hashed mers
    for (int p = 0; p < boundary; p++) {
      int start = p;
      int end = p + k;
      String curr = g.substring(start, end);
      hashmers_set.add(getHash(getCanonical(curr), hashType));
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

  // ----- MINIMIZER HELPER FUNCTIONS ------
  // Gets all minimizers from a given string
  HashSet<Integer> getAllMinimizers(String g, int window_size, int k, String hashType){
    // Get number of kmers and windows in this string
    int num_mers = g.length() - k + 1;
    int num_windows = g.length() - window_size + 1;

    // Initialize data structures to store kmers and minimizers
    int[] string_hashes = new int[num_mers];
    HashSet<Integer> minimizers = new HashSet<Integer>();

    // Get all the k-mer hashes from the input string
    for (int i = 0; i < num_mers; i++){
      string_hashes[i] = getHash(getCanonical(g.substring(i, i+k)), hashType);
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
}
