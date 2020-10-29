/*
* Generate the screen for input genomes
*/

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

public class ScreenGenerator {
  // User-specified Variables
  public int targetMatches;
  public int readLen;
  public double readErr;
  public String genomeFolder;
  public int k;
  public int numGenomes;
  public int[] sketch_size;
  public String[] genomeNames;
  public ArrayList<HashSet<String>> sketch;
  public ArrayList<HashSet<Integer>> sketch_hash;
  public String screenType;
  public int window;

  public static void main(String[] args) throws Exception
  {

    // Constructor

  }

  // ----- UTILITY FUNCTIONS ------
  // Compute the sketch size, given a genome, read lengths, error rates and # of target matches
  int getSketchSize(int genomeLength, int readLength, double readError, int targetMatches, int k){
    double unaffectedChance = 1.0 - readError;
    // System.out.println(unaffectedChance);
    double multiplier = Math.pow(unaffectedChance, k);
    // System.out.println(multiplier);

    double genomeLengthD = (double) genomeLength;
    double readLengthD = (double) readLength;
    double targetMatchesD = (double) targetMatches;

    double numerator = (targetMatchesD * genomeLengthD);
    // System.out.println(numerator);
    double denominator = (readLengthD * multiplier);
    // System.out.println(denominator);

    double num_sketches = numerator/denominator;

    return (int) num_sketches;
  }

  // Generates the reverse complement of a DNA sequence
  String reverseComplement(String sequence){
    String reversed_tmp = sequence.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").toUpperCase();
    String reversed = new StringBuffer(reversed_tmp).reverse().toString();
    return reversed;
  }

  int getHash(String seq){
    return seq.hashCode();
  }

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

  // Loads all .fasta files in the genomes directory
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

  // MINHASH HELPER FUNCTIONS

  HashSet<Integer> getMinHashes(String g, int sketch_size, int k){

    int gl = g.length();
    int boundary = gl - k;
    HashSet<Integer> hashmers_set = new HashSet<Integer>();


    for (int p = 0; p < boundary; p++) {
      int start = p;
      int end = p + k;
      String curr = g.substring(start, end);
      hashmers_set.add(getHash(getCanonical(curr)));
    }
    ArrayList<Integer> hashmers = new ArrayList<Integer>(hashmers_set);
    Collections.sort(hashmers);

    HashSet<Integer> minhashvals = new HashSet<Integer>();

    for (int q = 0; q < sketch_size; q++) {
      minhashvals.add(hashmers.get(q));
    }

    return minhashvals;
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


  // ----- MINIMIZER TRACKING HELPER FUNCTIONS ------
  // Helper Function to save the minimizers from the genome
  void saveMinimizers(int order, ArrayList<Integer> minimizer_list) throws Exception
  {
    String filename = "./tmp/ZYMO/calc/minimizers_hash_double10_" + order + ".txt";

    // String minimizer_string = String.join(",", minimizer_list);
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i <= minimizer_list.size() - 1; i++)
    {
      int num = minimizer_list.get(i);
      sb.append(num);
      sb.append(",");
    }
    String minimizer_string = sb.toString();

    PrintWriter out = new PrintWriter(new File(filename));
    out.println(minimizer_string);
    out.close();
  }
}
